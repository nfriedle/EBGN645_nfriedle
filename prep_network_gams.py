"""
Prepare a pixel network for the GAMS land-vs-treatment MILP.

Builds the inputs for this model (one node per candidate pixel i):

    Nin(i)  = e(i) + sum(j drains into i, Nout(j))
    R(i)   <= s(i) * Nin(i)
    R(i)   <= s(i) * M(i) * I(i)          I(i) binary: buy pixel i
    Nout(i) = Nin(i) - R(i)
    sum(o, M(o) - Nout(o)) + Q >= T       o = candidates draining to the outlet
    min  sum(i, costA(i) * I(i)) + ct * Q

What it computes for each candidate pixel i:
    M(i)     kg N/yr flowing through i with no purchases (routed n_total_export)
    D(i)     drainage area through i (ha)
    s(i)     removal share if bought = 1 - exp(-c * a / D(i)),
             a = pixel area (ha), c calibrated so a wetland covering
             --design-ratio of its drainage area removes --efficiency
    e(i)     N entering the network at i that does not pass through another
             candidate = M(i) - sum(M(j) for candidates j draining into i)
    cost(i)  land price x pixel area ($), and costA(i) = cost x --annualize
    down     links: candidate j drains (next) into candidate i

Candidates are pixels with drainage area between --min-drainage-ha and
--max-drainage-ha, a valid price and nitrogen flowing through them. If there
are more than --max-candidates, the ones with the lowest stand-alone $ per kg
(costA / (s x M)) are kept.

Inputs:
  --flow-dir   flow_direction.tif from the InVEST NDR workspace (D8, 0-7)
  --export     n_total_export.tif from the InVEST NDR workspace (kg/ha/yr)
  --price      PLACES land value raster, natural log of $ per ha
               (add --price-is-dollars if it is already $ per ha)

Outputs (in --out-dir):
  n_through_kg_yr.tif   kg N/yr flowing through each pixel
  drainage_area_ha.tif  drainage area of each pixel (ha)
  candidates.csv        one row per candidate pixel with all parameters
  links.csv             from_id, to_id (from drains into to)
  candidates.gpkg       candidate points for QGIS
  network_data.gms      GAMS $include file: sets i, j, down, o and parameters
  summary.txt           settings and totals

Example:
  python prep_network_gams.py \
      --flow-dir outputs/intermediate_outputs/flow_direction.tif \
      --export outputs/n_total_export.tif \
      --price places_iowa_utm15n.tif \
      --max-candidates 2000 --out-dir network_run1
"""
import argparse
import csv
import os

import numpy as np
import pygeoprocessing
import pygeoprocessing.routing
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()

HA_PER_M2 = 1e-4
NODATA = -1.0

# pygeoprocessing D8 codes -> (row step, column step)
D8_STEPS = {0: (0, 1), 1: (-1, 1), 2: (-1, 0), 3: (-1, -1),
            4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1)}


def is_nodata(a, nodata):
    bad = ~np.isfinite(a)
    if nodata is not None:
        bad |= np.isclose(a, nodata)
    return bad


def align_to(src, ref_info, dst):
    pygeoprocessing.warp_raster(
        src, ref_info['pixel_size'], dst, 'near',
        target_bb=ref_info['bounding_box'],
        target_projection_wkt=ref_info['projection_wkt'])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--flow-dir', required=True)
    ap.add_argument('--export', required=True)
    ap.add_argument('--price', required=True)
    ap.add_argument('--price-is-dollars', action='store_true')
    ap.add_argument('--efficiency', type=float, default=0.5,
                    help='removal share at the design ratio (calibrates c)')
    ap.add_argument('--design-ratio', type=float, default=0.01,
                    help='wetland area / drainage area giving --efficiency')
    ap.add_argument('--annualize', type=float, default=0.04,
                    help='multiplies one-time land cost to get $ per yr')
    ap.add_argument('--min-drainage-ha', type=float, default=40.0)
    ap.add_argument('--max-drainage-ha', type=float, default=2000.0)
    ap.add_argument('--max-candidates', type=int, default=2000)
    ap.add_argument('--out-dir', default='network_run')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    work = os.path.join(args.out_dir, '_work')
    os.makedirs(work, exist_ok=True)

    fd_info = pygeoprocessing.get_raster_info(args.flow_dir)
    px_w, px_h = fd_info['pixel_size']
    a_ha = abs(px_w * px_h) * HA_PER_M2
    c_rate = -np.log(1.0 - args.efficiency) / args.design_ratio
    log_price = not args.price_is_dollars

    # 1. Align export and price to the flow-direction grid
    export_al = os.path.join(work, 'export_aligned.tif')
    align_to(args.export, fd_info, export_al)
    price_al = os.path.join(work, 'price_aligned.tif')
    align_to(args.price, fd_info, price_al)

    # 2. Export as kg per pixel per year
    exp_nd = pygeoprocessing.get_raster_info(export_al)['nodata'][0]
    weight = os.path.join(work, 'export_kg_per_pixel.tif')

    def to_kg(e):
        out = e.astype(np.float64) * a_ha
        out[is_nodata(e, exp_nd)] = 0.0
        return out

    pygeoprocessing.raster_calculator(
        [(export_al, 1)], to_kg, weight, gdal.GDT_Float64, NODATA)

    total_export = 0.0
    for _, block in pygeoprocessing.iterblocks((weight, 1)):
        total_export += float(block[block > 0].sum())

    # 3. Route: N through each pixel (M) and drainage area (D)
    n_through = os.path.join(args.out_dir, 'n_through_kg_yr.tif')
    pygeoprocessing.routing.flow_accumulation_d8(
        (args.flow_dir, 1), n_through, weight_raster_path_band=(weight, 1))
    acc_px = os.path.join(work, 'flow_accum_pixels.tif')
    pygeoprocessing.routing.flow_accumulation_d8((args.flow_dir, 1), acc_px)
    acc_nd = pygeoprocessing.get_raster_info(acc_px)['nodata'][0]
    drain_ha = os.path.join(args.out_dir, 'drainage_area_ha.tif')

    def to_ha(a):
        out = a.astype(np.float64) * a_ha
        out[is_nodata(a, acc_nd)] = NODATA
        return out

    pygeoprocessing.raster_calculator(
        [(acc_px, 1)], to_ha, drain_ha, gdal.GDT_Float64, NODATA)
    print(f'Routed nitrogen and drainage area ({args.out_dir})')

    # 4. Find eligible pixels and their stand-alone $ per kg
    nt_nd = pygeoprocessing.get_raster_info(n_through)['nodata'][0]
    pr_nd = pygeoprocessing.get_raster_info(price_al)['nodata'][0]
    lo, hi = args.min_drainage_ha, args.max_drainage_ha
    nt_ds, dr_ds, pr_ds = gdal.Open(n_through), gdal.Open(drain_ha), gdal.Open(price_al)
    nt_b, dr_b, pr_b = (d.GetRasterBand(1) for d in (nt_ds, dr_ds, pr_ds))

    cols = {k: [] for k in ('row', 'col', 'M', 'D', 'price_ha', 'screen')}
    for off in pygeoprocessing.iterblocks((n_through, 1), offset_only=True):
        xo, yo, w, h = off['xoff'], off['yoff'], off['win_xsize'], off['win_ysize']
        nt = nt_b.ReadAsArray(xo, yo, w, h).astype(np.float64)
        dr = dr_b.ReadAsArray(xo, yo, w, h).astype(np.float64)
        pr = pr_b.ReadAsArray(xo, yo, w, h).astype(np.float64)
        ok = (~is_nodata(nt, nt_nd) & (nt > 0) & (dr >= lo) & (dr <= hi)
              & ~is_nodata(pr, pr_nd) & (pr > 0))
        if not ok.any():
            continue
        r, c = np.nonzero(ok)
        M, D = nt[r, c], dr[r, c]
        p_ha = np.exp(pr[r, c]) if log_price else pr[r, c]
        s = 1.0 - np.exp(-c_rate * a_ha / D)
        screen = p_ha * a_ha * args.annualize / (s * M)
        for k, v in (('row', r + yo), ('col', c + xo), ('M', M), ('D', D),
                     ('price_ha', p_ha), ('screen', screen)):
            cols[k].append(v)
    nt_ds = dr_ds = pr_ds = nt_b = dr_b = pr_b = None

    if not cols['row']:
        raise SystemExit('No eligible pixels. Check the drainage limits and '
                         'that the price raster covers the watershed.')
    cand = {k: np.concatenate(v) for k, v in cols.items()}
    n_elig = cand['row'].size
    if n_elig > args.max_candidates:
        keep = np.argsort(cand['screen'], kind='stable')[:args.max_candidates]
        cand = {k: v[keep] for k, v in cand.items()}
    n = cand['row'].size
    print(f'{n_elig:,} eligible pixels; keeping {n:,} candidates')

    # 5. Link candidates: walk downstream from each to the next candidate
    fd_ds = gdal.Open(args.flow_dir)
    fd = fd_ds.GetRasterBand(1).ReadAsArray()
    fd_ds = None
    n_rows, n_cols = fd.shape
    idx_of = {int(r) * n_cols + int(c): k
              for k, (r, c) in enumerate(zip(cand['row'], cand['col']))}
    next_cand = {}  # flat pixel index -> next candidate downstream (-1: none)
    down = np.full(n, -1, dtype=np.int64)

    for k in range(n):
        r, c = int(cand['row'][k]), int(cand['col'][k])
        path, result = [], -1
        while True:
            step = D8_STEPS.get(int(fd[r, c]))
            if step is None:
                break
            r, c = r + step[0], c + step[1]
            if not (0 <= r < n_rows and 0 <= c < n_cols):
                break
            flat = r * n_cols + c
            if flat in idx_of:
                result = idx_of[flat]
                break
            if flat in next_cand:
                result = next_cand[flat]
                break
            path.append(flat)
        for flat in path:
            next_cand[flat] = result
        down[k] = result
    del fd, next_cand

    # 6. Parameters: e(i), s(i), cost(i)
    M = cand['M']
    upstream_M = np.zeros(n)
    for k in range(n):
        if down[k] >= 0:
            upstream_M[down[k]] += M[k]
    e = M - upstream_M
    e[e < 0] = 0.0  # guards rounding only; routing makes this non-negative
    s = 1.0 - np.exp(-c_rate * a_ha / cand['D'])
    cost = cand['price_ha'] * a_ha
    cost_a = cost * args.annualize
    ids = [f'c{k + 1}' for k in range(n)]

    gt = gdal.Open(args.flow_dir).GetGeoTransform()
    xs = gt[0] + (cand['col'] + 0.5) * gt[1]
    ys = gt[3] + (cand['row'] + 0.5) * gt[5]

    # 7. Write outputs
    fields = ['id', 'x', 'y', 'drainage_area_ha', 'M_kg_yr', 'e_kg_yr', 's',
              'price_usd_per_ha', 'cost_usd', 'cost_annual_usd_yr',
              'standalone_usd_per_kg', 'drains_into', 'drains_to_outlet']
    rows = []
    for k in range(n):
        rows.append({
            'id': ids[k], 'x': round(float(xs[k]), 1), 'y': round(float(ys[k]), 1),
            'drainage_area_ha': round(float(cand['D'][k]), 3),
            'M_kg_yr': float(f'{M[k]:.6g}'), 'e_kg_yr': float(f'{e[k]:.6g}'),
            's': float(f'{s[k]:.6g}'),
            'price_usd_per_ha': round(float(cand['price_ha'][k]), 2),
            'cost_usd': round(float(cost[k]), 4),
            'cost_annual_usd_yr': round(float(cost_a[k]), 4),
            'standalone_usd_per_kg': float(f"{cand['screen'][k]:.6g}"),
            'drains_into': ids[down[k]] if down[k] >= 0 else '',
            'drains_to_outlet': int(down[k] < 0)})

    with open(os.path.join(args.out_dir, 'candidates.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.out_dir, 'links.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['from_id', 'to_id'])
        for k in range(n):
            if down[k] >= 0:
                w.writerow([ids[k], ids[down[k]]])

    srs = osr.SpatialReference()
    srs.ImportFromWkt(fd_info['projection_wkt'])
    gpkg = os.path.join(args.out_dir, 'candidates.gpkg')
    if os.path.exists(gpkg):
        os.remove(gpkg)
    vds = ogr.GetDriverByName('GPKG').CreateDataSource(gpkg)
    lyr = vds.CreateLayer('candidates', srs, ogr.wkbPoint)
    for fname in fields:
        ftype = ogr.OFTString if fname in ('id', 'drains_into') else (
            ogr.OFTInteger if fname == 'drains_to_outlet' else ogr.OFTReal)
        lyr.CreateField(ogr.FieldDefn(fname, ftype))
    for row in rows:
        feat = ogr.Feature(lyr.GetLayerDefn())
        for fname in fields:
            feat.SetField(fname, row[fname])
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint_2D(row['x'], row['y'])
        feat.SetGeometry(pt)
        lyr.CreateFeature(feat)
    vds = None

    def gms_param(name, desc, vals):
        lines = [f'Parameter {name}(i) "{desc}"', '/']
        lines += [f'  {ids[k]} {vals[k]:.10g}' for k in range(n)]
        lines.append('/;')
        return '\n'.join(lines)

    with open(os.path.join(args.out_dir, 'network_data.gms'), 'w') as f:
        f.write('* Generated by prep_network_gams.py\n')
        f.write('Set i "candidate pixels"\n/\n')
        f.write('\n'.join(f'  {x}' for x in ids))
        f.write('\n/;\nAlias (i, j);\n\n')
        f.write('Set down(i,j) "pixel i drains next into candidate j"\n/\n')
        f.write('\n'.join(f'  {ids[k]}.{ids[down[k]]}'
                          for k in range(n) if down[k] >= 0))
        f.write('\n/;\n\n')
        f.write('Set o(i) "candidates that drain to the outlet with no '
                'candidate below them"\n/\n')
        f.write('\n'.join(f'  {ids[k]}' for k in range(n) if down[k] < 0))
        f.write('\n/;\n\n')
        f.write(gms_param('e', 'N entering the network at i (kg N per yr)', e) + '\n\n')
        f.write(gms_param('M', 'N through i with no purchases (kg N per yr)', M) + '\n\n')
        f.write(gms_param('s', 'removal share if i is bought (unitless)', s) + '\n\n')
        f.write(gms_param('costA', 'annualized land cost of i ($ per yr)', cost_a) + '\n\n')
        f.write(f'Scalar totalExport "total N export, whole raster (kg N per yr)" '
                f'/ {total_export:.10g} /;\n')

    with open(os.path.join(args.out_dir, 'summary.txt'), 'w') as f:
        f.write(f'pixel area (ha): {a_ha}\n')
        f.write(f'efficiency {args.efficiency} at design ratio {args.design_ratio} '
                f'-> c = {c_rate:.4f}\n')
        f.write(f'annualize: {args.annualize}\n')
        f.write(f'drainage range (ha): {lo} - {hi}\n')
        f.write(f'eligible pixels: {n_elig}\ncandidates kept: {n}\n')
        f.write(f'candidates draining to outlet: {int((down < 0).sum())}\n')
        f.write(f'total N export, whole raster (kg/yr): {total_export:.6g}\n')
        f.write(f'max possible abatement if all bought, ignoring stacking '
                f'(kg/yr): {float((s * M).sum()):.6g}\n')
    print(f'Wrote candidates.csv, links.csv, candidates.gpkg, network_data.gms '
          f'and summary.txt to {args.out_dir}')


if __name__ == '__main__':
    main()
