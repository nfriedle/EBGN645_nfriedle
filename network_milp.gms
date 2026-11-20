$title Utility nitrogen target: buy upstream pixels or treat at the plant

$include network_data.gms

* ================= SETTINGS: change these =================
Scalar
    ct       "plant treatment cost ($ per kg N)"                   / 45 /
    T        "total N to remove, land + plant (kg N per yr)"       / 45000 /
    w        "plant-equivalent value of 1 kg removed upstream"     / 0.006 /;

$set solvetype mip
* ==========================================================

Binary   Variable buy(i)  "1 if pixel i is bought";
Positive Variable
    Nin(i) "N arriving at pixel i (kg N per yr)"
    R(i)   "N removed at pixel i (kg N per yr)"
    Q      "N removed at the plant (kg N per yr)";
Variable TC "total annual cost ($ per yr)";

Equation obj, inflow(i), cap1(i), cap2(i), target;

obj..       TC =e= sum(i, costA(i)*buy(i)) + ct*Q;
inflow(i).. Nin(i) =e= e(i) + sum(j$down(j,i), Nin(j) - R(j));
cap1(i)..   R(i) =l= s(i)*Nin(i);
cap2(i)..   R(i) =l= s(i)*M(i)*buy(i);
target..  w*sum(o, M(o) - (Nin(o) - R(o))) + Q =g= T;

Model network / all /;
Option optcr = 0.02, reslim = 120;

Solve network using %solvetype% minimizing TC;

Parameter rep(*) "results";
rep('total_cost')        = TC.l;
rep('land_cost')         = sum(i, costA(i)*buy.l(i));
rep('plant_cost')        = ct*Q.l;
rep('N_land_kg')         = sum(o, M(o) - (Nin.l(o) - R.l(o)));
rep('N_plant_kg')        = Q.l;
rep('pixels_bought')     = sum(i, buy.l(i));
rep('max_land_N_kg')     = sum(i, s(i)*M(i));
rep('shadow_usd_per_kg') = target.m;
Display rep;

Parameter bought(i,*) "pixels bought";
bought(i,'share')$(buy.l(i) > 0.01)      = buy.l(i);
bought(i,'removed_kg')$(buy.l(i) > 0.01) = R.l(i);
bought(i,'cost')$(buy.l(i) > 0.01)       = costA(i)*buy.l(i);

File fb / bought.csv /;
put fb 'id,share,removed_kg,cost_usd_yr' /;
loop(i$(buy.l(i) > 0.01),
    put i.tl:0 ',' buy.l(i):0:4 ',' R.l(i):0:3 ',' (costA(i)*buy.l(i)):0:2 /);
putclose fb;