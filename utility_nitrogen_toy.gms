Set p "parcels" / p1*p10 /;

Parameter r(p) "nitrogen retained by parcel (kg N per yr)"
    / p1 4000, p2 3500, p3 3000, p4 2800, p5 2500,
      p6 2200, p7 2000, p8 1800, p9 1600, p10 1500 /;

Scalar
    landcost "cost to buy one parcel ($ per yr)"     / 40000 /
    ct       "plant treatment cost ($ per kg N)"     / 15 /
    T        "required nitrogen reduction (kg N per yr)" / 40000 /;

Positive Variable
    x(p) "share of parcel p bought"
    Q    "nitrogen treated at plant (kg N per yr)";

Variable TC "total cost ($ per yr)";

Equation obj, target, maxbuy(p);
obj..       TC =e= sum(p, landcost*x(p)) + ct*Q;
target..    sum(p, r(p)*x(p)) + Q =g= T;
maxbuy(p).. x(p) =l= 1;

Model utility / all /;
Solve utility using lp minimizing TC;
Display x.l, Q.l, TC.l;

