F := Rationals();
R<x> := PolynomialRing(F);

num := 256*(x^6 + 3*x^4 + 3*x^2 + 1);
den := x^4;
f := num/den;

print "Belyi map and factorization:";
print f;
print Factorization(Numerator(f - 1728));
print Factorization(Numerator(f));
print Factorization(Denominator(f));

F := Rationals();
R<q> := PuiseuxSeriesRing(F);
jq := jInvariant(q);
S<x> := PolynomialRing(R);

num := 256*(x^6 + 3*x^4 + 3*x^2 + 1);
den := x^4;

print "All roots:";
dif := num - den*jq;
rts := Roots(dif, R);
print rts;

print "x and y on intersection:";
r := rts[1][1];
rp := -r*(r^2 + 1/4);
s := Sqrt(rp);
print r;
print s;

R<x> := PolynomialRing(Rationals());
E := EllipticCurve(x*(x^2 + (1/4)));
print "Intersection:";
print E;
print "j-invariants and conductor of intersection:";
print jInvariant(E);
print Conductor(E);

isog := TwoIsogeny(E ! [0, 0]);
Ep := Codomain(isog);
print "Curve X (Gamma):";
print Ep;
print "j-invariants and conductor of X (Gamma):";
print jInvariant(Ep);
print Conductor(Ep);

r := -r;
rp := (r^2 + 1/4) / r;
sp := (r^2*s - 1/4*s) / r^2;

rp := 16*rp;
sp := -64*sp;
print "x and y pushed forward and scaled:";
print rp;
print sp;
print "Check 0:";
print sp^2 - (rp^3 - 16^2*rp);
