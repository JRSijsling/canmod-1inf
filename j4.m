QQ := Rationals();
R<q> := PuiseuxSeriesRing(QQ, 30);
jq := jInvariant(q);
S<x> := PolynomialRing(R);

d := -27;
num := x^3 + d;
den := 1;

dif := (-1728/d)*num - den*jq;
rts := Roots(dif, R);

r := rts[1][1];
s := Sqrt(r^3 + d);

print "x and y pushed forward:";
print r;
print s;
print "Check 0:";
print s^2 - (r^3 + d);

r := 4*r;
s := 8*s;
d := 64*d;

R<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 + d);

print "Curve X (Gamma):";
print E;
print "j-invariants and conductor of X (Gamma):";
print jInvariant(E);
print Conductor(E);

print "x and y pushed forward and scaled:";
print r; print s;
print "Check 0:";
print s^2 - (r^3 + d);
