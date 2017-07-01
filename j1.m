F := Rationals();
R<x> := PolynomialRing(F);

num := -x^6 + 30*x^5 - 315*x^4 + 1300*x^3 - 1575*x^2 + 750*x - 125;
den := x;
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

num := -x^6 + 30*x^5 - 315*x^4 + 1300*x^3 - 1575*x^2 + 750*x - 125;
den := x;

print "All roots:";
dif := num - den*jq;
rts := Roots(dif, R);
print rts;

print "x and y on intersection:";
r := rts[2][1];
rp := -r*(r^2 - 22*r + 125);
s := Sqrt(rp);
print r;
print s;

R<x> := PolynomialRing(Rationals());
X := HyperellipticCurve(-x*(x^2 - 22*x + 125));
E := EllipticCurve(X);
print "Intersection:";
print E;
print "j-invariants and conductor of intersection:";
print jInvariant(E);
print Conductor(E);

E := EllipticCurve(x*(x^2 + 22*x + 125));
isog := TwoIsogeny(E ! [0, 0]);
Ep := Codomain(isog);
print "Curve X (Gamma):";
print Ep;
print "j-invariants and conductor of X (Gamma):";
print jInvariant(Ep);
print Conductor(Ep);

print "Pushforward isogeny:";
print isog;

print "x and y pushed forward:";
r := -r;
rp := (r^2 + 22*r + 125)/r;
sp := (r^2*s - 125*s) / r^2;
print rp; print sp;
print "Check 0:";
print sp^2 - (rp^3 - 44*rp^2 - 16*rp);
