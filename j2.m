F := Rationals();
R<x> := PolynomialRing(F);

num := 55296*x^9 - 497664*x^8 + 1492992*x^7 - 1575936*x^6 + 497664*x^5 - 746496*x^4 + 41472*x^3 - 124416*x^2 - 6912;
den := x^12 - 12*x^11 + 54*x^10 - 96*x^9 - 27*x^8 + 324*x^7 - 276*x^6 - 288*x^5 + 432*x^4 + 64*x^3 - 192*x^2;
//print num/6912;
f := num/den;

print "Belyi map and factorization:";
print f;
print Factorization(Numerator(f - 1728));
print Factorization(Numerator(f));
print Factorization(Denominator(f));

QQ := Rationals();
R<q> := PuiseuxSeriesRing(QQ);
jq := jInvariant(q);
S<x> := PolynomialRing(R);

num := 55296*x^9 - 497664*x^8 + 1492992*x^7 - 1575936*x^6 + 497664*x^5 - 746496*x^4 + 41472*x^3 - 124416*x^2 - 6912;
den := x^12 - 12*x^11 + 54*x^10 - 96*x^9 - 27*x^8 + 324*x^7 - 276*x^6 - 288*x^5 + 432*x^4 + 64*x^3 - 192*x^2;

print "All roots:";
dif := num - den*jq;
rts := Roots(dif, R);
print rts;

print "x and y on intersection:";
r := rts[5][1];
rp := (r - 3)*(r - 2)*r*(r + 1);
sp := Sqrt(-rp);
print r;
print sp;

X := HyperellipticCurve(-x*(x + 1)*(x - 2)*(x - 3));
P0 := X ! [0, 0];
P := X ! [r, sp];
E, XtoE := EllipticCurve(X, P0);
inv := Inverse(XtoE);

/*
print E;
print XtoE;
print jInvariant(E);
print inv(XtoE(P) + XtoE(P));
*/

Q := XtoE(P) + XtoE(P);
rp := Q[1]; sp := Q[2];
rp *:= 4; sp *:= 8;
print "x and y pushed forward:";
print rp; print sp;
print "Check y^2 = x^3 - 4*x^2 - 384*x - 2304:";
print sp^2 - (rp^3 - 4*rp^2 - 384*rp - 2304);
