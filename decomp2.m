/* Decomposing a cover encountered in case II */

S := Sym(12);
gens := [ S |
(1, 2)(3, 4)(5, 10)(6, 9)(7, 8)(11, 12),
(1, 3, 6)(2, 5, 4)(7, 9, 10)(8, 11, 12),
(1, 4)(2, 6, 7, 12, 8, 10)(3, 5, 9)
];
sigma2, sigma3, sigmainf := Explode(gens);
print "Permutation triple:";
print gens;

G := sub<S | gens>;
print "Size of monodromy group:";
print #G;
print "Size of its center:";
print #Center(G);

H := Stabilizer(G, 1);
Lat := SubgroupLattice(G);
print "Subgroup lattice:";
print Lat;

print "Check that we have group number 14:";
test, conj := IsConjugate(G, H, Lat[14]);
Hp := H^conj; K := Lat[22];
print "Check that we have can conjugate into the larger group number 22:";
print Hp.1 in K and Hp.2 in K;

print "Coset action for Y --> X:";
YtoX := CosetAction(G, K);
print YtoX(sigma2); print YtoX(sigma3); print YtoX(sigmainf);

/* Determine cosets for the step Y -> X */
lcs := [ G ! 1 ];
repeat
    cands := [ sigma2*lc : lc in lcs ] cat [ sigma3*lc : lc in lcs ] cat [ sigmainf*lc : lc in lcs ];
    stop := true;
    for cand in cands do
        isnew := true;
        for lc in lcs do
            if lc^(-1)*cand in K then
                isnew := false;
                break;
            end if;
        end for;
        if isnew then
            Append(~lcs, cand);
            stop := false;
        end if;
    end for;
until stop;
d := #lcs;
print d;

print "1728, single ramification:";
for lc in lcs do
    print lc;
    print lc^(-1)*sigma2*lc in K;
end for;

print "0, single ramification:";
for lc in lcs do
    print lc;
    print lc^(-1)*sigma3*lc in K;
end for;

print "inf, single ramification:";
for lc in lcs do
    print lc;
    print lc^(-1)*sigmainf*lc in K;
end for;

print "Coset action for Z --> Y:";
ZtoY := CosetAction(K, H);
print ZtoY(sigma2^2); print ZtoY((lcs[4]^(-1)*sigma2*lcs[4])^2);
print ZtoY(sigma3^3); print ZtoY(lcs[4]^(-1)*sigma3*lcs[4]);
print ZtoY(sigmainf); print ZtoY((lcs[4]^(-1)*sigmainf*lcs[4])^3);

/*
Outcome should be:
At the two points of Y over 0 nothing happens and we get 1^3
At the point of Y over 1 with index 3 we also get 1^3
At the point of Y over 1 with index 1 we get 3^1
At the point of Y over inf with index 3 we get 2^1 1^1
At the point of Y over inf with index 1 we get 2^1 1^1
*/
