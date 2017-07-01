load "ingredients.m";

/* Create generators of Gamma */
R<t> := PolynomialRing(Rationals());
K := SplittingField((t^2 - 2)*(t^2 - 3));
w2 := Roots(t^2 - 2, K)[1][1];
w3 := Roots(t^2 - 3, K)[1][1];
w6 := w2*w3;
lambda := (w2 + w6)/2;

M := Matrix(K, [[1, lambda], [1, lambda^(-1)]]);
w := Matrix(K, [[2*w3, 3*w2]]);
v := w*M^(-1);

alpha := Matrix(K, [[lambda, 0], [0, lambda^(-1)]]);
beta := Matrix(K, [[v[1,1], w2], [w2, v[1,1]]]);
gamma := -alpha*beta*alpha^(-1)*beta^(-1);

/* Verifications (omitted)
print Determinant(alpha);
print Determinant(beta);
print Determinant(gamma);
print Eigenvalues(gamma);
*/

B := [ IdentityMatrix(K, 2), alpha^2, beta^2, alpha^2*beta^2 ];
gens := [ alpha^2, beta^2, beta*gamma*beta^(-1), beta*alpha*gamma*alpha^(-1)*beta^(-1), gamma, alpha*gamma*alpha^(-1) ];

/* Checking which elements admit a rational expression in the order ZZ [Gamma^(2)]
print [ MatrixInBasis(gen, B) : gen in gens ];
print MatrixInBasis(alpha, B);
print MatrixInBasis(beta, B);
print MatrixInBasis(alpha*beta, B);
*/

/* Create corresponding algebra and trivialize */
seqs := [ ];
for b1 in B do
    seq := [ ];
    for b2 in B do
        Append(~seq, Eltseq(MatrixInBasis(b1*b2, B)));
    end for;
    Append(~seqs, seq);
end for;
A := Algebra< Rationals(), 4 | seqs >;
A := AssociativeAlgebra(A);
test, Q, toQ := IsQuaternionAlgebra(A);
test, M, toM := IsMatrixRing(Q : Isomorphism := true);

OOQ := QuaternionOrder([ toQ(A ! Eltseq(MatrixInBasis(gen, B))) : gen in gens ]);
BQ := Basis(OOQ);
BM := [ toM(b) : b in BQ ];
OOM := Order(Integers(), BM);

print "Basis of original order:";
print BM;
print "Index of original order:";
print Discriminant(OOQ);
print "Eichler?";
print IsEichler(OOQ);

BQs := Enlarge(BQ, 2);
BQ := BQs[1];
OOQ := QuaternionOrder(BQ);

BQs := Enlarge(BQ, 3);
BQ := BQs[1];
OOQ := QuaternionOrder(BQ);

/* Conjugate back into SL_2 (ZZ) and M_2 (ZZ) */
stop := false;
repeat
    D := [-3..3];
    w := Matrix(Rationals(), [[Random(D), Random(D)], [Random(D), Random(D)]]);
    if Determinant(w) ne 0 then
        gensM := [ w*toM(toQ(A ! Eltseq(MatrixInBasis(gen, B))))*w^(-1) : gen in gens ];
        if &and[ &and[ IsIntegral(c) : c in Eltseq(b) ] : b in gensM ] then
            BM := [ w*toM(b)*w^(-1) : b in BQ ];
            stop := true;
        end if;
    end if;
until stop;
OOM := Order(Integers(), BM);

print "Basis of enlarged order:";
print BM;
print "Index of enlarged order:";
print Discriminant(OOQ);
print "Eichler?";
print IsEichler(OOQ);

print "Permutation triple:";
tup := PermutationTriple(OOM);
sigma2, sigma3, sigmainf := Explode(tup);
print tup;
print sigmainf*sigma3*sigma2;

/* Some invariants */
S := Parent(tup[1]);
G := sub<S | tup>;
print "Size of monodromy group:";
print #G;
print "Size of its center:";
print #Center(G);

/* Calculations with subgroup lattice */
H := Stabilizer(G, 1);
Lat := SubgroupLattice(G);
print "Subgroup lattice:";
print Lat;

print "Check that we have group number 51:";
print IsConjugate(G, H, Lat[51]);

print "Permutation representations for intermediate covers:";

/*
f51 := CosetAction(G, Lat[51]);
print "---";
print f51(sigma2); print f51(sigma3); print f51(sigmainf);

f78 := CosetAction(G, Lat[78]);
print "---";
print f78(sigma2); print f78(sigma3); print f78(sigmainf);

f84 := CosetAction(G, Lat[84]);
print "---";
print f84(sigma2); print f84(sigma3); print f84(sigmainf);
*/

f70 := CosetAction(G, Lat[70]);
print "---";
print f70(sigma2); print f70(sigma3); print f70(sigmainf);

f81 := CosetAction(G, Lat[81]);
print "---";
print f81(sigma2); print f81(sigma3); print f81(sigmainf);
