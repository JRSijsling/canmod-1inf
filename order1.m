load "ingredients.m";

/* Create generators of Gamma */
R<t> := PolynomialRing(Rationals());
K<r> := NumberField(t^2 - 5);
lambda := (1 + r)/2;

M := Matrix(K, [[1, lambda], [1, lambda^(-1)]]);
w := Matrix(K, [[2*r, 5]]);
v := w*M^(-1);

alpha := Matrix(K, [[lambda, 0], [0, lambda^(-1)]]);
beta := Matrix(K, [[r, 2], [2, r]]);
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
gens := gens cat [ alpha*beta ];

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

OOQ := QuaternionOrder([ toQ(A ! Eltseq(MatrixInBasis(gen, B))) : gen in Append(gens, alpha*beta) ]);
BQ := Basis(OOQ);
BM := [ toM(b) : b in BQ ];
OOM := Order(Integers(), BM);

print "Basis of original order:";
print BM;
print "Index of original order:";
print Discriminant(OOQ);
print "Eichler?";
print IsEichler(OOQ);

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
