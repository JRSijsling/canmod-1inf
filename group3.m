load "ingredients.m";

L2 := Lattice(Matrix(QQ, [
[ 1,  0,  0,  1],
[ 0,  2, -2,  0],
[ 2, -1,  1, -2],
[ 0,  1,  3,  0]
]));
B2 := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(L2) ];

lists :=
[
    [[ 3, -4],
    [-2,  3]],

    [[-1, -4],
    [ 2,  7]],

    [[ 25,  16],
    [-36, -23]],

    [[  97,   72],
    [-128,  -95]],

    [[ 9,  8],
    [-8, -7]],

    [[ 49,  64],
    [-36, -47]],

    [[-1, -6],
    [ 1,  5]]
];

gensZZ := [ Matrix(ZZ, list) : list in lists ];
G := SL(2, quo< Integers() | 8 >);
gensred := [ G ! list : list in lists ];
H := sub< G | gensred >;

print "Cardinality of SL_2 (ZZ / 8 ZZ) and subgroup:";
print #G;
print #H;
print "Generators of subgroup:";
print Generators(H);

repeat
    done := false;
    Mnew := RandomGL2(3);
    B2new := [ Mnew*b*Mnew^(-1) : b in B2 ];
    ML2new := Matrix([ Eltseq(b) : b in B2new ]);
    L2new := Lattice(ML2new);
    if L2new eq L2 then
        tests := [ ];
        Append(~tests, GCD([ Integers() ! c : c in Eltseq(Mnew) ]) eq 1);
        Append(~tests, Determinant(Mnew) gt 0);
        Append(~tests, not (IsInLattice(Eltseq(Mnew), L2) and (Determinant(Mnew) eq 1)));
        // Strict or not
        Append(~tests, IsScalar(Mnew^2));
        //Append(~tests, IsInLattice(Eltseq(Mnew^2), L2));
        // Avoid wrong quotient
        Append(~tests, not (Determinant(Mnew) eq 1));
        if &and(tests) then
            print "Matrix of order 2 defining involution and determinant:";
            print Mnew;
            print Determinant(Mnew);
            done := true;
        end if;
    end if;
until done;
