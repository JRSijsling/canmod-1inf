load "order1.m";

B0 := BM;
B1 :=
[
Matrix(QQ, [
    [1, 0],
    [0, 0]]),
Matrix(QQ, [
    [0, 1],
    [0, 0]]),
Matrix(QQ, [
    [0, 0],
    [1, 0]]),
Matrix(QQ, [
    [0, 0],
    [0, 1]])
];
B5 :=
[
Matrix(QQ, [
    [1, 0],
    [0, 0]]),
Matrix(QQ, [
    [0, 1],
    [0, 0]]),
Matrix(QQ, [
    [0, 0],
    [5, 0]]),
Matrix(QQ, [
    [0, 0],
    [0, 1]])
];

/*
Ms := [ Matrix(QQ, [ [1, 0], [0, 1]]) ];
Bs := [ B1 ];
repeat
    MLs := [ Matrix([ Eltseq(b) : b in B ]) : B in Bs ];
    Ls := [ Lattice(ML) : ML in MLs ];
    wnew := RandomGL2(3);
    Bnew := [ wnew*b1*wnew^(-1) : b1 in B1 ];
    MLnew := Matrix([ Eltseq(b) : b in Bnew ]);
    Lnew := Lattice(MLnew);
    Lint := &meet(Ls) meet Lnew;
    inds := [ Index(L, Lint) : L in Ls ];
    if &and[ ind ne 1 : ind in inds ] then
        if &and[ IsInLattice(Eltseq(b), Lint) : b in B0 ] then
            //print inds;
            Append(~Ms, wnew);
            Append(~Bs, Bnew);
            //print Lint;
        end if;
    end if;
until #Ms eq 2;
//until false;
*/

/*
M0 := Matrix([ Eltseq(b) : b in B0 ]);
L0 := Lattice(M0);
repeat
    L2 := SmallOrder(2, 2^4);
until (L2 meet L0) eq L0;
print L2;
*/

/*
repeat
    done := false;
    wnew := RandomGL2(3);
    Bnew := [ wnew*b5*wnew^(-1) : b5 in B5 ];
    MLnew := Matrix([ Eltseq(b) : b in Bnew ]);
    Lnew := Lattice(MLnew);
    if Lnew eq Lint then
        if Lnew meet L2 eq Lint meet L2 then
            print wnew;
            print Lnew meet L2;
            done := true;
        end if;
    end if;
until done;
*/

M2 := Matrix(QQ, [
[ 1,  0,  0,  1],
[ 1,  0,  0, -1],
[ 0,  1,  1,  1],
[ 0,  1, -1,  1]
]);
M5 := Matrix(QQ, [
[ 1,  0,  0,  0],
[ 0,  1,  0,  0],
[ 0,  0,  5,  0],
[ 0,  0,  0,  1]
]);

L2 := Lattice(M2); L5 := Lattice(M5);
Ls := L2 meet L5;
print "Lattice at 2:";
print L2;
print "Lattice at 5:";
print L5;
print "Intersection:";
print Ls;

print "Level:";
prod := 1;
for M in [ M2, M5 ] do
  Ns := [ Matrix(Rationals(), [ [ KroneckerDelta(i, j) : j in [1..4] ] ])*M^(-1) : i in [1..4] ];
  prod *:= LCM(&cat[ [ Denominator(c) : c in Eltseq(N) ] : N in Ns ]);
end for;
print prod;

print "Checking conjugacy with original lattice...";
repeat
    done := false;
    wnew := RandomGL2(2);
    Bnew := [ wnew*b*wnew^(-1) : b in B0 ];
    MLnew := Matrix([ Eltseq(b) : b in Bnew ]);
    Lnew := Lattice(MLnew);
    if Lnew eq Ls then
        done := true;
    end if;
until done;
print "done";

Bs := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(Ls) ];
B2 := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(L2) ];
B5 := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(L5) ];

wnew := Matrix(Rationals(), [ [0, -1], [5, 0] ]);
Bsnew := [ wnew*b*wnew^(-1) : b in Bs ];
MLsnew := Matrix([ Eltseq(b) : b in Bsnew ]);
Lsnew := Lattice(MLsnew);
B2new := [ wnew*b*wnew^(-1) : b in B2 ];
ML2new := Matrix([ Eltseq(b) : b in B2new ]);
L2new := Lattice(ML2new);
B5new := [ wnew*b*wnew^(-1) : b in B5 ];
ML5new := Matrix([ Eltseq(b) : b in B5new ]);
L5new := Lattice(ML5new);

print "Involution and determinant:";
print wnew;
print Determinant(wnew);
print "Check true:";
print Lsnew eq Ls;
//print L2new eq L2;
//print L5new eq L5;

/*
// This does not work:
L4 := Lattice(Matrix(QQ, [
[ 1,  0,  0,  0],
[ 0,  1,  0,  0],
[ 0,  0,  4,  0],
[ 0,  0,  0,  1]
]));
repeat
    done := false;
    wnew := RandomGL2(12);
    B2new := [ wnew*b*wnew^(-1) : b in B2 ];
    ML2new := Matrix([ Eltseq(b) : b in B2new ]);
    L2new := Lattice(ML2new);
    test := (L2new meet L4) eq L4;
    if test then
        print wnew;
        done := true;
    end if;
until done;
*/

/*
// This does not work:
L20 := Lattice(Matrix(QQ, [
[ 1,  0,  0,  0],
[ 0,  1,  0,  0],
[ 0,  0, 20,  0],
[ 0,  0,  0,  1]
]));
repeat
    done := false;
    wnew := RandomGL2(24);
    B20new := [ wnew*b*wnew^(-1) : b in B20 ];
    ML20new := Matrix([ Eltseq(b) : b in B20new ]);
    L20new := Lattice(ML20new);
    test := (L20new meet Ls) eq L20new;
    if test then
        print wnew;
        done := true;
    end if;
until done;
*/
