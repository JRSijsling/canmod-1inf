load "order2.m";

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

M0 := Matrix([ Eltseq(b) : b in B0 ]);
L0 := Lattice(M0);

M3 := Matrix(QQ, [
[ 1,  0,  0,  0],
[ 0,  1,  0,  0],
[ 0,  0,  3,  0],
[ 0,  0,  0,  1]
]);
/*
repeat
    L2 := SmallOrder(2, 2^6);
until (L2 meet L3) eq L0;
print L2;
*/
M2 := Matrix(QQ, [
[ 1,  0,  0,  1],
[ 1,  0,  0, -1],
[ 0,  0,  2,  0],
[ 0,  2,  1,  0]
]);

L2 := Lattice(M2);
L3 := Lattice(M3);
Ls := L2 meet L3;

print "Lattice at 2:";
print L2;
print "Lattice at 3:";
print L3;
print "Intersection:";
print Ls;

print "Level:";
prod := 1;
for M in [ M2, M3 ] do
  Ns := [ Matrix(Rationals(), [ [ KroneckerDelta(i, j) : j in [1..4] ] ])*M^(-1) : i in [1..4] ];
  prod *:= LCM(&cat[ [ Denominator(c) : c in Eltseq(N) ] : N in Ns ]);
end for;
print prod;

B2 := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(L2) ];
B3 := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(L3) ];
Bs := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(Ls) ];

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

/*
repeat
    done := false;
    wnew := RandomGL2(4);
    Bsnew := [ wnew*b*wnew^(-1) : b in Bs ];
    MLsnew := Matrix([ Eltseq(b) : b in Bsnew ]);
    Lsnew := Lattice(MLsnew);
    if Lsnew eq Ls then
        tests := [ ];
        Append(~tests, GCD([ Integers() ! c : c in Eltseq(wnew) ]) eq 1);
        Append(~tests, Determinant(wnew) gt 0);
        Append(~tests, not (IsInLattice(Eltseq(wnew), Ls) and (Determinant(wnew) eq 1)));
        // Strict or not
        Append(~tests, IsScalar(wnew^2));
        //Append(~tests, IsInLattice(Eltseq(wnew^2), Ls));
        // Avoid wrong quotient
        //Append(~tests, not (Determinant(wnew) eq 1));
        if &and(tests) then
            print wnew;
            print Determinant(wnew);
            //done := true;
        end if;
    end if;
until done;
*/

/* Note: these matrices do not commute, but w2*w3 and w6 differ by an element
 * of the order */
w6 := Matrix(Rationals(), [
[ 0,  2],
[-3,  0]
]);
w2 := Matrix(Rationals(), [
[-2, -2],
[ 3,  2]
]);
w3 := Matrix(Rationals(), [
[-3, -4],
[ 3,  3]
]);

wnew := w2;
print "Involution and determinant:";
print wnew;
print Determinant(wnew);
print "Fixing of sublattices:";
B2new := [ wnew*b*wnew^(-1) : b in B2 ];
ML2new := Matrix([ Eltseq(b) : b in B2new ]);
L2new := Lattice(ML2new);
print L2new eq L2;
B3new := [ wnew*b*wnew^(-1) : b in B3 ];
ML3new := Matrix([ Eltseq(b) : b in B3new ]);
L3new := Lattice(ML3new);
print L3new eq L3;

wnew := w3;
print "Involution and determinant:";
print wnew;
print Determinant(wnew);
print "Fixing of sublattices:";
B2new := [ wnew*b*wnew^(-1) : b in B2 ];
ML2new := Matrix([ Eltseq(b) : b in B2new ]);
L2new := Lattice(ML2new);
print L2new eq L2;
B3new := [ wnew*b*wnew^(-1) : b in B3 ];
ML3new := Matrix([ Eltseq(b) : b in B3new ]);
L3new := Lattice(ML3new);
print L3new eq L3;

wnew := w6;
print "Involution and determinant:";
print wnew;
print Determinant(wnew);
print "Fixing of sublattices:";
B2new := [ wnew*b*wnew^(-1) : b in B2 ];
ML2new := Matrix([ Eltseq(b) : b in B2new ]);
L2new := Lattice(ML2new);
print L2new eq L2;
B3new := [ wnew*b*wnew^(-1) : b in B3 ];
ML3new := Matrix([ Eltseq(b) : b in B3new ]);
L3new := Lattice(ML3new);
print L3new eq L3;

