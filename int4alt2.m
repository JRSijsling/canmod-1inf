load "order4alt2.m";

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
    wnew := RandomGL2(12);
    Bnew := [ wnew*b1*wnew^(-1) : b1 in B1 ];
    MLnew := Matrix([ Eltseq(b) : b in Bnew ]);
    Lnew := Lattice(MLnew);
    Lint := &meet(Ls) meet Lnew;
    inds := [ Index(L, Lint) : L in Ls ];
    //if inds[1] eq 9 then
    if &and[ ind ne 1 : ind in inds ] then
        if &and[ IsInLattice(Eltseq(b), Lint) : b in B0 ] then
            print inds;
            Append(~Ms, wnew);
            Append(~Bs, Bnew);
            print Lint;
        end if;
    end if;
    //end if;
until false;
*/

M0 := Matrix([ Eltseq(b) : b in B0 ]);
L0 := Lattice(M0);

M2 := Matrix(QQ, [
[ 1,  0,  0,  1],
[ 1,  0,  0, -1],
[ 0,  1,  1,  1],
[ 0,  1, -1,  1]
]);
/*
repeat
    L3 := SmallOrder(3, 3^4);
until (L3 meet L0) eq L0;
print L3;
*/
M3 := Matrix(QQ, [
[ 1,  0,  0,  1],
[ 0,  1, -1,  0],
[ 1,  0,  0, -2],
[ 0,  1,  2,  0]
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
    wnew := RandomGL2(3);
    Bsnew := [ wnew*b*wnew^(-1) : b in Bs ];
    MLsnew := Matrix([ Eltseq(b) : b in Bsnew ]);
    Lsnew := Lattice(MLsnew);
    test := (Lsnew meet L0) eq L0;
    if test then
        done := true;
    end if;
until done;
print "done";

/*
repeat
    done := false;
    wnew := RandomGL2(12);
    Bsnew := [ wnew*b*wnew^(-1) : b in Bs ];
    MLsnew := Matrix([ Eltseq(b) : b in Bsnew ]);
    Lsnew := Lattice(MLsnew);
    if Determinant(wnew) eq 1 then
    if Lsnew eq Ls then
        tests := [ ];
        Append(~tests, GCD([ Integers() ! c : c in Eltseq(wnew) ]) eq 1);
        Append(~tests, Determinant(wnew) gt 0);
        Append(~tests, not (IsInLattice(Eltseq(wnew), Ls) and (Determinant(wnew) eq 1)));
        // Strict or not
        //Append(~tests, IsScalar(wnew^2));
        Append(~tests, IsInLattice(Eltseq(wnew^2), Ls));
        // Avoid wrong quotient
        //Append(~tests, not (Determinant(wnew) eq 1));
        if &and(tests) then
            print wnew;
            print Determinant(wnew);
            //done := true;
        end if;
    end if;
    end if;
until done;
*/

/* Add two known commutators: */
comm1 := Matrix(Rationals(), [
[ 1,  1],
[ 1,  2]
]);
comm2 := Matrix(Rationals(), [
[ 1, -1],
[-1,  2]
]);

print comm1;
print comm2;
print IsInLattice(Eltseq(comm1), Ls);
print IsInLattice(Eltseq(comm2), Ls);
print IsInLattice(Eltseq(comm2*comm1^(-1)), Ls);
print IsInLattice(Eltseq(comm2^(-1)*comm1), Ls);
print IsInLattice(Eltseq(comm1^2), Ls);
print IsInLattice(Eltseq(comm2^2), Ls);
