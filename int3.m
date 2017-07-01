load "order3.m";

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

Ms := Matrix(Rationals(), [
[ 1,  0,  0,  1],
[ 0,  2, -2,  0],
[ 2, -1,  1, -2],
[ 0,  1,  3,  0]
]);
Ls := Lattice(Ms);
Bs := [ Matrix(Rationals(), 2, 2, Eltseq(b)) : b in Basis(Ls) ];

print "Better lattice:";
print Ls;

print "Level:";
prod := 1;
for M in [ Ms ] do
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

print "Finding an involution:";
repeat
    done := false;
    wnew := RandomGL2(3);
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
        Append(~tests, not (Determinant(wnew) eq 1));
        if &and(tests) then
            print wnew;
            print wnew^(-1);
            print Determinant(wnew);
            done := true;
        end if;
    end if;
until done;
/*
[-1 -3]
[ 1  1]
2
*/
