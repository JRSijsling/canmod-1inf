load "order1.m";
/* EXTREMELY rough version */

function MatrixInBasis(M, Bs);

MBs := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(b) ] : b in Eltseq(B) ] : B in Bs ]);
MM := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(m) ] : m in Eltseq(M) ] ]);
return Matrix(Solution(MBs, MM));

end function;


function Act(M, z)

a, b, c, d := Explode(Eltseq(M));
return (a*z + b)/(c*z + d);

end function;


function ReduceToZero(z, genpairs)

H := UpperHalfPlane(); iH := H.1;
dist := Distance(H ! z, iH); hom := Matrix(Rationals(), [[0, 0]]);
while dist gt 0.001 do
    for genpair in genpairs do
        gen, homgen := Explode(genpair);
        znew := Act(gen, z); distnew := Distance(H ! znew, iH);
        if distnew lt dist then
            //print distnew;
            z := znew; dist := distnew; hom +:= homgen;
            break;
        end if;
    end for;
end while;
return hom;

end function;


function ElementOfPrimeNorm(B)

D := [-5..5];
while true do
    x := &+[ Random(D)*b : b in B];
    det := Determinant(x);
    if det gt 5 and det lt 33 and IsPrime(Integers() ! det) then
        return x;
    end if;
end while;

end function;


function InOrder(x, B);

w := MatrixInBasis(x, B);
return &and[ IsIntegral(c) : c in Eltseq(w) ];

end function;


function IsEquivalent(x, y, pi, B)

d := y^-1 * x;
test := InOrder(d, B) and InOrder(pi*d*pi^(-1), B);
return test;

end function;


function CosetsAndConjugates(ki, gens, pi, B)

c := Matrix([[1,0],[0,1]]);
k := [ ];
cosets := [ ];
repeat
    cprod := c;
    repeat
        Append(~cosets, cprod);
        //cprod := cprod * ki;
        //cconj := cprod * c^(-1);
        cprod := ki * cprod;
        cconj := c^(-1) * cprod;
        cdone := InOrder(cconj, B) and InOrder(pi*cconj*pi^(-1), B);
    until cdone;
    Append(~k, pi*cconj*pi^(-1));
    //Append(~k, pi^(-1)*cconj*pi);

    done := true;
    if #cosets lt Determinant(pi) + 1 then
        for gen in gens do
            //ccand := c*gen;
            ccand := gen*c;
            if &and[ not IsEquivalent(ccand, coset, pi, B) : coset in cosets ] then
                c := ccand;
                done := false;
                break;
            end if;
        end for;
    end if;
until done;
return cosets, k;

end function;


gens := gensM;
gens := gens cat [ gen^(-1) : gen in gens ];
gens := gens cat &cat[ [ gen1 * gen2 : gen1 in gens ] : gen2 in gens ] cat &cat[ &cat[ [ gen1 * gen2 * gen3 : gen1 in gens ] : gen2 in gens ] : gen3 in gens ];

homs := [
Matrix(Rationals(), [[2,0]]),
Matrix(Rationals(), [[0,2]]),
Matrix(Rationals(), [[0,0]]),
Matrix(Rationals(), [[0,0]]),
Matrix(Rationals(), [[0,0]]),
Matrix(Rationals(), [[0,0]]),
Matrix(Rationals(), [[1,1]])
];
homs := homs cat [ -hom : hom in homs ];
homs := homs cat &cat[ [ hom1 + hom2 : hom1 in homs ] : hom2 in homs ] cat &cat[ &cat[ [ hom1 + hom2 + hom3 : hom1 in homs ] : hom2 in homs ] : hom3 in homs ];

genpairs := [ [* gens[i], homs[i] *] : i in [1..#gens] ];
CC<iCC> := ComplexField();
/*
print ReduceToZero(iCC, genpairs);
print ReduceToZero(Act(gens[1], iCC), genpairs);
print ReduceToZero(Act(gens[3] * gens[6], iCC), genpairs);
*/

B := BM;
ks := [ gensM[1], gensM[7] ];
Ms := [ homs[1], homs[7] ];

ps := [ 7, 11, 13, 17, 19, 23, 29, 31 ];
ps_done := [ ];

print "Beware of this extremely rickety implementation, it may not even get off the ground...";
while true do
    pi := ElementOfPrimeNorm(B);
    p := Determinant(pi);
    if not p in ps_done then
        Append(~ps_done, p);
        print p;
        rows := [ ];
        for k in ks do
            cosets,prods := CosetsAndConjugates(k, gens, pi, B);
            //print "Number of cosets:";
            //print #cosets;
            row := &+[ -ReduceToZero(Act(prod, iCC), genpairs) : prod in prods ];
            print MatrixInBasis(row, Ms);
        end for;
    end if;
end while;

// 7: 2, 11: 0, 13: 2, 17: -6, 19:-4, 23: 6, 29: 6, 31 : -4
