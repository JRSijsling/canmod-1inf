function RandomGL2(B);

D := [-B..B];
while true do
    M := Matrix(Rationals(), [ [ Random(D), Random(D) ], [ Random(D), Random (D) ] ]);
    if Determinant(M) ne 0 then
        return M;
    end if;
end while;

end function;


function IsInLattice(v, L);

B := Basis(L);
MB := Matrix(Rationals(), [ Eltseq(c) : c in B ]);
Mv := Matrix(Rationals(), [ Eltseq(v) ]);
w := Matrix(Solution(MB, Mv));
return &and[ IsIntegral(c) : c in Eltseq(w) ];

end function;


function MatrixInBasis(M, Bs);

MBs := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(b) ] : b in Eltseq(B) ] : B in Bs ]);
MM := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(m) ] : m in Eltseq(M) ] ]);
return Matrix(Solution(MBs, MM));

end function;


function IsClosed(B, p, tup);

Bnew := B;
for i in [1..4] do
    if tup[i] ne 0 then
        b := &+[ tup[i]*B[i] : i in [1..4] ]/p;
        Bnew[i] := b;
        break;
    end if;
end for;

for b1 in Bnew do
    for b2 in Bnew do
        prod := Eltseq(MatrixInBasis(b1*b2, Bnew));
        if not &and[ IsIntegral(c) : c in prod ] then
            return false;
        end if;
    end for;
end for;
return true;

end function;


function Enlarge(B, p);

Bs := [ ];
for tup in CartesianPower([0..(p - 1)], 4) do
    if not &and[ c eq 0 : c in tup ] then
        if IsClosed(B, p, tup) then
            b := &+[ tup[i]*B[i] : i in [1..4] ]/p;
            Append(~Bs, Basis(QuaternionOrder(B cat [b])));
        end if;
    end if;
end for;
return Bs;

end function;


function MakesOrder(B);

for b1 in B do
    for b2 in B do
        prod := Eltseq(MatrixInBasis(b1*b2, B));
        if not &and[ IsIntegral(c) : c in prod ] then
            return false;
        end if;
    end for;
end for;
return true;

end function;


function SmallOrder(B, det0);

D := [-B..B];
while true do
    M := Matrix(Rationals(), 4, 4, [ Random(D) : i in [1..16] ]);
    if Determinant(M) ne 0 then
        L := Lattice(M);
        det := Determinant(L);
        if det eq det0 then
            return L;
        end if;
    end if;
end while;

end function;


function PermutationTriple(OOM);

gen2 := Matrix(Integers(), [[0,1],[-1,0]]);
geninf := Matrix(Integers(), [[1,0],[1,1]]);
gen3 := gen2^(-1)*geninf^(-1);

leftcosets := [ IdentityMatrix(Integers(), 2) ];
repeat
    cands := [ gen2*lc : lc in leftcosets ] cat [ gen3*lc : lc in leftcosets ] cat [ geninf*lc : lc in leftcosets ];
    stop := true;
    for cand in cands do
        isnew := true;
        for lc in leftcosets do
            if lc^(-1)*cand in OOM or -lc^(-1)*cand in OOM then
                isnew := false;
                break;
            end if;
        end for;
        if isnew then
            Append(~leftcosets, cand);
            stop := false;
        end if;
    end for;
until stop;
d := #leftcosets;

S := Sym(d);
seq2 := [ ]; seq3 := [ ]; seqinf := [ ];
for i:=1 to d do
    lci := leftcosets[i];
    prod2 := gen2 * lci; prod3 := gen3 * lci; prodinf := geninf * lci;
    for j:=1 to d do
        lcj := leftcosets[j];
        if lcj^(-1)*prod2 in OOM or -lcj^(-1)*prod2 in OOM then
            Append(~seq2, j);
        end if;
        if lcj^(-1)*prod3 in OOM or -lcj^(-1)*prod3 in OOM then
            Append(~seq3, j);
        end if;
        if lcj^(-1)*prodinf in OOM or -lcj^(-1)*prodinf in OOM then
            Append(~seqinf, j);
        end if;
    end for;
end for;

sigma2 := S ! seq2; sigma3 := S ! seq3; sigmainf := S ! seqinf;
return [ sigma2, sigma3, sigmainf ];

end function;
