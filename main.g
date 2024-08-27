ll:=3;
mm:=3;
nn:=4;;
LoadPackage("LINS");;
f := FreeGroup( "x", "y", "z" );;
g := f / [ f.1^ll, f.2^mm, f.3^nn, f.1*f.2*f.3 ];;
x := g.1;; y := g.2;; z:= g.3;;
L := LowIndexNormalSubs(g, 100);;


#finds parents of the subgroup L[i] in the list L
parents:=function(L,i)
    local pars,j;
    pars:=[];
    for j in [1..i-1] do
        if IsSubgroup(L[j],L[i]) then
            Add(pars,j);
        fi;
    od;
    return pars;
end;

#checks if H is torsion free. Assumes it is normal. 
IsTF:=function(H,x,y,z,ll,mm,nn)
    local i;
    for i in [1..ll-1] do
        if x^i in H then return false; fi;
    od;
    for i in [1..mm-1] do
        if y^i in H then return false; fi;
    od;
    for i in [1..nn-1] do
        if z^i in H then return false; fi;
    od;
    return true;
end;

#tfL is the torsion-free subset of L
tfL:=[];
for H in L do
    if IsTF(H,x,y,z,ll,mm,nn) then Add(tfL,H); fi;
od;


for i in [1..Length(tfL)] do
    Print(i," ", Index(g,tfL[i]), " ",parents(tfL,i) ,"\n");
od;