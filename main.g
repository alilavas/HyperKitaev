ll:=4;
mm:=4;
nn:=4;
LoadPackage("LINS");;
f := FreeGroup( "x", "y", "z" );;
g := f / [ f.1^ll, f.2^mm, f.3^nn, f.1*f.2*f.3 ];;
x := g.1;; y := g.2;; z:= g.3;;


#The following fixes a bug when looking for large indices
#LINS_MaxPGenerators should be around the max index we are looking for
MakeReadWriteGlobal("LINS_MaxPGenerators");
LINS_MaxPGenerators := 1000;;

L := LowIndexNormalSubs(g, 672);;


#finds the depth of the subgroup L[i]
depth:=function(pars,i)
    local p,d,dprime;

    if pars[i]=[] then return 0; fi;
    d:=0;
    for p in pars[i] do
        dprime:=depth(pars,p)+1;
        if d<dprime then d:=dprime; fi;
    od;

    return d;
end;

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
IsTF:=function(H,x,y,z)
    local i;
    for i in [1..Order(x)-1] do
        if x^i in H then return false; fi;
    od;
    for i in [1..Order(y)-1] do
        if y^i in H then return false; fi;
    od;
    for i in [1..Order(z)-1] do
        if z^i in H then return false; fi;
    od;
    return true;
end;

#rotate an element cyclically
RotateGen:=function(g,x,y,z)
    if g=x then return y; fi;
    if g=x^-1 then return y^-1; fi;
    if g=y then return z; fi;
    if g=y^-1 then return z^-1; fi;
    if g=z then return x; fi;
    if g=z^-1 then return x^-1; fi;
end;
RotateWord:=function(w,x,y,z)
    local i,cw;
    cw:=x^0;
    for i in [1..Length(w)] do
        cw:=cw*RotateGen(Subword( w, i, i ),x,y,z);
    od;
    return cw;
end;
#check if a subgroup is symmetric under cyclic rotation of its generators
IsCyclicSymmetric:=function(H,x,y,z)
    local i,gens,rotgens;
    gens:=GeneratorsOfGroup(H);
    rotgens:=[];
    for i in [1..Length(gens)] do
        rotgens[i]:=RotateWord(gens[i],x,y,z);
    od;
    return Group(rotgens)=H;
end;

#save the coset table matrix in a file
PrintCosetTableToFile:=function(g,H,filename)
    local M,r,elm;
    M:=CosetTable(g,H);
    PrintTo(filename,"");
    for r in M do
        for elm in r do
            AppendTo(filename,elm," ");
        od;
        AppendTo(filename,"\n");
    od;
end;

# tfL is the torsion-free subset of L
tfL:=[];
for H in L do
    if IsTF(H,x,y,z) then Add(tfL,H); fi;
od;


# pars[i] is the list of the parents of L[i], i.e. j is in pars[i] if L[i] is a subgroup of L[j].
pars:=[];
for i in [1..Length(tfL)] do
    Add(pars,parents(tfL,i));
    Print(i," ", Index(g,tfL[i]), " ",pars[i] ," ", depth(pars,i)+1, " ");
    if IsCyclicSymmetric(tfL[i],x,y,z) then Print("cyclic \n"); else  Print("non-cyclic \n"); fi;
od;


# Use the following to save the coset table 
PrintCosetTableToFile(g,tfL[1],"cosetTable");

# quit;