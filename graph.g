# For a given triangle group g and a torsion free normal subgroup H
# finds the trivalent three-colorable graph associated with it.
# returns the adjacency matrix with element 1,2,3 corresponding to the
# color that is assigned to an edege. 

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


ll:=3;
mm:=3;
nn:=3;
LoadPackage("LINS");;
f := FreeGroup( "x", "y", "z" );;
g := f / [ f.1^ll, f.2^mm, f.3^nn, f.1*f.2*f.3 ];;
x := g.1;; y := g.2;; z:= g.3;;

H:=Subgroup(g,[ (y*x^-1)^2, (y^-1*x)^2 ]);

if not IsSubgroup(g,H) then Print("ERROR: H is not a subgroup"); fi;
if not IsNormal(g,H) then Print("ERROR: H is not a normal"); fi;
if not IsTF(H,x,y,z) then Print("ERROR: H is not torsion free"); fi;

n:=2*Index(g,H);
adjMat:=NullMat(n,n);
ct:=CosetTable(g,H);
for i in [1 .. n/2] do
    adjMat[2*i,2*i-1]:=1;
    adjMat[2*i,2*ct[6,i]-1]:=2;
    adjMat[2*i,2*ct[3,i]-1]:=3;
od;

adjMat:=adjMat+TransposedMat(adjMat);

PrintArray(adjMat);


