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


ct:=CosetTable(g,H);


printMatrixToFile:=function(M,filename)
    local r,elm;
    for r in M do
        for elm in r do
            AppendTo(filename,elm," ");
        od;
        AppendTo(filename,"\n");
    od;
end;

printMatrixToFile(ct,"cosetTable");


n:=2*Index(g,H);
adjMat:=NullMat(n,n);
for i in [1 .. n/2] do
    adjMat[2*i,2*i-1]:=1;
    adjMat[2*i,2*ct[6,i]-1]:=2;
    adjMat[2*i,2*ct[3,i]-1]:=3;
od;
adjMat:=adjMat+TransposedMat(adjMat);
# PrintArray(adjMat);



redFace:=function(i,ct)
    local face,j;
    face:=[];
    Add(face,2*i);
    Add(face,2*ct[6,i]-1);
    j:=ct[4,ct[6,i]];
    while not j=i do
        Add(face,2*j);
        Add(face,2*ct[6,j]-1);
        j:=ct[4,ct[6,j]];
    od;
    return face;
end;

greenFace:=function(i,ct)
    local face,j;
    face:=[];
    Add(face,2*i-1);
    Add(face,2*i);
    j:=ct[3,i];
    while not j=i do
        Add(face,2*j-1);
        Add(face,2*j);
        j:=ct[3,j];
    od;
    return face;
end;

blueFace:=function(i,ct)
    local face,j;
    face:=[];
    Add(face,2*i);
    Add(face,2*i-1);
    j:=ct[5,i];
    while not j=i do
        Add(face,2*j);
        Add(face,2*j-1);
        j:=ct[5,j];
    od;
    return face;
end;

newFace:=function(F,f)
    local new, j;
    new:=true;
    for j in [1..Length(F)] do
        if Length(Intersection(f,F[j]))>2 then new:=false; break; fi;
    od;
    return new;
end;


n:=2*Index(g,H);

Faces:=[];;
for i in [1..Index(g,H)] do
    fs:=[redFace(i,ct),greenFace(i,ct),blueFace(i,ct)];
    for f in fs do
        if newFace(Faces,f) then Add(Faces,f); fi;
    od;
od;

Edges:=[];
for i in [1 .. n/2] do
    Add(Edges,[2*i,2*i-1]);
    Add(Edges,[2*i,2*ct[6,i]-1]);
    Add(Edges,[2*i,2*ct[3,i]-1]);
od;
for e in Edges do
    Sort(e);
od;




facesHaveRightIntersection:=function(F)
    #check if any two face have either zero or two points in common
    local i,j;
    for i in [1..Length(F)] do
        for j in [1.. i-1] do
            if not (Length(Intersection(F[i],F[j]))=0 or Length(Intersection(F[i],F[j]))=2) then
                Print(i," ",j);
                return false;
            fi;
        od;
    od;
    return true;
end;
facesHaveRightIntersection(Faces);

genus:=function(V,E,F)
    return (2-V+E-F)/2;
end;

Print("genus:",genus(n,Length(Edges),Length(Faces)));

### F --d2--> E --d1--> V

d2:=NullMat(Length(Edges),Length(Faces));
d1:=NullMat(n, Length(Edges));

faceToEdges:=function(f)
    local i,edges;
    for i in [1..Length(f)-1] do
        Add(edges,[f[i],f[i+1]]);
    od;
    Add(edges,[f[Length(f)],f[1]]);
    for e in edges

end;



