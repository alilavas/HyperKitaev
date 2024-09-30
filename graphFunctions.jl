"""
    adjMatX(ct)

returns the adjacancy matrix corresponding to ct, with elements 1,2,3 corresponding to
red,green and blue edges resepectively
"""
function adjMatX(ct)
    n=2*maximum(ct);
    adjMat = zeros(Int, n, n)
    for i in 1:(n÷2)
        adjMat[2*i, 2*i-1] = 1
        adjMat[2*i, 2*ct[6,i]-1] = 2
        adjMat[2*i, 2*ct[3,i]-1] = 3
    end
    adjMat += adjMat'
    return adjMat
end

"""
    redFace(i, ct)

return the red face next to vertex i, as an ordered list of vertices, going anti-clockwise

#Arguments:
-'ct::Matrix{Int}': the coset table
"""
function redFace(i, ct)
    face = Int[]
    push!(face, 2i)
    push!(face, 2ct[6,i]-1)
    j = ct[4,ct[6,i]]
    while j != i
        push!(face, 2j)
        push!(face, 2ct[6,j]-1)
        j = ct[4,ct[6,j]]
    end
    return face
end

"""
    greenFace(i, ct)

return the green face next to vertex i, as an ordered list of vertices, going anti-clockwise.

#Arguments:
-'ct::Matrix{Int}': the coset table
"""
function greenFace(i, ct)

    face = Int[]
    push!(face, 2i-1)
    push!(face, 2i)
    j = ct[3,i]
    while j != i
        push!(face, 2j-1)
        push!(face, 2j)
        j = ct[3,j]
    end
    return face
end


"""
    blueFace(i, ct)

return the blue face next to vertex i, as an ordered list of vertices, going anti-clockwise.

#Arguments:
-'ct::Matrix{Int}': the coset table
"""
function blueFace(i, ct)
    face = Int[]
    push!(face, 2i)
    push!(face, 2i-1)
    j = ct[5,i]
    while j != i
        push!(face, 2j)
        push!(face, 2j-1)
        j = ct[5,j]
    end
    return face
end



"""
    newFace(F, f)

check whether the face f is already in faces F
"""
function newFace(F, f)
    for ff in F
        if length(findall(in(ff),f)) > 2
            return false
        end
    end
    return true
end


"""
    facesHaveRightIntersection(F)

check if the set of faces F has the correct intersection with eachother
"""
function facesHaveRightIntersection(F)
    for i in eachindex(F)
        for j in 1:(i-1)
            intersect_length = length(findall(in(F[i]),F[j]))
            if !(intersect_length == 0 || intersect_length == 2)
                return false
            end
        end
    end
    return true
end



"""
    faces(ct)

return the list of faces, each face is represented by a list of
vertices around the face, orderd as moving anticlockwise around 
the face

#Arguments:
-'ct::Matrix{Int}': the coset table
"""
function faces(ct)
    n=2*maximum(ct)
    F = Vector{Vector{Int64}}()
    for i in 1:n÷2
        fs = [redFace(i, ct), greenFace(i, ct), blueFace(i, ct)]
        for f in fs
            if newFace(F, f)
                push!(F, f)
            end
        end
    end

    return F
    
end


"""
    edges(ct)

return the list of edges, each edges is represented by a 2 element
list holding the endpoint vertices

#Arguments:
-'ct::Matrix{Int}': the coset table
"""
function edges(ct)
    Edges = Vector{Vector{Int}}()
    n=2*maximum(ct)
    for i in 1:(n÷2)
        push!(Edges, sort([2i, 2i-1]))
        push!(Edges, sort([2i, 2ct[6,i]-1]))
        push!(Edges, sort([2i, 2ct[3,i]-1]))
    end
    return Edges
end



"""
    genus(ct)

computes the genus of the surface corresponding to ct

#Arguments:
-'ct::Matrix{Int}': the coset table
"""
function genus(ct)

    V=2*maximum(ct)
    E=length(edges(ct))
    F=length(faces(ct))
    dg=(2 - V + E - F)
    if dg % 2 != 0
        println("ERROR in genus calculation: 2g is not even!")
        return -1
    end
    return dg÷2
end


"""
    faceToEdges(f)

given an ordered list of vertices in a face f, returns the corresponding set of edges 
"""
function faceToEdges(f)
    k=length(f)
    edges=[sort([f[i], f[i%k+1]]) for i in 1:k]
    return edges;
end


"""
    dEdges(d2)

return the edges in the dual lattice, takes the boundary map from faces to edges as input
"""
function dEdges(d2)
    d2T=transpose(d2)
    dE=[Vector{Int}(undef,2) for e in axes(d2T, 2)]
    for e in axes(d2T, 2)
        de=sort(findall(isone,d2T[:,e]))
        dE[e]=de
    end
    return dE 
    
end


"""
        geometry(ct)

return the edges E, the dual edges dE, the faces F, the boundary map 
from edges to vertices d1, the boundary map from faces to edges d2, 
the graph g and the dual graph dg

"""
function geometry(ct)


    n = 2*maximum(ct)

    F = faces(ct)

    E = edges(ct)

    d2_mat=_d2(E,F)

    d1_mat=_d1(E)
    
    g = _graph(E)

    dE = dEdges(d2_mat)
    dg = _graph(dE)

    return E,dE, F,d1_mat,d2_mat,g,dg
end



"""
    _graph(E)

return the graph g given the set of edges E. g[u] will be the vertices neighboring u.

# Arguments:
- 'E::Vector{Vector{Int}}': Array of edges e, where each e=[u,v] represent an edge
                            between vertex u and v. Assumes e=sort(e)
"""
function _graph(E)
    n=maximum(maximum.(E))
    g = [Vector{Int}(undef,0) for i=1:n]
    for e in E
        push!(g[e[1]],e[2])
        push!(g[e[2]],e[1])
    end
     return g
end




"""
    connectedComponent(g,v)
    
compute the connected subgraph of the graph g connected to vertext v
"""
function connectedComponent(g,v)
    visited=zeros(UInt8,length(g))
    _visitNeighbors!(g,v,visited)
    return findall(isone,visited)
end



"""
    _visitNeighbors!(g,v,visited)

visit vertex v (meaning update the corresponding element in visited to 1), and then visit its neighbors on the graph g.
it is a helper function to implement connectedComponent function  recursively 
"""
function _visitNeighbors!(g,v,visited)
    if visited[v]==1
        return 
    end

    visited[v]=1
    for u in g[v]
        _visitNeighbors!(g,u,visited)
    end
    return
end



"""
    removeEdge!(g,e)
    
return the graph g with the edge e=[u,v] removed.
does nothing if there is no edge between v and u in g.

"""
function removeEdge!(g,e)
    u,v=e
    setdiff!(g[v],[u])
    setdiff!(g[u],[v])
    return 
end


"""
    _d2(E,F)
    
return the d2 boundary operator from faces F to edges E.

# Arguments
- 'E::Vector{Vector{Int}}' : Array of edges e, where each e=[u,v] represent an edge between vertex u and v. Assumes e=sort(e)
- 'F::Vector{Vector{Int}}' : Array of faces f, where each f=[u,v,w,...] represent the vertices around face f, anticlockwise ordered.
"""
function _d2(E,F)
    d2_mat=zeros(Int, length(E), length(F))
    for j in eachindex(F)
        f=F[j]
        for e in faceToEdges(f)
            i=findfirst(isequal(e),E)
            d2_mat[i,j]=1
        end
    end
    return d2_mat
end



"""
    _d1(E)
return the d1 boundary operator from edges E to the vertice.

# Arguments
- 'E::Vector{Vector{Int}}' : Array of edges e, where each e=[u,v] represent an edge between vertex u and v. Assumes e=sort(e)
"""
function _d1(E)

    n = maximum(maximum.(E))
    d1_mat=zeros(Int, n, length(E))
    for j in eachindex(E)
        d1_mat[E[j][1],j]=1
        d1_mat[E[j][2],j]=1
    end

    return d1_mat
end

"""
    checkBoundaryOps(d1,d2)

test whether d1 and d2 are valid boundary operators. Basically check d1*d2=0 mod 2
"""
function checkBoundaryOps(d1,d2)
    if sum(mod.(d1*d2,2))>0
        println("ERORR: boundary of boundary is not empty")
        return false
    end
    return true
end


"""
    shortestPaths(graph, source)

compute shortest paths to the source. returns 2 arrays, dist and prev. dist[u] is the shortest path from 
source to u. prev[u] is the previous vertex on the shortest path from source to u, which can be used to 
find the shortest path from the source to any vertex. 
"""
function shortestPaths(graph, source)
    n=length(graph)
    dist=[Inf for i=1:n]
    dist[source]=0
    prev=[Int[] for i=1:n]
    queue=[i for i=1:n]
    while length(queue)>0
        u=argmin(u->dist[u],queue)
        for v in graph[u]
            if dist[u]+1 < dist[v]
                dist[v]=dist[u]+1
                prev[v]=[u]
            elseif dist[u]+1 == dist[v]
                append!(prev[v],u)
            end
        end
        filter!(!isequal(u),queue)
    end

    return dist, prev
end


"""
    shortestNonContractableCycle(E,dE,g,dg,d1,d2, dist, prev)

find a shortest noncontractible cycle
"""
function shortestNonContractableCycle(E,dE,g,dg,d1,d2, dist, prev)
    n=length(g)
    L=Inf
    systol=Array{Int}(undef,0)
    for u in 1:n
        l=dist[u]
        for v in g[u]
            if (dist[v]==l) || (dist[v]==l-1 && v!=prev[u][1]) 
                c_seq= makeVertexSeq( prev, u, v)
                c=vertexSequenceToEdges(E,c_seq)
                if !isContractible(d1,d2,dg,dE,c)
                    if length(c)<L
                        systol=c
                        L=length(c)
                    end
                end
            end
        end
    end
    return systol
end




"""
    _shortestPathFromTree(prev,u)

find the shortest path from u to the source, given the prev graph structure

#Arguments
- 'prev::Vector{Vector{Int}}':  Array of arrays, where prev[u] is the parents of vertex u found during
                                the Dijkstra algorithem. The function only uses prev[u][1] members which
                                can be thought of to represent a minimal spanning tree 
"""
function _shortestPathFromTree(prev,u)
    path=[u]
    while prev[path[end]]!=[]
        append!(path,prev[path[end]][1])
    end
    return path
end


"""
    shortestPath(graph, source, target)

find the shortest path from source to the target on the given graph. 
returns an array of vertices starting with target and ending with source
"""
function shortestPath(graph, source, target)
    dist,prev=shortestPaths(graph,source);
    return _shortestPathFromTree(prev,target)
end

"""
    makeVertexSeq( prev, u, v)

return an array of vertices representing a loop composed of a shortest path from the 
source (root of prev) to u, from u to v, and then a shortest path from v to source. 
Assumes u and v are neighbors
"""
function makeVertexSeq( prev, u, v)
    gamma1=reverse!(_shortestPathFromTree(prev,u))
    gamma2=_shortestPathFromTree(prev,v)
    return [gamma1;gamma2]
end


"""
    vertexSequenceToEdges(E,seq)

turn a sequence of vertices, given by the array seq, into the array of corresponding edge indices.
"""
function vertexSequenceToEdges(E,seq)
    es=zeros(Int,length(seq)-1)
    for i=1:length(seq)-1
        es[i]=findfirst(isequal(sort([seq[i],seq[i+1]])),E)
    end
    return es
end


"""
    isContractible(d1,d2,dg,dE,cycle)

check if the cycle (array of edge indices)  is contractible or not, given the boundary maps d1 and d2 and
the dual graph dg and the dual edge array dE. cycle should be an array of edge indices, which 
form a simple closed loop (The function does not check whether it is true or not).
"""
function isContractible(d1,d2,dg,dE,cycle)
    parts=_parts(dg,dE,cycle)
    if length(parts)==1
        return false
    else
        for p in parts
            V1,E1=subgraph(d1,d2,p)
            if eulerChar(length(V1),length(E1),length(p))==1
                return true
            end
        end
        return false
    end
end


"""
    _parts(dg,dE,cycle)

return connected components in the dual graph after removing the dual edges 
in cycle (array of edge indices). 
"""
function _parts(dg,dE,cycle)
    dg_c=deepcopy(dg)
    for i in cycle
        removeEdge!(dg_c,dE[i])
    end
    f1,_=dE[cycle[1]]
    cc=connectedComponent(dg_c,f1)

    if length(cc)<length(dg)
        f2=findfirst(x-> !(x in cc), [i for i =1:length(dg)] )
        cc2=connectedComponent(dg_c,f2)


        @assert intersect(cc,cc2)==[]
        @assert length(union(cc,cc2))==length(dg)

        return [cc,cc2]
    else
        return  [cc]
    end
end



"""
    subgraph(d1,d2,F)

contruct the set of vertices and edges belonging to faces F based on boundary maps d1 and d2
"""
function subgraph(d1,d2,F)
    E=Array{Int}(undef,0)
    V=Array{Int}(undef,0)
    for f in F
        append!(E,findall(isone,d2[:,f]))
    end
    E=Set(E)
    for e in E
        append!(V,findall(isone,d1[:,e]))
    end
    V=Set(V)
    return V,E
end



"""
    eulerChar(v,e,f)

calculate the euler characteristics given the number of vertices v, number of edges e and the number of faces f
"""
function eulerChar(v,e,f)
    return v-e+f
end


"""
    compactifyEdges!(F,es)

return a new set of faces where each edge e in es is compactified to a point,
i.e. its endpoints are identified. 
"""
function compactifyEdges!(F,es)
    n=maximum(maximum.(F))
    offset=zeros(Int,n)

    for e in es
        u,v=sort(e)
        for f in F
            if u in f
                filter!(!isequal(v),f)
            elseif v in f
                v_idx=findfirst(isequal(v),f)
                f[v_idx]=u
            end
        end
        offset[v+1:end].+=1
    end
    for f in F
        for i in eachindex(f)
            f[i]=f[i]-offset[f[i]]
        end
    end
    return F
end


"""
    edgesFromFaces(F)

construct the set of edges for the set of faces F.
"""
function edgesFromFaces(F)
    E=Vector{Int}[]
    for f in F
        for i=1:length(f)-1
            e=sort([f[i],f[i+1]])
            if !(e in E)
                push!(E,e)                
            end
        end
        e=sort([f[end],f[1]])
        if !(e in E)
            push!(E,e)
        end
    end
    return E
end


"""
    geometryFromFace(F)

similar to geometry function, but it takes teh set of faces F as the input.
"""
function geometryFromFace(F)

    E = edgesFromFaces(F)

    d2_mat=_d2(E,F)

    d1_mat=_d1(E)
    
    g = _graph(E)

    dE = dEdges(d2_mat)
    dg = _graph(dE)

    return E,dE,d1_mat,d2_mat,g,dg
end


"""
    coloredEdges(E,coloredAdjMat,color)

find the subset of edges in E which have color color

#Arguments:
-'coloredAdjMat::Matrix{Int}':  The colored adjacancy matrix where [i,j] component is
                                equal to 1,2 or 3 if there is a red, green or blue edge
                                between i and j.
"""
function coloredEdges(E,coloredAdjMat,color)
    cE=Vector{Int}[]
    for e in E
        u,v=e
        if coloredAdjMat[u,v]==color
            push!(cE,e)
        end
    end
    return cE
end