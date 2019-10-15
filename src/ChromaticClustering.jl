
using SparseArrays
using LinearAlgebra
using StatsBase
using JuMP
using Gurobi

include("SyntheticGenerator.jl")
include("isolating_heuristic.jl")

gurobi_env = Gurobi.Env()

"""
CHROMECLUSOPT
Optimally solves the Chromatic Clustering Objective by solving and ILP.
Does NOT support edges that are colorless.

Input:
    EdgeList = the edge list of a Chromatic Clustering instance
    EdgeColors = the edge colors
    n = number of nodes in the problem

Output:
    OPT = The number of edge color mistakes made
    c   = The assigned color label of nodes
"""
function ChromeClusOpt_Graph(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64)

    k = maximum(EdgeColors)
    M = size(EdgeList,1)

    m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0, gurobi_env))

    # Variables for nodes and edges
    @variable(m, y[1:M])
    @variable(m, x[1:n,1:k],Bin)
    @objective(m, Min, sum(y[i] for i=1:M))
    #
    # @constraint(m,x .<= ones(n,k))
    # @constraint(m,x .>= zeros(n,k))
    # @constraint(m,y .<= ones(M))
    # @constraint(m,y .>= zeros(M))

    for i = 1:n
        @constraint(m, sum(x[i,j] for j = 1:k) == k-1)
    end

    for e = 1:M
        color = EdgeColors[e]
        edge = EdgeList[e,:]

        # For every node in the edge, there's a constraint for the
        # node-color variable
        for v = edge
            @constraint(m, y[e] >= x[v,color])
        end

    end
    start = time()
    JuMP.optimize!(m)
    runtime = time()-start

    # Return clustering and objective value
    X = JuMP.value.(x)
    OPT = JuMP.objective_value(m)

    # An inefficient way to extract the label, which works for now
    X2 = sparse(ones(n,k) - X)
    A,B,C = findnz(X2)
    p = sortperm(A)
    c = B[p]

    return OPT, c, runtime
end




"""
IsolatingCutLP

An LP-based implementation of the isolating cut heuristic.
"""
function IsolatingCutLP(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64,outputflag::Int64=0)
    k = maximum(EdgeColors)
    TwoWayCuts = zeros(n,k)
    M = length(EdgeColors)
    scores = zeros(k)

    # Solve the problem for index i
    for i = 1:k
        inds_i = findall(x->x == i,EdgeColors)  # get the i-colored edges
        non_i = setdiff(1:M,inds_i)             # get the non-i edges
        ecol = ones(Int64,M)                  # get a new coloring vector
        ecol[non_i] .= 2                      # now everything is a 2 if it's not color i

        LPval, X, runtime = ChromeClusLP_Graph(EdgeList,ecol,n,outputflag)
        C = rowmin(X)
        c = C[:,2]
        c_inds = findall(x->x==1,c)
        TwoWayCuts[c_inds,i] .= 1
        scores[i] = LPval
    end

    # Remove the most expensive cut
    a = findmax(scores)[2]
    keepinds = setdiff(1:k,a)

    twc = TwoWayCuts[:,keepinds]
    default = MajorityVote(EdgeList,EdgeColors,n,k)
    cdeg = get_color_degree(EdgeList,EdgeColors,n)
    colors = zeros(n)
    for i = 1:n
        options = keepinds[findall(!iszero,twc[i,:])]

        # if the node is in any of the returned cut sets, choose one at random
        # (in theory, if it is in more than one, then you can put it anywhere
        # and the theoretical approximation guarantee holds)
        if length(options) > 0
            colors[i] = options[rand(1:length(options))]
        elseif cdeg[a,i] > 0
            # otherwise, if the node is attached to the one other terminal,
            # cluster it there
            colors[i] = a
            colors[i] = default[i]
        else
            # in the last case, go with the majority vote default
            colors[i] = default[i]
        end
    end


    return round.(Int64,colors)
end


"""
CHROMECLUSLP
Solves the Chromatic Clustering Linear Programming Relaxation

Input:
    EdgeList = the edge list of a Chromatic Clustering instance
    EdgeColors = the edge colors
    n = number of nodes in the problem

Output:
    LPval = The LP relaxation objective score. This is a lower bound on OPT.
    X  = The distance labels from solving the objective
    runtime = Solver runtime
"""
function ChromeClusLP_Graph(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64,outputflag::Int64=0)

    k = maximum(EdgeColors)
    M = size(EdgeList,1)

    m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=outputflag, gurobi_env))

    # Variables for nodes and edges
    @variable(m, y[1:M])
    @variable(m, x[1:n,1:k])
    @objective(m, Min, sum(y[i] for i=1:M))

    @constraint(m,x .<= ones(n,k))
    @constraint(m,x .>= zeros(n,k))
    @constraint(m,y .<= ones(M))
    @constraint(m,y .>= zeros(M))

    for i = 1:n
        @constraint(m, sum(x[i,j] for j = 1:k) == k-1)
    end

    for e = 1:M
        color = EdgeColors[e]
        edge = EdgeList[e,:]

        # For every node in the edge, there's a constraint for the
        # node-color variable
        for v = edge
            @constraint(m, y[e] >= x[v,color])
        end

    end
    start = time()
    JuMP.optimize!(m)
    runtime = time()-start

    # Return clustering and objective value
    X = JuMP.value.(x)
    LPval= JuMP.objective_value(m)

    return LPval, X, runtime
end


## This can be either the optimal solution or not, and it takes hyperedges
function ChromeClusGeneral(EdgeList::Vector{Vector{Int64}},EdgeColors::Array{Int64,1},n::Int64,optimalflag::Bool= false,outputflag::Int64=0)

    k = maximum(EdgeColors)
    M = length(EdgeList)

    m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=outputflag, gurobi_env))

    # Variables for nodes and edges
    @variable(m, y[1:M])

    @objective(m, Min, sum(y[i] for i=1:M))

    if optimalflag
        @variable(m, x[1:n,1:k],Bin)
    else
        @variable(m, x[1:n,1:k])
        @constraint(m,x .<= ones(n,k))
        @constraint(m,x .>= zeros(n,k))
        @constraint(m,y .<= ones(M))
        @constraint(m,y .>= zeros(M))
    end

    for i = 1:n
        @constraint(m, sum(x[i,j] for j = 1:k) == k-1)
    end

    for e = 1:M
        color = EdgeColors[e]
        edge = EdgeList[e]

        # For every node in the edge, there's a constraint for the
        # node-color variable
        for v = edge
            @constraint(m, y[e] >= x[v,color])
        end

    end
    start = time()
    JuMP.optimize!(m)
    runtime = time()-start

    # Return clustering and objective value
    X = JuMP.value.(x)
    LPval= JuMP.objective_value(m)

    return LPval, X, runtime
end




## This can be either the optimal solution or not, and it takes hyperedges
# and weights if desired
function ChromeClusWeighted(EdgeList::Vector{Vector{Int64}},EdgeColors::Array{Int64,1},EdgeWeights::Array{Float64,1},n::Int64,optimalflag::Bool= false,outputflag::Int64=0)

    k = maximum(EdgeColors)
    M = length(EdgeList)

    m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=outputflag, gurobi_env))

    # Variables for nodes and edges
    @variable(m, y[1:M])

    @objective(m, Min, sum(EdgeWeights[i]*y[i] for i=1:M))

    if optimalflag
        @variable(m, x[1:n,1:k],Bin)
    else
        @variable(m, x[1:n,1:k])
        @constraint(m,x .<= ones(n,k))
        @constraint(m,x .>= zeros(n,k))
        @constraint(m,y .<= ones(M))
        @constraint(m,y .>= zeros(M))
    end

    for i = 1:n
        @constraint(m, sum(x[i,j] for j = 1:k) == k-1)
    end

    for e = 1:M
        color = EdgeColors[e]
        edge = EdgeList[e]

        # For every node in the edge, there's a constraint for the
        # node-color variable
        for v = edge
            @constraint(m, y[e] >= x[v,color])
        end

    end
    start = time()
    JuMP.optimize!(m)
    runtime = time()-start

    # Return clustering and objective value
    X = JuMP.value.(x)
    LPval= JuMP.objective_value(m)

    return LPval, X, runtime
end

"""
CHROMECLUSOBJ

Returns the number of mistakes made by a clustering in an instance of Chromatic
Clustering.

"""
function ChromeClusObj(EdgeList::Union{Array{Int64,2},Vector{Vector{Int64}}},EdgeColors::Array{Int64,1},c::Vector)
    n = length(c); Mistakes = 0
    for i = 1:size(EdgeList,1)
        if size(EdgeList,2) == 2
            edge = EdgeList[i,:]
        else
            edge = EdgeList[i]
        end
        for v in edge
            if c[v] != EdgeColors[i]
                Mistakes += 1
                break
            end
        end
    end
    return Mistakes
end

"""
SIMPLEROUND

A simple deterministic algorithm for rounding the LP relaxation of the Chromatic
Clustering objective. Guaranteed to be a 2-approximation for any problem.
"""
function SimpleRound(EdgeList::Union{Array{Int64,2},Vector{Vector{Int64}}},EdgeColors::Array{Int64,1},X::Array{Float64,2},LPval::Float64)
    C = rowmin(X)
    c = C[:,2]
    RoundScore = ChromeClusObj(EdgeList,EdgeColors,c)
    RoundRatio = RoundScore/LPval
    return c, RoundScore, RoundRatio
end


"""
Julia doesn't seem to have a nice way to extract the entry of the minimum value
in each row of a matrix and return it as an n x 2 matrix. So here it is.
"""
function rowmin(X::Array{Float64,2})
    n = size(X,1); Y = zeros(n,2)
    for i = 1:n
        g = findmin(X[i,:])
        Y[i,1] = g[1]; Y[i,2] = g[2]
    end
    return Y
end

function rowmax(X::Array{Float64,2})
    n = size(X,1)
    Y = zeros(n,2)
    for i = 1:n
        g = findmax(X[i,:])
        Y[i,1] = g[1]
        Y[i,2] = g[2]
    end
    return Y
end


"""
ChromaticBalls

Algorithm from Bonchi et al.
(Chromatic Correlation Clustering, ACM Trans. Knowl. Discov. Data 2015)
"""
function ChromaticBalls(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64)

    # A naive implementation that calls the find function over and over again

    # Convert to a graph first
    A = sparse(EdgeList[:,1],EdgeList[:,2],EdgeColors,n,n)

    if istriu(A)
        A = A+A'
    else
        @assert(issymmetric(A))
    end

    assigned = zeros(n)     # indicator for which nodes have been assigned to a cluster
    # undecided = ones(M)    # indicator for which edges haven't been decided
    node2color = zeros(n)   # We will additionally keep track of the colors of nodes
    cluster = zeros(n)
    clusnum = 1
    K = maximum(EdgeColors)
    # Make a copy, whose entries we'll constantly delete
    Anew = copy(A)

    I,J,EdgeColor = findnz(triu(Anew))
    EdgeList = [I J]
    M = size(EdgeList,1)

    while M > 0

        ind = rand(1:M)
        x = EdgeList[ind,1]
        y = EdgeList[ind,2]
        l = EdgeColors[ind]

        # Neighbors of x of color l
        Nx = findall(x->x == l,Anew[x,:])

        # Neighbors of y of color l
        Ny = findall(x->x == l,Anew[y,:])

        # Neighbors of both of color l
        Nboth = intersect(Nx, Ny)

        # form a new cluster
        newclus = [x;y; Nboth]

        # Given then the cluster label and color
        for v in newclus
            cluster[v] = clusnum
            node2color[v] = l
            assigned[v] = 1
        end
        clusnum += 1

        # Delete all edges involving nodes that were just clustered
        nonclus = setdiff(1:n,newclus)
        Anew[newclus,:] .= 0
        Anew[:,newclus] .= 0
        dropzeros!(Anew)

        # Find the new edge list (yes, yes, this is all very inefficient)
        I,J,EdgeColors = findnz(triu(Anew))
        EdgeList = [I J]
        M = size(EdgeList,1)

    end

    nonlabeled = findall(iszero,assigned)

    for v in nonlabeled
        cluster[v] = clusnum
        clusnum += 1
        # Get all the edges of node v, and assign it a color based on naive-majority
        Vv = nonzeros(A[v,:])
        if length(Vv) > 0
            node2color[v] = StatsBase.mode(Vv)
        else
            node2color[v] = rand(1:K)
        end
    end

    return round.(Int64,cluster), round.(Int64,node2color)
end


"""
LazyChromaticBalls

Another algorithm from Bonchi et al.
(Chromatic Correlation Clustering, ACM Trans. Knowl. Discov. Data 2015)

This code is not at all optimized.
"""
function LazyChromaticBalls(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64)

    # A naive implementation that calls the find function over and over again

    # Convert to a graph first
    A = sparse(EdgeList[:,1],EdgeList[:,2],EdgeColors,n,n)
    # A = A+A'
    if istriu(A)
        A = A+A'
    else
        @assert(issymmetric(A))
    end

    assigned = zeros(n)     # indicator for which nodes have been assigned to a cluster
    # undecided = ones(M)    # indicator for which edges haven't been decided
    node2color = zeros(n)   # We will additionally keep track of the colors of nodes
    cluster = zeros(n)
    clusnum = 1

    K = maximum(EdgeColors)

    # Make a copy, whose entries we'll constantly delete
    Anew = copy(A)

    I,J,EdgeColors = findnz(triu(Anew))
    EdgeList = [I J]
    M = size(EdgeList,1)

    time1 = 0
    time2 = 0

    while M > 0
        # Get the color-degree of each node
        start = time()
        ColorDegree = zeros(n,K)
        for e_ind = 1:M
            I = EdgeList[e_ind,1]
            J = EdgeList[e_ind,2]
            col = EdgeColors[e_ind]
            ColorDegree[I,col] += 1
            ColorDegree[J,col] += 1
        end
        DL = rowmax(ColorDegree)
        Delta = DL[:,1]
        Lambda = DL[:,2]
        time1 += time()-start

        # Sample a pivot node
        X = StatsBase.sample(1:n,aweights(Delta))
        colorx = round(Int64,Lambda[X])

        # Sample another edge y from its neighbors, proportional to the
        # number of edges of colorx that that node y has
        # neighbsX = findall(!iszero,Anew[X,:])

        # The algorithm pseudocode doesn't say to do this, but I believe
        # Bonchi et al really mean to sample a node that specifically shares a
        # colorx edge with x, otherwise, why both sampling a pivot node in this way?
        neighbsX = findall(x->x == colorx,Anew[X,:])

        # Consider only neighbors of x, and look only at the number of colorx edges they are indicent to
        # @show neighbsX, colorx
        yweights = aweights(ColorDegree[neighbsX,colorx])

        Y = StatsBase.sample(neighbsX,yweights)
        l = A[X,Y]
        @assert(l == colorx)

        # Neighbors of x and y of color l
        Nx = findall(x->x == l,Anew[X,:])
        Ny = findall(x->x == l,Anew[Y,:])
        Nboth = intersect(Nx, Ny)

        # All neighbors of y
        neighbsY = findall(x->x == colorx,Anew[Y,:])

        # Form a preliminary cluster
        newclus = [X;Y]
        cluster[X] = clusnum
        cluster[Y] = clusnum
        node2color[X] = l
        node2color[Y] = l
        assigned[X] = 1
        assigned[Y] = 1
        for t in Nboth
            push!(newclus,t)
            cluster[t] = clusnum
            assigned[t] = 1
            node2color[t] = l
        end

        start = time()
        # Now add new nodes into the cluster
        stilladding = true
        while stilladding


            stilladding = false
            Axn = A[X,newclus]
            Ayn = A[Y,newclus]
            C_toadd = Vector{Int64}()
            neighbsX_notin_C = setdiff(neighbsX,newclus)
            for Z in neighbsX_notin_C

                cs = A[X,Z]    # the color they share

                # From the current cluster, which nodes share a cs edge with X?
                local_friends_of_x = findall(x->x==cs,Axn)
                friends_of_x = newclus[local_friends_of_x]

                # From the current cluster, which nodes share a cs edge with Z?
                local_friends_of_z = findall(x->x==cs,A[Z,newclus])
                friends_of_z = newclus[local_friends_of_z]

                # Check if there is even one node in the overlap
                overlapXZ = length(intersect(friends_of_x,friends_of_z))
                if overlapXZ > 0
                    push!(C_toadd,Z)
                end
            end

            neighbsY_notin_C = setdiff(neighbsY,newclus)
            for Z in neighbsY_notin_C

                cs = A[Y,Z]    # the color they share

                # From the current cluster, which nodes share a cs edge with Y?
                local_friends_of_y = findall(x->x==cs,Ayn)
                friends_of_y = newclus[local_friends_of_y]

                # From the current cluster, which nodes share a cs edge with Z?
                local_friends_of_z = findall(x->x==cs,A[Z,newclus])
                friends_of_z = newclus[local_friends_of_z]

                # Check if there is even one node in the overlap
                overlapYZ = length(intersect(friends_of_y,friends_of_z))
                if overlapYZ > 0
                    push!(C_toadd,Z)
                end
            end

            if length(C_toadd) > 0
                stilladding = true
                for t in C_toadd
                    push!(newclus,t)
                    cluster[t] = clusnum
                    assigned[t] = 1
                    node2color[t] = l
                end
            end

        end
        time2 += time()-start
        clusnum += 1

        # Delete all edges involving nodes that were just clustered
        nonclus = setdiff(1:n,newclus)
        Anew[newclus,:] .= 0
        Anew[:,newclus] .= 0
        dropzeros!(Anew)

        # Find the new edge list (this is very inefficient)
        I,J,EdgeColors = findnz(triu(Anew))
        EdgeList = [I J]
        M = size(EdgeList,1)

    end

    nonlabeled = findall(iszero,assigned)

    for v in nonlabeled
        cluster[v] = clusnum
        clusnum += 1
        # Get all the edges of node v, and assign it a color based on naive-majority
        Vv = nonzeros(A[v,:])
        if length(Vv) > 0
            node2color[v] = StatsBase.mode(Vv)
        else
            node2color[v] = rand(1:K)
        end
    end

    # println("LCB: $time1 \t $time2")

    return round.(Int64,cluster), round.(Int64,node2color)
end

function MajorityVote(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64,k::Int64)
    return NaiveLabel(EdgeList,EdgeColors,n,k)
end

function MajorityVote(EdgeList::Vector{Vector{Int64}},EdgeColors::Array{Int64,1},n::Int64,k::Int64)
    return NaiveLabel(EdgeList,EdgeColors,n,k)
end

# This assumes every edge has only one color
function NaiveLabel(EdgeList::Array{Int64,2},EdgeColors::Array{Int64,1},n::Int64,k::Int64)
    naive_c = zeros(Int64,n)
    A = sparse(EdgeList[:,1],EdgeList[:,2],EdgeColors,n,n)
    @assert(istriu(A))
    A = A+A'
    for i = 1:n
        cneighbs = nonzeros(A[i,:])
        if length(cneighbs) > 0
            color = round(Int64,StatsBase.mode(cneighbs))
            if color > 10
                @show color, cneighbs
            end
            naive_c[i] = color
        else
            naive_c[i] = rand(1:k)
        end
    end

    return naive_c
end

function NaiveLabel(EdgeList::Vector{Vector{Int64}},EdgeColors::Vector{Int64},n::Int64,k::Int64,default::Int64=1)
    ColorDegree = zeros(k,n)
    M = length(EdgeList)
    for t = 1:length(EdgeList)
        edge = EdgeList[t]
        color = EdgeColors[t]
        for node in edge
            ColorDegree[color,node] += 1
        end
    end

    naive_c = zeros(Int64,n)
    for i = 1:n
        colorvec = ColorDegree[:,i]
        if maximum(colorvec) > 0
            naive_c[i] = round(Int64,argmax(colorvec)[1])
        else
            naive_c[i] = rand(1:default)
        end
    end
    return naive_c
end

function get_color_degree(EdgeList::Vector{Vector{Int64}},EdgeColors::Vector{Int64},n::Int64)
    k = round.(Int64,maximum(EdgeColors))
    ColorDegree = zeros(k,n)
    M = length(EdgeColors)
    for t = 1:M
        edge = EdgeList[t]
        color = EdgeColors[t]
        for node in edge
            ColorDegree[color,node] += 1
        end
    end
    return ColorDegree
end


function get_color_degree(EdgeList::Array{Int64,2},EdgeColors::Vector{Int64},n::Int64)
    k = round.(Int64,maximum(EdgeColors))
    ColorDegree = zeros(k,n)
    M = length(EdgeColors)
    for t = 1:M
        n1 = EdgeList[t,1]
        n2 = EdgeList[t,2]
        color = EdgeColors[t]
        ColorDegree[color,n2] += 1
        ColorDegree[color,n1] += 1
    end
    return ColorDegree
end

"""
RemoveColor

Remove all edges of a certain color in an instance of Chromatic Clustering.
This is useful when there is a very common type of edge that too many nodes
participate in, which makes the objective focus too much on satisfying that
objective

k is the color to remove

Note that this also gets rid of many nodes that may not be in any type color but
the removed color. If this is the case, then you need to also call RenameNodes
in order to re-numbe them to be from 1 to number of unique nodes.
"""
function RemoveColor(EdgeList::Vector{Vector{Int64}}, EdgeColor::Vector{Int64}, k::Int64)

    NewList = Vector{Vector{Int64}}()
    NewColor = Vector{Int64}()
    Nodes = Set{Int64}()
    for i = 1:length(EdgeColor)
        edge = EdgeList[i]
        if EdgeColor[i] != k
            push!(NewList,edge)
            push!(NewColor,EdgeColor[i])
            for node in edge
                push!(Nodes,node)
            end
        end
    end
    return NewList, NewColor, sort(unique(Nodes))
end

"""
Rename Nodes

Sometimes hyperedges are made up of a list of node labels, where the nodes labels
aren't 1 through k for some integer k. In this case, we'd like to relabel the nodes.
"""
function RenameNodes(EdgeList::Vector{Vector{Int64}})
    # First, get the list of unique node ids
    nodes = Set{Int64}()
    for i = 1:length(EdgeList)
        for node in EdgeList[i]
            push!(nodes,node)
        end
    end
    new2old = sort(unique(nodes))

    # Get a dictionary that maps from old labels to new labels
    old2new = Dict()
    for t = 1:length(new2old)
        old2new[new2old[t]] =  t
    end

    NewList = Vector{Vector{Int64}}()
    # Go through and rename the nodes in the hyperedges
    for i = 1:length(EdgeList)
        edge = EdgeList[i]
        EdgeVec = Vector{Int64}()
        for node in edge
            Jnew = old2new[node]
            push!(EdgeVec,Jnew)
        end
        push!(NewList,EdgeVec)
    end

    return NewList, old2new, new2old
end


function MaxHyperedgeSize(EdgeList::Vector{Vector{Int64}})
    M = length(EdgeList)
    msize = 0
    for j = 1:M
        mnew = length(EdgeList[j])
        if mnew > msize
            msize = mnew
        end
    end
    return msize
end

# Get the maximum k elements in a vector, and their locations in the vector.
function maxk(a, k)
    b = partialsortperm(a, 1:k, rev=true)
    b = vec(round.(Int64,b))
    return b, a[b]
end

"""
RemoveByTotalDegree

Remove nodes from a labeled hypergraph based on whether they participate in
too many nodes of a different color than the main color they participate in.
"""
function RemoveByTotalDegree(EdgeList::Vector{Vector{Int64}},EdgeColors::Vector{Int64},n::Int64,beta::Union{Float64,Int64})

    M = length(EdgeList)
    # Get the color degree for each node
    Cdeg = get_color_degree(EdgeList,EdgeColors,n)

    # For each node, compute the total degree (all edges of all colors),
    # minus the maximum degree for any one single edge type.
    # This lower bounds the number of mistakes that will be made at this node.
    Tdeg = vec(sum(Cdeg,dims = 1))
    Mdeg = vec(maximum(Cdeg,dims = 1))
    Mistake_Bound = Tdeg - Mdeg

    # Keep everything below the threshold of minimum mistakes.
    # This provides a map from new indices to their old indices.
    new2old = findall(x->x<=beta,Mistake_Bound)

    # Get a dictionary that maps from old labels to new labels
    old2new = Dict()
    for t = 1:length(new2old)
        old2new[new2old[t]] = t
    end

    NewList = Vector{Vector{Int64}}()
    NewColors = Vector{Int64}()
    NewListInds = Vector{Int64}()

    # Go through and rename the nodes in the hyperedges
    for i = 1:length(EdgeList)
        edge = EdgeList[i]              # the previous hyperedge
        EdgeVec = Vector{Int64}()       # initialize a new hyperedge
        color = EdgeColors[i]           # the hyperedge color

        # For each node in this old hyperedge
        for node in edge

            # If it is in the set of new nodes, then add it to the new edge
            Jnew = get(old2new,node,0)
            if Jnew > 0
                Jnew = old2new[node]
                push!(EdgeVec,Jnew)
            end
        end

        # If there are at least two nodes in the hyperedge, then keep it
        # in the new edge list
        if length(EdgeVec) > 1
            push!(NewList,EdgeVec)
            push!(NewColors,color)
            push!(NewListInds,i)
        end
    end

    n_new = length(new2old)

    Mnew = length(NewList)
    println("Previously: $M edges, $n nodes. Now: $Mnew edges, $n_new nodes")
    # Return the new Edgelist, color list, the new n,
    # and the set of old node labels that are in the new set (which will be
    # indexed starting from one in the returned hyperedge list), and a map
    # from old to new.
    return NewList, NewColors, n_new, new2old, old2new, NewListInds

end


"""
CountInclusion

Given a list of nodes and a list of edges, simply count the number of edges
that include all of the nodes in the node list.
"""
function CountInclusion(nodes::Vector{Int64}, EdgeList::Vector{Vector{Int64}})
    count = 0
    for edge in EdgeList
        remain = setdiff(nodes,edge)
        @show length(remain)
        if length(remain) == 0
            count += 1
        end
    end
    return count
end


"""
CountParticipation

Given a list of nodes and a list of edges, simply count the number of edges
that include at least one of the nodes in the node list.
"""
function CountParticipation(nodes::Vector{Int64}, EdgeList::Vector{Vector{Int64}})

    count = 0
    for edge in EdgeList
        remain = intersect(nodes,edge)
        if length(remain) > 0
            count += 1
        end
    end

    return count
end

"""
Calculate the ratio of nodes that contribute towards the objective.
"""
function NodeUse(EdgeList::Union{Array{Int64,2},Vector{Vector{Int64}}},EdgeColors::Array{Int64,1},c::Vector)
    n = length(c); Mistakes = 0
    Used = Set{Int64}()
    for i = 1:size(EdgeList,1)
        if size(EdgeList,2) == 2
            edge = EdgeList[i,:]
        else
            edge = EdgeList[i]
        end
        for v in edge
            if c[v] != EdgeColors[i]
                Mistakes += 1
                break
            end
        end
    end
    return Mistakes
end

"""
Convert a Labeled Hypergraph to a Labeled graph. For each pair of nodes that
share in at least one hyperedge, draw an edge with label corresponding to the
most common edge type they participate in together.
"""
function Hypergraph2Graph(EdgeList::Vector{Vector{Int64}},EdgeLabels::Vector{Int64})
    K = maximum(EdgeLabels)
    ELmap = Dict()

    for t = 1:length(EdgeList)
        edge = EdgeList[t]
        lab = EdgeLabels[t]
        m = length(edge)
        for i = 1:m
            ei = edge[i]
            for j = i+1:m
                ej = edge[j]
                if ei != ej
                    I = min(ei,ej)
                    J = max(ej,ei)
                    if haskey(ELmap,[I,J])
                        ELmap[[I,J]][lab]+=1
                    else
                        ELmap[[I,J]] = zeros(Int64,K)
                        ELmap[[I,J]][lab]+=1
                    end
                end
            end
        end
    end

    ks = collect(keys(ELmap))
    NewList = zeros(length(ks),2)
    NewLabels = zeros(Int64,length(ks))
    for i = 1:length(ks)
        NewList[i,:] = ks[i]
        NewLabels[i] = argmax(ELmap[ks[i]])

    end

    return round.(Int64,NewList), NewLabels
end
