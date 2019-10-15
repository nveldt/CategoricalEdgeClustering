## This can be either the optimal solution or not, and it takes hyperedges
# and weights if desired
function EdgeCatClusWeighted(EdgeList::Vector{Vector{Int64}},EdgeColors::Array{Int64,1},EdgeWeights::Array{Float64,1},n::Int64,optimalflag::Bool= false,outputflag::Int64=0)

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
IsolatingCutLP

An LP-based implementation of the isolating cut heuristic.
"""
function IsolatingCutLP(EdgeList::Vector{Vector{Int64}},EdgeColors::Vector{Int64},EdgeWeights::Vector{Float64},n::Int64,default,outputflag::Int64=0)
    k = maximum(EdgeColors)
    TwoWayCuts = zeros(n,k)
    M = length(EdgeColors)
    scores = zeros(k)

    # Solve the problem for index i
    for i = 1:k
        println("at $i of $k two-way problems for isocut")
        # The same edge may have multiple colors, so we need to be careful.
        NewColors = Vector{Int64}()
        NewList = Vector{Vector{Int64}}()
        NewWeights = Vector{Float64}()
        A_not_i = spzeros(n,n)            # keep track of weight of non-i edges
        for t = 1:length(EdgeList)
            edge = EdgeList[t]
            @assert(length(edge) == 2)
            lab = EdgeColors[t]
            weight = EdgeWeights[t]
            if lab == i
                push!(NewList,edge)
                push!(NewColors,1)
                push!(NewWeights,weight)
            else
                ei = min(edge[1],edge[2])
                ej = max(edge[1],edge[2])
                A_not_i[ei,ej] += weight/2
            end
        end

        # Now sum together all the non-i-labeled edges
        Is, Js, Vs = findnz(A_not_i)
        for tt = 1:length(Vs)
            push!(NewList,[Is[tt]; Js[tt]])
            push!(NewColors,2)
            push!(NewWeights,Vs[tt]/2)
            @assert(length(NewColors) == length(NewWeights))
        end

        LPval, X, runtime = EdgeCatClusWeighted(NewList,NewColors,NewWeights,n, true, outputflag)
        @assert(length(NewColors) == length(NewWeights))

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
    cdeg = get_color_degree(EdgeList,EdgeColors,n)
    colors = zeros(n)
    defaultcount = 0
    for i = 1:n
        # options = findall(!iszero,twc[i,:])
        # @show options, keepinds[options]
        options = keepinds[findall(!iszero,twc[i,:])]
        # if the node is in any of the returned cut sets, choose one at random
        # (in theory, if it is in more than one, then you can put it anywhere
        # and the theoretical approximation guarantee holds)

        if length(options) > 0
            colors[i] = options[rand(1:length(options))]
        elseif cdeg[a,i] > 0
            # otherwise, if the node is attached to the one other terminal,
            # cluster it there
            # colors[i] = a
            colors[i] = default[i]  # just kidding, it works better to just go with the default
                                    # uncomment the above and change if desired
            defaultcount += 1
        else
            # in the last case, go with the majority vote default
            colors[i] = default[i]
            defaultcount += 1
        end
    end
    defaultratio = defaultcount/n
    println("defaultratio= $defaultratio")
    return round.(Int64,colors)
end


"""
Convert a Labeled Hypergraph to a Labeled graph in a way that allows us to use
the isolating cut heuristic to get an approximation.

To do this, replace an r-node hyperedge with label l to an r-node clique where
every edge has label l and weight 2/r. Then, when this is labeled graph is
reduced to a multiway cut problem, we divide the weight by two, so that this
edge contribution ends up being 1/2*(2/r) = 1/r. This is the same as reducing
directly from a hypergraph to a multiway cut problem. I.e., we are doing:

Hypergraph Categorical Clustering --> Weight Graph Categorical clustering (edges adjusted by factor 2/r)
                                  --> Multiway cut problem (edges adjusted by 1/2 times previous 2/r factor)

Same as:
Hypergraph Categorical Clustering --> Multiway cut problem (edges adjusted by factor 1/r)

The reason we do it the first way is because it is easier to then just apply
previous algorithms that are already implemented.

"""
function IsocutHyper2Graph(EdgeList::Vector{Vector{Int64}},EdgeLabels::Vector{Int64})
    K = maximum(EdgeLabels)
    ELmap = Dict()

    for t = 1:length(EdgeList)
        edge = EdgeList[t]
        lab = EdgeLabels[t]
        r = length(edge)
        for i = 1:r
            ei = edge[i]
            for j = i+1:r
                ej = edge[j]
                I = min(ei,ej)
                J = max(ej,ei)
                if haskey(ELmap,[I;J])
                    ELmap[[I;J]][lab] = 2/r
                else
                    ELmap[[I;J]] = zeros(Float64,K)
                    ELmap[[I;J]][lab] += 2/r
                end
            end
        end
    end

    ks = collect(keys(ELmap))

    ## Now for each pair of nodes that is together in some hyperedge,
    # we write down weighted and labeled edges
    NewList = Vector{Vector{Int64}}()
    NewLabels = Vector{Int64}()
    NewWeights = Vector{Float64}()
    for i = 1:length(ks)            # For each pair sharing some edge

        pair = ks[i]                # get the pair of nodes
        @assert(length(pair) == 2)
        weightvec = ELmap[pair]     # get the weight of the edge associated with each color
        @assert(length(weightvec) == K)
        for k = 1:K
            wk = weightvec[k]       # what is the label-k weight?
            if wk > 0               # if > 0,
                push!(NewList,pair) # then we have a weighted edge with label k
                push!(NewLabels,k)
                push!(NewWeights,wk)
            end
        end
    end

    return NewList, NewLabels, NewWeights
end
