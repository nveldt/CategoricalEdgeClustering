# Here we include generators for a number of different types of synthetic
# labeled graph models

using StatsBase
using Random

"""
GenerateCCC_model

This is Algorithm 5 from Bonchi et al. 2015, for generating a synthetic
chromatic graph with ground truth clusters.

* Assign n nodes to K clusters, uniformly at random.
* Assign to each cluster a color from among L colors, chosen uniformly at random
(Importantly, there could be two distinct clusters of the same color, and the number
of colors does not equal the number of clusters).
* For a pair of nodes in the same cluster, assign an edge with probability pin.
    * With probability w, give that edge a random label
    * With probability 1-w, give the edge the color of the cluster
* For a pair of nodes in different clusters, assign an eddge with probability pout.
    * Assign it a random color.

Input:
    n = number of nodes
    L = number of labels
    K = number of clusters
    pin = probability of an edge between nodes in the same cluster
    pout = probability of an edge between nodes in different clusters
    w = probability that an in-cluster edge has the wrong color
"""
function GenerateCCC_model(n::Int64,L::Int64,K::Int64, pin::Float64, pout::Float64, w::Float64)

    Node1 = Vector{Int64}()
    Node2 = Vector{Int64}()
    EdgeColors = Vector{Int64}()

    node2cluster = rand(1:K,n)    # assign a node to a random cluster
    node2cluster, d, b = relabel(node2cluster)    # In case a cluster is empty, relabel
    K = round(Int64,maximum(node2cluster))        # Update the number of clusters there actually are

    cluster2color = rand(1:L,K)   # assign a cluster to a random color
    cluster2color, d, b = relabel(cluster2color)  # In case all colors not used, relabel
    L = round(Int64,maximum(cluster2color))       # Update the number of colors there actually are

    # Record the color of each node
    node2color = zeros(n)
    for k = 1:K
        inds = findall(x->x == k,node2cluster)  # get nodes in cluster k
        color = cluster2color[k]                # get the color of that cluster
        node2color[inds] .= color               # assign the color of each node
    end

    for i = 1:n
        color_i = node2color[i]
        for j = i+1:n

            # If nodes are in the same cluster
            if node2cluster[i] == node2cluster[j]
                clustercolor = node2color[i]
                @assert(node2color[i] == node2color[j])

                # Draw an edge at random
                if rand(1)[1] < pin
                    # If they have an edge, now choose a color
                    if rand(1)[1] < w

                        # With probability w, the color is random
                        color = rand(1:L)
                    else

                        # Otherwise, the color matches the cluster color
                        color = color_i
                    end

                    # Now add that edge
                    push!(Node1,i)
                    push!(Node2,j)
                    push!(EdgeColors,color)
                end
            else

                # If the node are in different clusters, we add a
                # randomly-colored edge with probability pout
                if rand(1)[1] < pout
                    color = rand(1:L)
                    push!(Node1,i)
                    push!(Node2,j)
                    push!(EdgeColors,color)
                end
            end
        end
    end

    EdgeList = [Node1 Node2]
    return EdgeList, EdgeColors, round.(Int64,node2cluster), round.(Int64,node2color)

end


"""
GenerateCCC_model_3

A generalization of the Bonchi et al model to 3-uniform hypergraphs.

* Assign n nodes to K clusters, uniformly at random.
* Assign to each cluster a color from among L colors, chosen uniformly at random
(Importantly, there could be two distinct clusters of the same color, and the number
of colors does not equal the number of clusters).
* For a TRIPLET of nodes in the same cluster, assign an edge with probability pin.
    * With probability w, give that edge a random label
    * With probability 1-w, give the edge the color of the cluster
* For a pair of nodes in different clusters, assign an eddge with probability pout.
    * Assign it a random color.

Input:
    n = number of nodes
    L = number of labels
    K = number of clusters
    pin = probability of an edge between nodes in the same cluster
    pout = probability of an edge between nodes in different clusters
    w = probability that an in-cluster edge has the wrong color
    stable = if true, then make the number of labels equal the number of clusters
"""
function GenerateCCC_model3(n::Int64,L::Int64,K::Int64, pin::Float64, pout::Float64, w::Float64, stable::Bool=false)

    EdgeList = Vector{Vector{Int64}}()
    EdgeColors = Vector{Int64}()

    node2cluster = rand(1:K,n)    # assign a node to a random cluster
    node2cluster, d, b = relabel(node2cluster)    # In case a cluster is empty, relabel
    K = round(Int64,maximum(node2cluster))        # Update the number of clusters there actually are

    if !stable
        cluster2color = rand(1:L,K)   # assign a cluster to a random color
        cluster2color, d, b = relabel(cluster2color)  # In case all colors not used, relabel
        L = round(Int64,maximum(cluster2color))       # Update the number of colors there actually are

        # Record the color of each node
        node2color = zeros(n)
        for k = 1:K
            inds = findall(x->x == k,node2cluster)  # get nodes in cluster k
            color = cluster2color[k]                # get the color of that cluster
            node2color[inds] .= color               # assign the color of each node
        end

    else
        node2color = node2cluster
    end

    for i = 1:n
        color_i = node2color[i]
        for j = i+1:n
            color_j = node2color[j]
            for k = j+1:n

                # If nodes are in the same cluster
                if node2cluster[i] == node2cluster[j] && node2cluster[k] == node2cluster[j]
                    clustercolor = color_i

                    # Draw an edge at random
                    if rand(1)[1] < pin
                        # If they have an edge, now choose a color

                        if rand(1)[1] < w
                            # With probability w, the color is random
                            color = rand(1:L)
                        else
                            # Otherwise, the color matches the cluster color
                            color = color_i
                        end

                        # Now add that edge
                        edge = [i;j;k]
                        push!(EdgeList,edge)
                        push!(EdgeColors,color)
                    end
                else

                    # If the node are in different clusters, we add a
                    # randomly-colored edge with probability pout
                    if rand(1)[1] < pout
                        color = rand(1:L)
                        edge = [i;j;k]
                        push!(EdgeList,edge)
                        push!(EdgeColors,color)
                    end
                end
            end
        end
    end

    return EdgeList, EdgeColors, round.(Int64,node2cluster), round.(Int64,node2color)
end


"""
In some cases algorithms return a bunch of labels for a set of objects, but
    labels do not go from 1 to K = (number of unique labels). This function
    relabels a vector of labels to go from 1 to K, and provides the map back.
"""
function relabel(L_old::Vector)

    current = 1
    n = length(L_old)

    # Get the set of unique old labels.
    # OldLabels[i] returns the old label associated with new label i
    OldLabels = sort(unique(L_old))

    K = length(OldLabels)

    L_new = zeros(n)
    for i = 1:K
        inds = findall(x->x==OldLabels[i],L_old)  # get the index set of OldLabel[i]
        L_new[inds] = i*ones(length(inds))
    end

    # Get a way to quickly map old labels to new labels
    Old2New = Dict()
    for i = 1:K
        Old2New[OldLabels[i]] = i
    end

    return L_new, Old2New, OldLabels
end
