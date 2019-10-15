using Statistics

"""
Given a clustering of a temporal graph, check the timestamps for the cluster.
"""
function ClusterCheck(EdgeList, EdgeTimes, EdgeLabels, c)

    # We need to first go through the EdgeList and find
    # which are in the same cluster. If they are cut, then
    # we keep track of the number that are cut
    cut = 0
    M = length(EdgeTimes)
    InnerEdges = zeros(M)
    for t = 1:M
        i = min(EdgeList[t,1],EdgeList[t,2])
        j = max(EdgeList[t,1],EdgeList[t,2])

        # If this edge is inside a cluster, label it with that cluster
        if c[i] == c[j]
            InnerEdges[t] = c[i]
        else
            # Otherwise, it is a cut edge
            cut += 1
        end
    end

    k = round(Int64,maximum(c))
    sdevs = zeros(k)
    Avgs = zeros(k)
    Sz = zeros(k)
    for i = 1:k

        # Get the set of inner edges
        inds = findall(x->x==i,InnerEdges)

        # Check their time stamps
        times = EdgeTimes[inds]

        # Find average and variance
        Avgs[i] = mean(times)
        sdevs[i] = std(times)
        Sz[i] = length(inds)
    end

    # Now calculate the average distance from the mean
    avgtimediff = 0
    for i = 1:k

        # Get the set of inner edges and their time stamps
        inds = findall(x->x==i,InnerEdges)
        times = EdgeTimes[inds]
        for j = 1:length(times)
            avgtimediff += abs(times[j]-Avgs[i])   # time difference between the interior edge, and it's average point
        end

    end
    avgtimediff = avgtimediff/length(InnerEdges)

    return sdevs, Avgs, Sz, cut, avgtimediff
end

function Ncut(A,c)
    d = sum(A,dims = 2)
    ncut = 0
    for i = 1:maximum(c)
        inds = findall(x->x==i,c)
        vol = sum(d[inds])
        Asub = A[inds,inds]
        edges = sum(nonzeros(Asub))/2
        cut = vol-edges
        if vol > 0
            ncut+=cut/vol
        end
    end
    return ncut
end
# All stats please

function AllStats(EdgeList, EdgeColors, EdgeTimes, LPval, n, c)
    M = length(EdgeColors)
    obj = ChromeClusObj(EdgeList,EdgeColors,c)
    ratio = obj/LPval
    sat = 1-obj/M

    sdevs, Avgs, Sz, cut, avgtimediff = ClusterCheck(EdgeList, EdgeTimes, EdgeColors, c)

    return ratio, sat, cut, avgtimediff
end

# Count the number of interactions two people have
function CountEdgeWeight(EdgeList,n)

    A = spzeros(n,n)
    for t = 1:size(EdgeList,1)
        i = min(EdgeList[t,1],EdgeList[t,2])
        j = max(EdgeList[t,1],EdgeList[t,2])
        A[i,j] += 1
    end

    return A

end
