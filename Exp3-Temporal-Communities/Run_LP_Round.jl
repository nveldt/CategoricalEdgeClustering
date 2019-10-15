include("../src/EdgeCatClusAlgs.jl")
include("Temporal_Cluster_Subroutines.jl")
using MAT
mat = matread("../data/CollegeMsg/college-msg.mat")
EdgeList = round.(Int64,mat["EdgeList"])
EdgeList_Array = EdgeList
EdgeTime = mat["EdgeTime"]

# Normalize time stamps so they start at zero and are measured in hours
EdgeTime .-= minimum(EdgeTime)
EdgeTime = EdgeTime./3600

NewList = Vector{Vector{Int64}}()
for t = 1:length(EdgeTime)
    push!(NewList,EdgeList[t,:])
end
EdgeList = NewList

## Get an undirected, not time stamped version of the dataset.
A = CountEdgeWeight(EdgeList_Array,n)
A = A+A'
for i= 1:n
    A[i,i] = 0
end
A = dropzeros(A)
matwrite("A_sym.mat",Dict("A"=>A))

##
n = 1899
mat = matread("Graclus_clusterings_2to40.mat")
c_grac = round.(Int64,mat["Cmets"])
M = length(EdgeList)

## Change the EdgeTime into EdgeColors by choosing number of pieces to sort it in
ks = 2:40
lps = zeros(n,length(ks))

##
for i = 1:length(ks)
    k = ks[i]
    mintime = minimum(EdgeTime)
    maxtime = maximum(EdgeTime)
    total = length(EdgeTime)

    edgespercolor = round(Int64,floor(total/k))

    ## Define discrete labels by binning time stamps.
    # Form equally sized bins.
    EdgeColors = zeros(length(EdgeTime))
    for i = 1:k
        EdgeColors[((i-1)*edgespercolor+1):i*edgespercolor] .= i
    end

    extra = mod(total,k)
    EdgeColors[total-extra+1:total] .= k
    EdgeColors = round.(Int64,EdgeColors)

    # save the edge labels if desired
    matwrite("EdgeColors/EdgeColors_k_$k.mat",Dict("EdgeColors"=>EdgeColors))

    ## LP-round
    start = time()
    LPval, X, runtime = EdgeCatClusGeneral(EdgeList,EdgeColors,n,false,1)
    LPtime = round(time()-start,digits=2)
    lp_c, lp_mistakes, lp_ratio = SimpleRound(EdgeList,EdgeColors,X,LPval)
    lps[:,i] = lp_c

end

matwrite("Output/lp_clusterings_2to40.mat", Dict("ks" => collect(ks),"lps" => lps))
