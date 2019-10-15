using MAT

data = matread("Output/lp_clusterings_2to40.mat")
lps = data["lps"]

include("Temporal_Cluster_Subroutines.jl")

mat = matread("Output/Graclus_clusterings_2to40.mat")
c_grac = round.(Int64,mat["Cmets"])

using MAT
mat = matread("../data/CollegeMsg/college-msg.mat")
EdgeList = round.(Int64,mat["EdgeList"])
EdgeList_Array = EdgeList
EdgeTime = mat["EdgeTime"]
EdgeTime .-= minimum(EdgeTime)
EdgeTime = EdgeTime./3600

## Extract clusterings and measure normalized cut and avg time diff
ks = 2:40

ncuts = zeros(length(ks),2)
Diffs = zeros(length(ks),2)

#
for i = 1:length(ks)
    k = ks[i]
    ec = matread("EdgeColors/EdgeColors_k_$k.mat")
    EdgeColors = ec["EdgeColors"]

    ## LP-round stats for each cluster
    # (average time stamp, standard deviation of inner edge timestamps, cluster sizes, and time difference)
    lp_c = lps[:,i]
    sdevs_lp, Avgs_lp, Sz_lp, lp_cut, lp_diff = ClusterCheck(EdgeList_Array, EdgeTime,EdgeColors, lp_c)

    # Graclus stats for each cluster
    # (average time stamp, standard deviation of inner edge timestamps, cluster sizes, and time difference)
    cg = c_grac[:,i]
    sdevs_grac, Avgs_grac, Sz_grac, grac_cut, grac_diff = ClusterCheck(EdgeList_Array, EdgeTime,EdgeColors, cg)

    ncuts[i,:] = [Ncut(A,lp_c),Ncut(A,cg)]
    Diffs[i,:] = [lp_diff, grac_diff]
end

##

xs = 2:40
x_label = "k = number of clusters"
y_label = "1/k *(Norm. Cut) "
l_place = :bottomright
s1 = 300
s2 = 250
ms = 5
lw = 2
title = ""
plot(xs,ncuts[:,1]./ks, title = title,
    labels = "LP-round",
    grid = false, size = (s1,s2),
    xlabel = x_label, xlim = [0,40], ylabel = y_label, legend = l_place,
    linewidth = lw,markerstrokewidth = 0,markershape = :circle,
    color = :blue, markersize = ms)

plot!(xs,ncuts[:,2]./ks,linewidth = lw, labels = "Graclus",markershape = :circle,
markerstrokewidth = 0,color = :black, markersize = ms)


savefig("Plots/Normcuts.pdf")

## Plots the Average Time Difference

x_label = "k = number of clusters"
y_label = "Avg Time Difference (hrs)"
l_place = :topright
s1 = 300
s2 = 250
ms = 5
lw = 2
title = ""
plot(xs,Dists[:,1], title = title,
    labels = "LP-round",
    grid = false, size = (s1,s2),markershape = :circle,
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw,markerstrokewidth = 0,
    color = :blue, markersize = ms)

plot!(xs,Dists[:,2],linewidth = lw, linestyle = :solid, labels = "Graclus",
markerstrokewidth = 0,color = :black, markersize = ms,markershape = :circle)

savefig("Plots/Timedev.pdf")
