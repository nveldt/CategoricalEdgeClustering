using MAT
using Plots
using Statistics

# Load results
name = "Graph_Output/Synthetic_Large_K_Graph_Results.mat"
plotname = "Plots/Graph_Large_K_Accuracy.pdf"
plotname2 = "Plots/Graph_Large_K_ARI.pdf"

mat = matread(name)

lcb_accs = median(mat["lcballs_accs"],dims=2)
cb_accs = median(mat["cballs_accs"],dims=2)
maj_accs = median(mat["maj_accs"],dims=2)
lp_accs = median(mat["lp_accs"],dims=2)
lcb_aris = median(mat["lcballs_aris"],dims=2)
cb_aris = median(mat["cballs_aris"],dims=2)
maj_aris = median(mat["maj_aris"],dims=2)
lp_aris = median(mat["lp_aris"],dims=2)
Ks = mat["Ks"]

## Accuracy scores plot
l_place = :best
x_label = "K = Number of Clusters"
y_label = "Node Label Accuracy"
s1 = 300
s2 = 250
ms = 5
lw = 2
title = ""
plot(Ks,lp_accs, title = title,
    labels = "LP", grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw, markerstrokewidth = 0,
    color = :blue, markershape = :circle, markersize = ms)

plot!(Ks,maj_accs,
    labels = "MV",linewidth = lw, markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)

plot!(Ks,lcb_accs,
    labels = "LCB",linewidth = lw, markerstrokewidth = 0,
    color = :green, markershape = :circle, markersize = ms)

plot!(Ks,cb_accs,
    labels = "CB",linewidth = lw, markerstrokewidth = 0,
    color = :purple, markershape = :circle, markersize = ms)

savefig(plotname)

## ARI scores plot
x_label = "K = Number of Clusters"
y_label = "ARI"
s1 = 300
s2 = 250
ms = 5
lw = 2
title = ""
plot(Ks,lp_aris, title = title,
    labels = "LP", grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw, markerstrokewidth = 0,
    color = :blue, markershape = :circle, markersize = ms)

plot!(Ks,maj_aris,
    labels = "MV",linewidth = lw, markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)

plot!(Ks,lcb_aris,
    labels = "LCB",linewidth = lw, markerstrokewidth = 0,
    color = :green, markershape = :circle, markersize = ms)

plot!(Ks,cb_aris,
    labels = "CB",linewidth = lw, markerstrokewidth = 0,
    color = :purple, markershape = :circle, markersize = ms)

savefig(plotname2)
