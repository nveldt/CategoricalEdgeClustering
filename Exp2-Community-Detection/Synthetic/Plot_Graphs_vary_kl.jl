using Plots
using MAT

# First Experiment
name = "Graph_Output/Synthetic_vary_kl_Graph_Results.mat"
plotname = "Plots/Graph_vary_kl.pdf"

mat = matread(name)

lcballs_accs = mat["lcballs_accs"]
cballs_accs = mat["cballs_accs"]
maj_accs = mat["maj_accs"]
lp_accs = mat["lp_accs"]
Ls = mat["Ls"]

lp_median = median(lp_accs,dims=2)
lc_median = median(cballs_accs,dims=2)
lcb_median = median(lcballs_accs,dims=2)
maj_median = median(maj_accs,dims=2)

xs = Ls

## Label accuracy plot
x_label = "K = Number of Clusters"
y_label = "Node Label Accuracy"

l_place = :best
s1 = 300
s2 = 250
ms = 5
lw = 2
title = ""
plot(xs,lp_median, title = title,
    labels = "LP",
    ylim = [0,1],grid = false, size = (s1,s2), markerstrokewidth = 0,
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw,
    color = :blue, markershape = :circle, markersize = ms)

plot!(xs,maj_median,
    labels = "MV",linewidth = lw, markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)

plot!(xs,lcb_median,
labels = "LCB",    linewidth = lw, markerstrokewidth = 0,
color = :green, markershape = :circle, markersize = ms)

plot!(xs,lc_median,
labels = "CB",    linewidth = lw,markerstrokewidth = 0,
color = :purple, markershape = :circle, markersize = ms)

savefig(plotname)
