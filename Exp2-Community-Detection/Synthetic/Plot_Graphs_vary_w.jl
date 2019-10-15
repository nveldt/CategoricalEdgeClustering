using Plots
using MAT

# First Experiment
name = "Graph_Output/Synthetic_vary_w_Graph_Results.mat"
plotname = "Plots/Graph_vary_w.pdf"

mat = matread(name)

lcballs_accs = mat["lcballs_accs"]
cballs_accs = mat["cballs_accs"]
maj_accs = mat["maj_accs"]
lp_accs = mat["lp_accs"]
ws = mat["ws"]

lp_median = StatsBase.median(lp_accs,dims=2)
lc_median = StatsBase.median(cballs_accs,dims=2)
lcb_median = StatsBase.median(lcballs_accs,dims=2)
maj_median = StatsBase.median(maj_accs,dims=2)

## Label accuracy plot
x_label = "w = Pr(in-cluster edge is wrong)"
y_label = "Node Label Accuracy"

l_place = :best
s1 = 300
s2 = 250
ms = 5
lw = 2
n = 1000
title = ""
plot(alphas,lp_median, title = title,
    labels = "LP", markerstrokewidth = 0,
    ylim = [0,1],grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw,
    color = :blue, markershape = :circle, markersize = ms)

plot!(alphas,n_median, markerstrokewidth = 0,
    labels = "MV",linewidth = lw,
    color = :red, markershape = :circle, markersize = ms)

plot!(alphas,lcb_median, markerstrokewidth = 0,
labels = "LCB",    linewidth = lw,
color = :green, markershape = :circle, markersize = ms)

plot!(alphas,lc_median,markerstrokewidth = 0,
labels = "CB",    linewidth = lw,
color = :purple, markershape = :circle, markersize = ms)

savefig(plotname)
