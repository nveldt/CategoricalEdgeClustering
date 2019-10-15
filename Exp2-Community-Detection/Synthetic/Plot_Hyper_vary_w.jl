using Plots
using MAT

name = "Hypergraph_Output/Hypergraph_all_vary_w_results.mat"
plotname = "Plots/Hypergraph_vary_w.pdf"
mat = matread(name)
ws = mat["ws"]
maj_accs = mat["maj_accs"]
lp_accs = mat["lp_accs"]
lcb_accs = mat["lcb_accs"]
cb_accs = mat["cb_accs"]

lps = median(lp_accs,dims = 2)
majs = median(maj_accs, dims = 2)
lcbs = median(lcb_accs, dims = 2)
cbs = median(cb_accs,dims = 2)


## Label accuracy plot
x_label = "w = Pr(in-cluster edge is wrong)"
y_label = "Node Label Accuracy"
l_place = :left
s1 = 300
s2 = 250
ms = 5
lw = 2
title = ""
plot(ws,lps, title = title,
    labels = "LP",
    ylim = [0,1],grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw, markerstrokewidth = 0,
    color = :blue, markershape = :circle, markersize = ms)

plot!(ws,majs,
    labels = "MV",linewidth = lw, markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)

plot!(ws,lcbs,
    labels = "LCB",linewidth = lw, markerstrokewidth = 0,
    color = :green, markershape = :circle, markersize = ms)

plot!(ws,cbs,
    labels = "CB",linewidth = lw, markerstrokewidth = 0,
    color = :purple, markershape = :circle, markersize = ms)



savefig(plotname)
