##

mat = matread("Output/Email_Results.mat")
lp_accs = median(mat["lp_accs"],dims = 2)
maj_accs = median(mat["naive_accs"],dims = 2)
cb_accs = median(mat["cballs_accs"],dims = 2)
lcb_accs = median(mat["lcballs_accs"],dims = 2)

## Label accuracy plot
using Plots
x_label = "w = Pr(edge label is wrong)"
y_label1 = "Node Label Accuracy"

l_place = :false
ann = [.5, .5, text("Hi")]
s1 = 300
s2 = 250
ms = 5
lw = 2
alphas = 0.0:.05:0.75
title = ""
plot(alphas,lp_accs, title = title,
    labels = "LP",
    ylim = [0,1],grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label1, legend = l_place,
    linewidth = lw,markerstrokewidth = 0,
    color = :blue, markershape = :circle, markersize = ms)

annotate!(.6, .85, text("LP",font(9,:blue)))

plot!(alphas,maj_accs,
    labels = "MV",linewidth = lw,markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)
annotate!(.7, .62, text("MV",font(9,:red)))

plot!(alphas,lcb_accs,
    labels = "LazyCB",    linewidth = lw,markerstrokewidth = 0,
    color = :green, markershape = :circle, markersize = ms)
annotate!(.18, .67, text("LCB",font(9,:green)))

plot!(alphas,cb_accs,
    labels = "CB",    linewidth = lw,markerstrokewidth = 0,
    color = :purple, markershape = :circle, markersize = ms)
annotate!(.2, .27, text("CB",font(9,:purple)))


savefig("Email_Plot.pdf")
