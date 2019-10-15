using MAT
using JLD
data = load("../../data/JLD_Files/Walmart_NodeLabels.jld")
EdgeList = data["EdgeList"]
Categories = data["Categories"]
Departments = data["Departments"]
NodeLabels = data["NodeDepartment"]
n = data["n"]

## Load data for maj and LP. We ran these in two separate experiments,
# so we extract results separately
ws = .05:.1:.45
length(ws)
mat = matread("Output/Departments_w_0.45.mat")
lp_stats = mat["lp_stats"]
maj_stats = mat["maj_stats"]

ws = 0.0:.1:.4
mat2 = matread("Output/Departments_w_0.5.mat")
lp_stats2 = mat2["lp_stats"]
maj_stats2 = mat2["maj_stats"]

## Need to inter-weave solutions
majs = zeros(11)
lps = zeros(11)
for i = 0:4
    majs[2*i+1] = maj_stats2[3,i+1]
    lps[2*i+1] = lp_stats2[4,i+1]
    majs[2*i+2] = maj_stats[3,i+1]
    lps[2*i+2] = lp_stats[4,i+1]
end
majs[11] = maj_stats[3,end]
lps[11] = lp_stats[4,end]

## Load data for LCB and CB
ws = .00:0.05:.5

cbs = zeros(11)
lcbs = zeros(11)
for i = 1:length(ws)
    w = ws[i]
    data = matread("Output/Departments_w_$w"*"_CCC.mat")
    lcb_c = data["lcb_c"]
    cb_c = data["cb_c"]
    maj_stats = data["maj_stats"]
    lcbs[i] = 1-count(!iszero,lcb_c-NodeLabels)/n
    cbs[i] = 1-count(!iszero,cb_c-NodeLabels)/n
end

using Plots
x_label = "w = Pr(edge label is wrong)"
y_label1 = "Node Label Accuracy"
l_place = :false
s1 = 300
s2 = 250
ms = 5
lw = 2
ws = 0.0:.05:0.5
title = ""
plot(ws,lps, title = title, titlefontsize=10,
    labels = "LP", grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label1, legend = l_place,
    linewidth = lw,markerstrokewidth = 0,
    color = :blue, markershape = :circle, markersize = ms)
annotate!(.4, .37, text("LP",font(9,:blue)))

plot!(ws,majs,
    labels = "MV",linewidth = lw,markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)
annotate!(.35, .29, text("MV",font(9,:red)))
plot!(ws,lcbs,
    labels = "LCB",linewidth = lw,markerstrokewidth = 0,
    color = :green, markershape = :circle, markersize = ms)
annotate!(.1, .27, text("LCB",font(9,:green)))
plot!(ws,cbs,
    labels = "CB",linewidth = lw,markerstrokewidth = 0,
    color = :purple, markershape = :circle, markersize = ms)
annotate!(.1, .19, text("CB",font(9,:purple)))

savefig("Walmart-Products-Plot.pdf")
