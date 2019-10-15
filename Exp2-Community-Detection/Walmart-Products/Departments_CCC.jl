include("../../src/EdgeCatClusAlgs.jl")
using MAT
using JLD
data = load("../../data/JLD_Files/Walmart_NodeLabels.jld")
EdgeList = data["EdgeList"]
Categories = data["Categories"]
Departments = data["Departments"]
NodeDepartment = data["NodeDepartment"]
NodeCategories = data["NodeLabels"]
n = data["n"]

##
NodeLabels = NodeDepartment
ws = 0.0:0.05:0.5
lps = zeros(n,length(ws))
majs = zeros(n,length(ws))
cbs = zeros(n,length(ws))
lcbs = zeros(n,length(ws))

lp_stats = zeros(4,length(ws))
cb_stats = zeros(4,length(ws))
lcb_stats = zeros(4,length(ws))
maj_stats = zeros(3,length(ws))
K = maximum(NodeLabels)
M = length(EdgeList)
for ii= 1:length(ws)
    w = ws[ii]

    # Having already run experiments with the LP, just extract the LPval
    mat = matread("Output/Departments_w_$w.mat")
    LPval = mat["LPval"]
    EdgeColors = mat["EdgeColors"]

    M = length(EdgeColors)
    NewList, NewLabels = Hypergraph2Graph(EdgeList,EdgeColors)

    maj_c = MajorityVote(EdgeList,EdgeColors,n,K)
    majs[:,ii] = maj_c
    maj_obj = EdgeCatClusObj(EdgeList,EdgeColors,maj_c)
    maj_sat = 1-maj_obj/M
    maj_approx = maj_obj/LPval
    maj_acc = 1-count(!iszero,maj_c-NodeLabels)/n
    maj_stats[:,ii] = [maj_sat;maj_approx;maj_acc]

    println("Running Chromatic Balls")
    start = time()
    cb_cluster, cb_c = ChromaticBalls(NewList,NewLabels,n)
    cb_run = time()-start
    cbs[:,ii] = cb_c
    cb_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,cb_c)
    cb_approx = cb_mistakes/LPval
    cb_sat = 1-cb_mistakes/M
    cb_acc = 1-count(!iszero,cb_c-NodeLabels)/n
    cb_stats[:,ii] = [cb_run;cb_sat;cb_approx;cb_acc]

    println("Running Lazy Chromatic Balls")
    start = time()
    lcb_cluster, lcb_c = LazyChromaticBalls(NewList,NewLabels,n)
    lcb_run = time()-start
    lcbs[:,ii] = lcb_c
    lcb_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,lcb_c)
    lcb_approx = lcb_mistakes/LPval
    lcb_sat = 1-lcb_mistakes/M
    lcb_acc = 1-count(!iszero,lcb_c-NodeLabels)/n
    lcb_stats[:,ii] = [lcb_run;lcb_sat;lcb_approx;lcb_acc]

    matwrite("Output/Departments_w_$w"*"_CCC.mat", Dict("maj_stats"=>maj_stats,
    "EdgeColors" => EdgeColors,"lcb_c"=>lcb_c, "cb_c"=>cb_c))

end
