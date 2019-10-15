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
ws = .05:.1:.5 #ws = 0.0:.1:.5
lps = zeros(n,length(ws))
majs = zeros(n,length(ws))
lp_stats = zeros(4,length(ws))
maj_stats = zeros(3,length(ws))
M = length(EdgeList)
for ii= 1:length(ws)
    w = ws[ii]

    # Assign synthetic edge labels
    EdgeColors = zeros(Int64,length(EdgeList))
    K = maximum(NodeLabels)
    correct_edges = 0
    for t = 1:length(EdgeList)
        EdgeColors[t] = rand(1:K)
        edge = EdgeList[t]
        labs = NodeLabels[edge]
        maxl = maximum(labs)
        minl = minimum(labs)
        if maxl == minl && rand(1)[1] < 1-w
            EdgeColors[t] = round(Int64,maxl)
            correct_edges += 1
        end
    end

    ## LP-based algorithms
    start = time()
    LPval, X, runtime = EdgeCatClusGeneral(EdgeList,EdgeColors,n,false,1)
    lp_run = time()-start

    lp_c, RoundScore, lp_approx = SimpleRound(EdgeList,EdgeColors,X,LPval)
    lp_approx = RoundScore/LPval
    lp_c = round.(Int64,lp_c)
    lps[:,ii] = lp_c
    lp_sat = 1-RoundScore/M
    lp_acc = 1-count(!iszero,lp_c-NodeLabels)/n
    lp_stats[:,ii] = [lp_run;lp_sat;lp_approx;lp_acc]
    ##
    maj_c = MajorityVote(EdgeList,EdgeColors,n,K)
    majs[:,ii] = maj_c
    maj_obj = EdgeCatClusObj(EdgeList,EdgeColors,maj_c)
    maj_sat = 1-maj_obj/M
    maj_approx = maj_obj/LPval
    maj_acc = 1-count(!iszero,maj_c-NodeLabels)/n
    maj_stats[:,ii] = [maj_sat;maj_approx;maj_acc]

    matwrite("Output/Departments_w_$w"*".mat", Dict("maj_stats"=>maj_stats,
    "lp_stats" => lp_stats, "EdgeColors" => EdgeColors, "LPval" => LPval))

end
