using JLD
using MAT
include("../src/EdgeCatClusAlgs.jl")
include("../src/lp_isocut.jl")

datasets = ["Brain","MAG-10","Cooking","DAWN", "Walmart-Trips"]

numdata = length(datasets)
dataset_stats = zeros(numdata,4)
lp_stats = zeros(numdata,4)
maj_stats = zeros(numdata,4)
cb_stats = zeros(numdata,4)
lcb_stats = zeros(numdata,4)
iso_stats = zeros(numdata,4)

## Load all data from jld files and run all algorithms
#
# Note: for Walmart and MAG-10, some algorithms take a very long time.
# In practice we ran experiments for several algorithms and datasets separately.
# Output is stored in various .mat files in the "Output" folder
for i = 1:length(datasets)
    dataset = datasets[i]

    data = load("../data/JLD_Files/"*dataset*".jld")
    EdgeColors = data["EdgeColors"]
    EdgeList = data["EdgeList"]
    n = data["n"]
    M = length(EdgeColors)
    msize = MaxHyperedgeSize(EdgeList)
    k = maximum(EdgeColors)
    println("Hypergraph: "*dataset*" has $M edges $n nodes, and $msize maximum order, $k colors")

    dataset_stats[i,:] = [n,M,msize,k]

    ## Solve the Chromatic Clustering Objective
    start = time()
    LPval, X, runtime = EdgeCatClusGeneral(EdgeList,EdgeColors,n,false,1)
    lp_run = round(time()-start,digits=2)

    # Round the clustering
    C = rowmin(X)
    lp_c = C[:,2]
    lp_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,lp_c)
    lp_ratio = lp_mistakes/LPval
    lp_edgesat = 1 - lp_mistakes/M
    lp_stats[i,:] = [lp_mistakes, lp_ratio, lp_edgesat, lp_run]

    println("Running Majority")
    ## Majority vote
    start = time()
    maj_c = MajorityVote(EdgeList,EdgeColors,n,k)
    maj_run = time()-start
    maj_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,maj_c)
    maj_ratio = maj_mistakes/LPval
    maj_edgesat = 1 - maj_mistakes/M
    maj_stats[i,:] = [maj_mistakes, maj_ratio, maj_edgesat, maj_run]

    println("Running IsoCut")
    ## Isolating Cut Heuristic
    start = time()
    NewList, NewColors, NewWeights = IsocutHyper2Graph(EdgeList,EdgeColors)
    iso_c = IsolatingCutLP(NewList,NewColors,NewWeights,n,maj_c,0)
    iso_run = time()-start
    iso_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,iso_c)
    iso_ratio = iso_mistakes/LPval
    iso_edgesat = 1-iso_mistakes/M
    iso_stats[i,:] = [iso_mistakes, iso_ratio, iso_edgesat, iso_run]

    # Convert the hypergraphs to graphs for the Chromatic CC algorithms
    NewList, NewLabels = Hypergraph2Graph(EdgeList,EdgeColors)

    println("Running Chromatic Balls")
    ## Chromatic Balls Algorithm
    start = time()
    cb_cluster, cb_c = ChromaticBalls(NewList,NewLabels,n)
    cb_run = time()-start
    cb_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,cb_c)
    cb_ratio = cb_mistakes/LPval
    cb_edgesat  = 1-cb_mistakes/M
    cb_stats[i,:] = [cb_mistakes, cb_ratio, cb_edgesat, cb_run]

    println("Running Lazy Chromatic Balls")
    ## Lazy Chromatic Balls Algorithm
    start = time()
    lcb_cluster, lcb_c = LazyChromaticBalls(NewList,NewLabels,n)
    lcb_run = time()-start
    lcb_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,lcb_c)
    lcb_ratio = lcb_mistakes/LPval
    lcb_edgesat  = 1-lcb_mistakes/M
    lcb_stats[i,:] = [lcb_mistakes, lcb_ratio, lcb_edgesat, lcb_run]

    matwrite("Output/"*dataset*"_results.mat", Dict("LPval"=>LPval, "X" => X,
     "lp_run" => lp_run, "lp_c" => lp_c, "maj_c" => maj_c, "lp_mistakes" => lp_mistakes,
     "lp_ratio" => lp_ratio, "lp_edgesat" => lp_edgesat, "maj_mistakes" => maj_mistakes,
     "iso_edgesat" => iso_edgesat, "iso_ratio" => iso_ratio, "iso_c" => iso_c, "iso_run"=> iso_run,
     "cb_edgesat" => cb_edgesat, "cb_ratio"=> cb_ratio, "cb_c" => cb_c,
     "lcb_edgesat" => lcb_edgesat, "lcb_ratio"=> lcb_ratio, "lcb_c" => lcb_c,
     "cb_run"=>cb_run, "lcb_run"=>lcb_run,
     "maj_ratio" => maj_ratio, "maj_edgesat" => maj_edgesat,"maj_run"=>maj_run))

end

matwrite("Output/All_datasets_results.mat", Dict("lp_stats" => lp_stats,
"maj_stats" => maj_stats, "datasets" => datasets,
"dataset_stats"=>dataset_stats, "iso_stats" => iso_stats,
"cb_stats"=>cb_stats, "lcb_stats"=> lcb_stats))
