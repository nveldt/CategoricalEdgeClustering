include("../../src/EdgeCatClusAlgs.jl")
using Clustering
using MAT

println("\n")

n = 1000
Lbound = 15
Kbound = Lbound
ws = 0.0:.05:.75
numws = length(ws)
numtimes = 5
p =  0.005
q = 0.0001

## Store the ARI scores and the color identification accuracy for each method
maj_accs = zeros(numws,numtimes)
lp_accs = zeros(numws,numtimes)
cb_accs = zeros(numws,numtimes)
lcb_accs = zeros(numws,numtimes)

## Store runtimes, just in case we end up caring
lp_run = zeros(numws,numtimes)
cb_run = zeros(numws,numtimes)
lcb_run = zeros(numws,numtimes)

for ii = 1:numws
    w = ws[ii]

    # Save the output of each method for this value of w
    lps = zeros(n,numtimes)
    majs = zeros(n,numtimes)
    cbs = zeros(n,numtimes)
    lcbs = zeros(n,numtimes)
    allaccs = zeros(4,numtimes)

    for times = 1:numtimes

        EdgeList, EdgeColors, node2cluster, node2color = GenerateCCC_model3(n,Lbound,Kbound,p,q,w,true)

        # Reduce to a graph for the chromatic correlation clustering algorithms
        k = maximum(EdgeColors)
        M = length(EdgeColors)
        CCCList, CCCLabels = Hypergraph2Graph(EdgeList,EdgeColors)

        println("$k $M")

        ## LP-based algorithms
        start = time()
        LPval, X, runtime = EdgeCatClusGeneral(EdgeList,EdgeColors,n,false,1)
        lp_run[ii,times] = time()-start

        matwrite("Hypergraphs/Hypergraph_vary_w_$w"*"_$times.mat",
        Dict("EdgeList" => EdgeList, "EdgeColors" => EdgeColors,
        "node2cluster" => node2cluster, "node2color" => node2color,
        "LPval" => LPval,"CCCList"=>CCCList, "CCCLabels"=> CCCLabels))

        lp_c, RoundScore, RoundRatio = SimpleRound(EdgeList,EdgeColors,X,LPval)
        lp_c = round.(Int64,lp_c)
        lps[:,times] = lp_c
        lpacc = 1-count(!iszero,lp_c-node2color)/n
        lp_accs[ii,times] = lpacc

        ## Majority Vote Algorithm
        maj_c = MajorityVote(EdgeList,EdgeColors,n,k)
        majs[:,times] = maj_c
        majacc = 1-count(!iszero,maj_c-node2color)/n
        maj_accs[ii,times] = majacc

        println("Running Chromatic Balls")
        start = time()
        cb_cluster, cb_c = ChromaticBalls(CCCList,CCCLabels,n)
        cbs[:,times] = cb_c
        cbacc = 1-count(!iszero,cb_c-node2color)/n
        cb_accs[ii,times] = cbacc

        println("Running Lazy Chromatic Balls")
        start = time()
        lcb_cluster, lcb_c = LazyChromaticBalls(CCCList,CCCLabels,n)
        lcb_run[ii,times] = time()-start
        lcbs[:,times] = lcb_c
        lcbacc =  1-count(!iszero,lcb_c-node2color)/n
        lcb_accs[ii,times] = lcbacc
        allaccs[:,times] = [lpacc;majacc;lcbacc;cbacc]
        println("$lpacc $majacc $lcbacc $cbacc")
    end

    # Save the output for this w
    matwrite("Hypergraph_Output/Hypergraph_all_w_$w.mat",
    Dict("majs" => majs, "lps" => lps, "allaccs"=>allaccs, "lcbs"=>lcbs, "cbs"=>cbs))

end

## Save Overall Output
matwrite("Hypergraph_Output/Hypergraph_all_vary_w_results.mat",
Dict("n" => n, "Kbound" => Kbound, "Lbound" => Lbound, "q" => q, "p" => p,
"ws" => collect(ws),"maj_accs" => maj_accs, "lp_accs" => lp_accs,
"lp_run" => lp_run, "lcb_accs"=>lcb_accs,
"cb_accs"=>cb_accs, "lcb_run"=>lcb_run, "cb_run"=> cb_run))
