include("../../src/EdgeCatClusAlgs.jl")
using Clustering
using MAT

println("\n")

n = 1000
w = 0.2
numtimes = 5
p =  0.005
q = 0.0001
Ls = 5:5:50
numLs = length(Ls)
## Store the ARI scores and the color identification accuracy for each method
maj_accs = zeros(numLs,numtimes)
lp_accs = zeros(numLs,numtimes)
lcb_accs = zeros(numLs,numtimes)
cb_accs = zeros(numLs,numtimes)

## Store runtimes, just in case we end up caring
lp_run = zeros(numLs,numtimes)
lcb_run = zeros(numLs,numtimes)
cb_run = zeros(numLs,numtimes)

for ii = 1:numLs
    L = Ls[ii]

    # Save the output of each method for this value of L
    lps = zeros(n,numtimes)
    majs = zeros(n,numtimes)
    lcbs = zeros(n,numtimes)
    cbs = zeros(n,numtimes)

    allaccs = zeros(4,numtimes)


    for times = 1:numtimes
    # times = 1

        EdgeList, EdgeColors, node2cluster, node2color = GenerateCCC_model3(n,L,L,p,q,w,true)

        # Get actual number of colors for this problem
        k = maximum(EdgeColors)
        M = length(EdgeColors)

        println("$k $M")

        ## LP-based algorithms
        start = time()
        LPval, X, runtime = EdgeCatClusGeneral(EdgeList,EdgeColors,n,false,1)
        lp_run[ii,times] = time()-start

        matwrite("Graphs/Hyper3unif_L_$L"*"_$times.mat",
        Dict("EdgeList" => EdgeList, "EdgeColors" => EdgeColors,
        "node2cluster" => node2cluster, "node2color" => node2color,
        "LPval" => LPval))

        lp_c, RoundScore, RoundRatio = SimpleRound(EdgeList,EdgeColors,X,LPval)
        lp_c = round.(Int64,lp_c)
        lps[:,times] = lp_c
        lpacc = 1-count(!iszero,lp_c-node2color)/n
        lp_accs[ii,times] = lpacc

        ## Majority Vote
        maj_c = MajorityVote(EdgeList,EdgeColors,n,k)
        majs[:,times] = maj_c
        majacc = 1-count(!iszero,maj_c-node2color)/n
        maj_accs[ii,times] = majacc

        NewList, NewLabels = Hypergraph2Graph(EdgeList,EdgeColors)
        println("Running Chromatic Balls")
        start = time()
        cb_cluster, cb_c = ChromaticBalls(NewList,NewLabels,n)
        cb_run[ii,times] = time()-start
        cbs[:,times] = cb_c
        cbacc = 1-count(!iszero,cb_c-node2color)/n
        cb_accs[ii,times] = cbacc

        println("Running Lazy Chromatic Balls")
        start = time()
        lcb_cluster, lcb_c = LazyChromaticBalls(NewList,NewLabels,n)
        lcb_run[ii,times] = time()-start
        lcbs[:,times] = lcb_c
        lcbacc =  1-count(!iszero,lcb_c-node2color)/n
        lcb_accs[ii,times] = lcbacc
        allaccs[:,times] = [lpacc;majacc;lcbacc;cbacc]
        println("$lpacc $majacc $lcbacc $cbacc")
    end

    # Save the output for this w
    matwrite("Output/Hyper3unif_results_L_$L.mat",
    Dict("majs" => majs, "allaccs"=>allaccs, "lcbs" =>lcbs, "cbs" =>cbs, "lps" => lps, "w" => w))

end

## Save Overall Output
matwrite("Output/Hyper3unif_vary_kl_results.mat",
Dict("n" => n, "w"=>w, "q" => q, "p" => p,
"Ls" => collect(Ls),"maj_accs" => maj_accs, "lp_accs" => lp_accs,
"lp_run" => lp_run, "lcb_run"=>lcb_run, "lcb_accs"=>lcb_accs,
"cb_run"=>cb_run, "cb_accs"=>cb_accs))
