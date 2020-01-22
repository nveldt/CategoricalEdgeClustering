include("../../src/EdgeCatClusAlgs.jl")
using Clustering
using MAT

println("\n")

# Experiment in which we vary the number of nodes

n = 1000
p =  0.05
q = 0.01
Ls = 5:5:50
w = 0.2
numLs = length(Ls)
numtimes = 5

## Store the color identification accuracy for each method
maj_accs = zeros(numws,numtimes)
lp_accs = zeros(numws,numtimes)
cballs_accs = zeros(numws,numtimes)
lcballs_accs = zeros(numws,numtimes)

## Store runtimes (runtime for majority vote is negligible)
lp_run = zeros(numws,numtimes)
cballs_run = zeros(numws,numtimes)
lcballs_run = zeros(numws,numtimes)


for ii = 1:numLs
    Lbound = Ls[ii]
    Kbound = Lbound

    # Save the output of each method for this value of w
    lps = zeros(n,numtimes)
    majs = zeros(n,numtimes)
    cballs = zeros(n,numtimes)
    lcballs = zeros(n,numtimes)

    for times = 1:numtimes

        EdgeList, EdgeColors, node2cluster, node2color = GenerateCCC_model(n,Lbound,Kbound,p,q,w)

        matwrite("Graphs/Synthetic_VaryKL_Graph_$times"*"_L_$Lbound.mat",
        Dict("EdgeList" => EdgeList, "EdgeColors" => EdgeColors,
        "node2cluster" => node2cluster, "node2color" => node2color))

        # Get actual number of colors for this problem
        k = maximum(EdgeColors)

        ## LP-based algorithms
        start = time()
        LPval, X, runtime = EdgeCatClusLP_Graph(EdgeList,EdgeColors,n)
        lp_run[ii,times] = time()-start

        lp_c, RoundScore, RoundRatio = SimpleRound(EdgeList,EdgeColors,X,LPval)
        lp_c = round.(Int64,lp_c)
        lps[:,times] = lp_c
        lp_accs[ii,times] = 1-count(!iszero,lp_c-node2color)/n

        ## Majority Vote Algorithm
        maj_c = MajorityVote(EdgeList,EdgeColors,n,k)
        majs[:,times] = maj_c
        maj_accs[ii,times] = 1-count(!iszero,maj_c-node2color)/n

        # Chromatic Balls and Lazy-Chromatic-Balls
        start = time()
        cballs_cluster, cballs_color = ChromaticBalls(EdgeList,EdgeColors,n)
        cballs_run[ii,times] = time()-start
        cballs[:,times] = cballs_color
        cballs_accs[ii,times] = 1-count(!iszero,cballs_color-node2color)/n

        start = time()
        lcballs_cluster, lcballs_color = LazyChromaticBalls(EdgeList,EdgeColors,n)
        lcballs_run[ii,times] = time()-start
        lcballs[:,times] = lcballs_color
        lcballs_accs[ii,times] = 1-count(!iszero,lcballs_color-node2color)/n

    end

    # Save the output for this L
    matwrite("Graph_Output/Synthetic_VaryKL_Graph_L_$Lbound.mat",Dict("lcballs" => lcballs,
    "cballs" => cballs, "majs" => majs, "lps" => lps, "w" => w))

end

## Save Overall Output
matwrite("Graph_Output/Synthetic_VaryKL_Graph_Results.mat",Dict("Ls" => collect(Ls),
"lcballs_accs" => lcballs_accs, "cballs_run" => cballs_run,
"cballs_accs" => cballs_accs, "maj_accs" => maj_accs,"lp_accs" => lp_accs)
