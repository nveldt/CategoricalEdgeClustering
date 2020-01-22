include("../../src/EdgeCatClusAlgs.jl")
using MAT

# This is designed to be a parameter regime taken from Bonchi et al.
# (Chromatic Correlation Clustering, Bonchi et al. 2015). The number of
# clusters (K) is much larger than the number of labels (L).

println("\n")

n = 1000
p =  0.5
q = 0.03
L = 20
Ks = 50:25:200
w = 0.5
numKs = length(Ks)
numtimes = 5
@show numtimes

## Store the ARI scores and the label identification accuracy for each method
maj_aris = zeros(numKs,numtimes)
maj_accs = zeros(numKs,numtimes)
lp_aris = zeros(numKs,numtimes)
lp_accs = zeros(numKs,numtimes)
cballs_aris = zeros(numKs,numtimes)
cballs_accs = zeros(numKs,numtimes)
lcballs_aris = zeros(numKs,numtimes)
lcballs_accs = zeros(numKs,numtimes)

## Store runtimes
lp_run = zeros(numKs,numtimes)
cballs_run = zeros(numKs,numtimes)
lcballs_run = zeros(numKs,numtimes)

for ii = 1:numKs
    Lbound = L
    Kbound = Ks[ii]

    # Save the output of each method for this value of w
    lps = zeros(n,numtimes)
    majs = zeros(n,numtimes)

    # Labelings and clusterings are distinct for LCB and CB
    cballs_lab = zeros(n,numtimes)
    lcballs_lab = zeros(n,numtimes)
    cballs_clus = zeros(n,numtimes)
    lcballs_clus = zeros(n,numtimes)

    for times = 1:numtimes

        EdgeList, EdgeColors, node2cluster, node2color = GenerateCCC_model(n,Lbound,Kbound,p,q,w)

        matwrite("Graphs/Synthetic_Large_K_Graph_$times"*"_K_$Kbound.mat",
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
        lp_aris[ii,times] = ari(lp_c, node2cluster)
        lp_accs[ii,times] = 1-count(!iszero,lp_c-node2color)/n


        ## Majority Vote Algorithm
        maj_c = MajorityVote(EdgeList,EdgeColors,n,k)
        majs[:,times] = maj_c
        maj_aris[ii,times] = ari(maj_c, node2cluster)
        maj_accs[ii,times] = 1-count(!iszero,maj_c-node2color)/n

        # Chromatic Balls and Lazy-Chromatic-Balls return a labeling AND a clustering
        start = time()
        cballs_cluster, cballs_color = ChromaticBalls(EdgeList,EdgeColors,n)
        cballs_run[ii,times] = time()-start
        cballs_clus[:,times] = cballs_cluster
        cballs_lab[:,times] = cballs_color
        cballs_aris[ii,times] = ari(cballs_cluster, node2cluster)
        cballs_accs[ii,times] = 1-count(!iszero,cballs_color-node2color)/n

        start = time()
        lcballs_cluster, lcballs_color = LazyChromaticBalls(EdgeList,EdgeColors,n)
        lcballs_run[ii,times] = time()-start
        lcballs_clus[:,times] = lcballs_cluster
        lcballs_lab[:,times] = lcballs_color
        lcballs_aris[ii,times] = ari(lcballs_cluster, node2cluster)
        lcballs_accs[ii,times] = 1-count(!iszero,lcballs_color-node2color)/n

        println("$k \t\t $Kbound")

    end

    # Save the output for this K
    matwrite("Graph_Output/Synthetic_Large_K_Graph_$Kbound.mat",
    Dict("lcballs_clus" => lcballs_clus, "lcballs_lab" => lcballs_lab,
    "cballs_clus" => cballs_clus, "cballs_lab" => cballs_lab,
    "majs" => majs, "lps" => lps, "w" => w))

end

## Save Overall Output
matwrite("Graph_Output/Synthetic_Large_K_Graph_Results.mat",
Dict("n" => n, "Ks" => collect(Ks), "qout" => q, "p" => p,
"w" => w, "lcballs_run" => lcballs_run, "lcballs_aris" => lcballs_aris,
"lcballs_accs" => lcballs_accs, "cballs_run" => cballs_run,
"cballs_aris" => cballs_aris,"cballs_accs" => cballs_accs,
"maj_aris" => maj_aris, "maj_accs" => maj_accs,
"lp_aris" => lp_aris, "lp_accs" => lp_accs, "lp_run" => lp_run))
