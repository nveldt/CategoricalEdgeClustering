include("../../src/EdgeCatClusAlgs.jl")
using Clustering
using MAT

println("\n")
mat = matread("../../data/Email-eu-core/Email_Graph.mat")
A = mat["A"]
A = triu(A)
c = vec(round.(Int64,mat["truth"]))
n = size(A,1)

I,J,V = findnz(A)
K = maximum(c)
k = K

ws = 0.0:.05:0.75
numws = length(ws)
numtimes = 1
node2color = c
node2cluster = c

## Store the label identification accuracy for each method
maj_accs = zeros(numws,numtimes)
lp_accs = zeros(numws,numtimes)
cballs_accs = zeros(numws,numtimes)
lcballs_accs = zeros(numws,numtimes)

## Store runtimes
lp_run = zeros(numws,numtimes)
cut_run = zeros(numws,numtimes)
cballs_run = zeros(numws,numtimes)
lcballs_run = zeros(numws,numtimes)

for ii = 1:numws
    w = ws[ii]

    # Save the output of each method for this value of w
    lps = zeros(n,numtimes)
    majs = zeros(n,numtimes)
    cballs = zeros(n,numtimes)
    lcballs = zeros(n,numtimes)

    for times = 1:numtimes

        # Assign synthetic edge labels
        for t = 1:length(I)
            V[t] = rand(1:K)
            # Assign the right color edge only with probability beta, if they are in the same cluster
            if c[I[t]] == c[J[t]] && rand(1)[1] < 1-w
                V[t] = round(Int64,c[J[t]])
            end
        end

        ##
        EdgeList = [I J]
        EdgeColors = round.(Int64,V)

        ## LP-based algorithms
        start = time()
        LPval, X, runtime = EdgeCatClusLP_Graph(EdgeList,EdgeColors,n,0)
        lp_run[ii,times] = time()-start

        lp_c, RoundScore, RoundRatio = SimpleRound(EdgeList,EdgeColors,X,LPval)
        lp_c = round.(Int64,lp_c)
        lps[:,times] = lp_c
        lp_accs[ii,times] = 1-count(!iszero,lp_c-node2color)/n

        ## Majority Vote
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

        println("$k \t\t $w")

    end

    # Save the output for this w
    matwrite("Output/Email_w_$w.mat",
    Dict("lcballs" => lcballs, "cballs" => cballs,
    "majs" => majs "lps" => lps, "w" => w))

end

## Save Overall Output
matwrite("Output/Email_Results.mat",
Dict("ws" => collect(ws), "maj_accs" => maj_accs, "lp_accs" => lp_accs,
"lp_run" => lp_run, "lcballs_run" => lcballs_run,"lcballs_accs" => lcballs_accs,
 "cballs_run" => cballs_run,"cballs_accs" => cballs_accs))
