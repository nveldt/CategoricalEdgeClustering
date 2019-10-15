# Get the list of nodes IDs that are actually present in the hyperedes
filename = "walmart-all/walmart-hyperedges.tsv"
AllNodes = Set{Int64}()
AllColors = Set{Int64}()
f = open(filename)
lines = readlines(f)[2:end]
triptype = Set{Int64}()
close(f)
for i = 1:length(lines)
    a = split(lines[i])
    tt = parse(Int64,a[3])
    push!(triptype,tt)

    nodes = split(a[1],",")
    for j = nodes
        push!(AllNodes,parse(Int64,j))
    end
end

new2old = sort(unique(AllNodes))
m = length(lines)
n = length(new2old)
N = maximum(new2old)
M = length(unique(triptype))
println("$m hyperedges, $M colors, $n unique nodes, $N largest node name")

# Create a map from old indices, to a new set of indices from 1 to n = number of unique nodes
old2new = Dict()
for t = 1:length(new2old)
    old2new[new2old[t]] =  t
end

##
# Now rename all the hyperedges
#
EdgeList = Vector{Vector{Int64}}()
EdgeSizes = Vector{Int64}()
for i = 1:length(lines)
    a = split(lines[i])
    nodes = split(a[1],",")
    EdgeVec = Vector{Int64}()
    for j = nodes
        J = parse(Int64,j)
        Jnew = old2new[J]
        push!(EdgeVec,Jnew)
    end
    push!(EdgeList,EdgeVec)
    push!(EdgeSizes,length(EdgeVec))
end

## Take a look at edge sizes
using Plots
histogram( EdgeSizes, bins = 25)

## Get Node Labels for the set of unique nodes that show up in hyperedges
filename = "walmart-all/walmart-nodes.tsv"
f = open(filename)
lines = readlines(f)[2:end]
close(f)

oldID2label = zeros(Int64,N)
for i = 1:length(lines)
    pair = parse.(Int64,split(lines[i],"\t"))
    @assert(pair[1] == i)
    oldID2label[i] = pair[2]
end

NodeLabels = oldID2label[new2old]

## Look at node label histogram
histogram( NodeLabels, bins= 68)

## Get Node Label List
filename = "walmart-all/node_cats.txt"
f = open(filename)
lines = readlines(f)[1:end]
close(f)

Categories = Vector{String}()
for i = 1:length(lines)
    pair = split(lines[i],":")
    push!(Categories,pair[2])
end

## Check number of nodes in each category
numtype = zeros(68)
for j = 1:68
    inds = findall(x->x==j,NodeLabels)
    num = length(inds)
    numtype[j] = num

    println(Categories[j]*" ($j): $num products")
end

## The original 68 categories are too broad. We will refine them.
filename = "walmart-all/Ten_Departments.txt"
f = open(filename)
lines = readlines(f)[1:end]
close(f)

old2new_label = zeros(Int64,68)
Departments = Vector{String}()
for i = 1:length(lines)
    pair = split(lines[i],":")
    oldset = parse.(Int64,split(pair[2],","))
    for j in oldset
        old2new_label[j] = i
    end
    push!(Departments,pair[1])
end

# Now put each node in its department category
NodeDepartment = zeros(Int64,n)
for i = 1:n
    NodeDepartment[i] = old2new_label[NodeLabels[i]]
end

# Check the histogram of the new department-based labels
numdeps = maximum(NodeDepartment)
histogram(NodeDepartment,bins = numdeps)

##
numtype = zeros(numdeps)
for j = 1:numdeps
    inds = findall(x->x==j,NodeDepartment)
    num = length(inds)
    numtype[j] = num
    println(Departments[j]*" ($j): $num products")
end

# Get several different subsets of the data for analysis
clothes_inds = findall(x->x==1,NodeDepartment)
clothes_labels = NodeLabels[clothes_inds]
food_inds = findall(x->x==7,NodeDepartment)
food_labels = NodeLabels[food_inds]
home_inds = findall(x->x==3,NodeDepartment)
home_labels = NodeLabels[home_inds]
EO_inds = findall(x->x==2,NodeDepartment)
EO_labels = NodeLabels[EO_inds]
health = findall(x->x==8,NodeDepartment)
keep = findall(x->x<numdeps,NodeDepartment)
KeepLabels = NodeDepartment[keep]

save("../JLD_Files/Walmart_NodeLabels.jld","EdgeList", EdgeList,
"Categories", Categories, "Departments", Departments, "n", n,
"NodeLabels",NodeLabels, "NodeDepartment",NodeDepartment)
