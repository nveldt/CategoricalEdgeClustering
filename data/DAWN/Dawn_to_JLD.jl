filename = "colored-DAWN.txt"
AllNodes = Vector{Int64}()
AllColors = Vector{Int64}()
f = open(filename)
lines = readlines(f)
close(f)
for i = 1:length(lines)
    a = split(lines[i])
    color = parse(Int64,a[2])
    push!(AllColors,color)

    nodes = split(a[1],",")
    for j = nodes
        push!(AllNodes,parse(Int64,j))
    end
end
new2old = sort(unique(AllNodes))
m = length(lines)
n = length(new2old)
N = maximum(new2old)
M = length(unique(AllColors))
println("$m hyperedges, $M colors, $n unique nodes, $N largest node name")

old2new = Dict()
for t = 1:length(new2old)
    old2new[new2old[t]] =  t
end

# Now store the hyperedges under the new names
EdgeList = Vector{Vector{Int64}}()
EdgeColors = zeros(Int64,m)

for i = 1:length(lines)
    a = split(lines[i])
    color = parse(Int64,a[2])
    EdgeColors[i] = color

    nodes = split(a[1],",")
    EdgeVec = Vector{Int64}()
    for j = nodes
        J = parse(Int64,j)
        Jnew = old2new[J]
        push!(EdgeVec,Jnew)
    end
    push!(EdgeList,EdgeVec)
end

# Now read in the node IDs and their names
f = open("colored-DAWN-nodes.txt")
lines = readlines(f)[2:end]
close(f)

NodeNames = Vector{String}(UndefInitializer(),n)
for lind = 1:length(lines)
    l = lines[lind]
    idname = split(l,"\t")
    id = parse(Int64,idname[1])
    name = idname[2]
    if haskey(old2new,id)
        newID = old2new[id]
        NodeNames[newID] = name
    end
end

EdgeList, EdgeColors, n, NodeNames


Labels = ["1: DISCHARGED HOME",
"2: RELEASED TO POLICE/JAIL",
"3: REFERRED TO DETOX/TREATMENT",
"4: ICU/CRITICAL CARE",
"5: SURGERY",
"6: CHEMICAL DEPENDENCY/DETOX, PSYCHIATRIC UNIT",
"7: OTHER INPATIENT UNIT",
"8: TRANSFERRED",
"9: LEFT AGAINST MEDICAL ADVICE",
"10: DIED"]

save("../JLD_Files/DAWN.jld","EdgeList", EdgeList, "EdgeColors", EdgeColors, "NodeNames", NodeNames, "Labels", Labels, "n",n)
