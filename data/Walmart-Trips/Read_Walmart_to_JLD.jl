filename = "colored-walmart.txt"
AllNodes = Set{Int64}()
AllColors = Set{Int64}()
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

##
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

EdgeColors, blank, blank2 = relabel(EdgeColors)
EdgeColors = round.(Int64,EdgeColors)
##

save("../JLD_Files/Walmart-Trips.jld","EdgeList", EdgeList, "EdgeColors", EdgeColors, "n", n)
