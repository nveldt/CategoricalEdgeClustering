using JSON
using JLD
cooking = JSON.parsefile("train.json")

# Extract the Ingredients
Ingredients = Set{String}()
Cuisines = Set{String}()
for entry in cooking
    for i in entry["ingredients"]
        push!(Ingredients,i)
    end
    push!(Cuisines,entry["cuisine"])
end

Ingredients = collect(Ingredients)
Cuisines = collect(Cuisines)
EdgeList = Vector{Vector{Int64}}()
EdgeColors = Vector{Int64}()
# Go from ingredient name to node number
CuNum = Dict()
IngNum = Dict()
for i = 1:length(Ingredients)
    IngNum[Ingredients[i]] = i

    if i <= length(Cuisines)
        CuNum[Cuisines[i]] = i
    end
end

for recipe in cooking

    # Put the nodes in an edgelist
    edgevec = Vector{Int64}()
    for i in recipe["ingredients"]
        i_num = IngNum[i]
        push!(edgevec,i_num)
    end
    push!(EdgeList,edgevec)
    push!(EdgeColors,CuNum[recipe["cuisine"]])
end
Cuisine2Label = CuNum
Ingredient2num= IngNum
n = length(Ingredients)

save("../JLD_Files/Cooking.jld","EdgeList", EdgeList, "EdgeColors", EdgeColors, "Cuisine2Label", Cuisine2Label, "Ingredient2num", Ingredient2num,"Ingredients",Ingredients, "Cuisines", Cuisines, "n", n)
