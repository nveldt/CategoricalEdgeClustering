using MAT
using JLD
data = load("../data/JLD_Files/Cooking.jld")
EdgeColors_1 = data["EdgeColors"]
EdgeList_1 = data["EdgeList"]
Cuisine2Label_1 = data["Cuisine2Label"]
Cuisines_1 = data["Cuisines"]
Ingredient2num_1 = data["Ingredient2num"]
Ingredients_1 = data["Ingredients"]
n_1 = length(Ingredients_1)

## Choose a value of beta and a cuisine number (from 1 to 20)
# This will extract all of the differences in recipes and ingredients between
# LP-round and Majority Vote, saving the output to a file.
beta = 70
cuis = 15
tag = ""
printatmost = 30    # maximum number of ingredients you want to print to the file
output_mat = "Output/Beta_$beta"*"_stats.mat"
output_txt = "Cluster_Summary/Beta_$beta"*"_"*Cuisines_1[cuis]*"_"*tag*".txt"
EdgeList, EdgeColors, n, newlist1, old2new, NewListInds = RemoveByTotalDegree(EdgeList_1,EdgeColors_1,n_1,beta)
Ingredients = Ingredients_1[newlist1]
EdgeList, old2new, newlist2 = RenameNodes(EdgeList)
Ingredients = Ingredients[newlist2]
n = length(Ingredients)
M = length(EdgeColors)
k = maximum(EdgeColors)
Cdeg = get_color_degree(EdgeList,EdgeColors,n)

EdgeMap = Dict()
for j = 1:length(NewListInds)
    EdgeMap[NewListInds[j]] = j
end

mat = matread(output_mat)
lp_c = mat["lp_c"]
maj_c = mat["maj_c"]
lp_cuisine = mat["lp_cuisine"]
maj_cuisine = mat["maj_cuisine"]

# We'll save all the output data
lp_cuisine = zeros(9,20)
maj_cuisine = zeros(9,20)
LP_ing_list = zeros(n,20)
Maj_ing_list = zeros(n,20)
Shared_ing_list = zeros(n,20)
Cuisine_Ings = zeros(20)
Cuisine_Recipes = zeros(20)
Shared_Active_Counts = zeros(20)
Shared_Used_Counts = zeros(20)
open(output_txt, "w") do f
        write(f,"Beta = $beta: $n nodes and $M hyperedges\n")
end
println("Beta = $beta: $n nodes and $M hyperedges")

# Wasted ingredients are those nodes which could be thrown out without
# changing the number of satisfied recipes. Count these for each cuisine
# and for each method
LP_wasted = zeros(20)
Maj_wasted = zeros(20)
Shared_wasted = zeros(20)
All_shared = zeros(20)

# Record the number of ingredients and recipes total for this cuisine
cuis_deg = Cdeg[cuis,:]
cuis_ingredients = findall(!iszero,cuis_deg)
cuis_recipes = findall(x->x==cuis, EdgeColors)

# Get indices for the ingredients each algorithm placed in this cuisine
lp_cuis_inds = findall(x->x==cuis,lp_c)
maj_cuis_inds = findall(x->x==cuis,maj_c)
shared_cuis_inds = intersect(lp_cuis_inds,maj_cuis_inds)

lp_unique = setdiff(lp_cuis_inds,shared_cuis_inds)
maj_unique = setdiff(maj_cuis_inds,shared_cuis_inds)

LP_Used = Set{Int64}()
Maj_Used = Set{Int64}()
Shared_Used = Set{Int64}()

lp_recipes = Set{Int64}()
maj_recipes = Set{Int64}()
shared_recipes = Set{Int64}()

for i = 1:length(cuis_recipes)

    # Get i-th recipe from this cuisine. Check total number of ingredients
    r_num = cuis_recipes[i]
    recipe = EdgeList[r_num]
    ing_count = length(recipe)
    lp_has = intersect(lp_cuis_inds,recipe)
    maj_has = intersect(maj_cuis_inds,recipe)
    shared_has = intersect(shared_cuis_inds,recipe)

    if length(lp_has) == ing_count       # LP method has all recipes
        push!(lp_recipes,r_num)
        for ingred in lp_has
            push!(LP_Used,ingred)
        end
    end
    if length(maj_has) == ing_count      # Do the same for Majority Vote
        push!(maj_recipes,r_num)
        for ingred in maj_has
            push!(Maj_Used,ingred)
        end
    end
    if length(shared_has) == ing_count   # Do the same for the intersection
        push!(shared_recipes,r_num)
        for ingred in shared_has
            push!(Shared_Used,ingred)
        end
    end
end

LP_unique_used = setdiff(LP_Used,Shared_Used)
Maj_unique_used = setdiff(Maj_Used,Shared_Used)
LP_unique_recipes = setdiff(lp_recipes,shared_recipes)
Maj_unique_recipes = setdiff(maj_recipes,shared_recipes)

lp_used = length(unique(LP_Used))
maj_used = length(unique(Maj_Used))
shared_used = length(unique(Shared_Used))
lunique = length(lp_unique)
luu = length(LP_unique_used)

additional = setdiff(LP_unique_used,lp_unique)
more = length(additional)

open(output_txt, "a") do f
    write(f,"There are $lunique ingredients that LP includes and MV does not\n")
end

if lunique < printatmost

    lpos = collect(lp_unique)
    open(output_txt, "a") do f
        write(f,"\nThe LP includes $lunique ingredients that Majority does not:\n")
    end
    for jj = 1:lunique
        nodeid = lpos[jj]
        ingre = Ingredients[nodeid]
        open(output_txt, "a") do f
            write(f,ingre*", ")
        end
    end
end

# Do the same for used ingredients
if more < printatmost
    lpos = collect(additional)
    open(output_txt, "a") do f
        write(f,"\nThe LP uses $more additional ingredients that Majority does not use:\n")
    end
    for jj = 1:more
        nodeid = lpos[jj]
        ingre = Ingredients[nodeid]
        open(output_txt, "a") do f
            write(f,ingre*", ")
        end
    end
end

## Recipes that are unique to LP-round
open(output_txt, "a") do f
    write(f,"\nLP can make the following recipes that Majority Vote cannot.\n")
end

lur = unique(LP_unique_recipes)
for j = 1:length(lur)
    r = lur[j]

    recipe = EdgeList[r]
    # println("\nRecipe $j")
    # for ingredient in recipe
    #     println(Ingredients[ingredient])
    # end

    recipe = EdgeList_1[NewListInds[r]]
    println("\nRecipe $j with common ingredients:")
    open(output_txt, "a") do f
        write(f,"\n\nRecipe $j with common ingredients:\n ")
    end
    for ingredient in recipe
        println(Ingredients_1[ingredient]*" $ingredient")

        open(output_txt, "a") do f
            write(f,Ingredients_1[ingredient]*", ")
        end
    end
    # open(output_txt, "a") do f
    #     write(f,"\$\\}\$ ")
    # end

end
