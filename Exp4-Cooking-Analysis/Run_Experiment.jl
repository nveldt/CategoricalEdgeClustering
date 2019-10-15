include("../../src/EdgeCatClusAlgs.jl")
using MAT
using JLD
data = load("../../data/JLD_Files/Cooking.jld")
EdgeColors_1 = data["EdgeColors"]
EdgeList_1 = data["EdgeList"]
Cuisine2Label_1 = data["Cuisine2Label"]
Cuisines_1 = data["Cuisines"]
Ingredient2num_1 = data["Ingredient2num"]
Ingredients_1 = data["Ingredients"]
n_1 = length(Ingredients_1)

##
# Remove ingredients that show up in too many dishes
betas = 0:10:200

numtimes = length(betas)
LP_stats = zeros(4,numtimes)
Maj_stats = zeros(3,numtimes)
lp_total_wasted = zeros(numtimes)
maj_total_wasted = zeros(numtimes)
Problem_stats = zeros(2,numtimes)
Total_shared = zeros(numtimes)
printatmost = 30
output_final = "Output/All_stats.mat"

for ii = 1:length(betas)
   beta = betas[ii]

   output_txt = "Output/Beta_$beta"*"_stats.txt"
   output_mat = "Output/Beta_$beta"*"_stats.mat"

   # Remove all ingredients that are in more than beta dishes that ARE NOT in
   # the cuisine that the ingredient appears most frequently in.
   EdgeList, EdgeColors, n, newlist1, old2new = RemoveByTotalDegree(EdgeList_1,EdgeColors_1,n_1,beta)
   Ingredients = Ingredients_1[newlist1]

   # Relabel nodes, since after refinement some nodes have degree zero
   EdgeList, old2new, newlist2 = RenameNodes(EdgeList)
   Ingredients = Ingredients[newlist2]
   ingredient_subset = newlist1[newlist2]       # gives indices of ORIGINAL data that is considered here
   @assert(Ingredients==Ingredients_1[ingredient_subset])
   n = length(Ingredients)
   M = length(EdgeColors)
   Problem_stats[:,ii] = [n,M]
   k = maximum(EdgeColors)
   Cdeg = get_color_degree(EdgeList,EdgeColors,n)

   # Run the LP algorithm
   start = time()
   LPval, X, runtime = EdgeCatClusGeneral(EdgeList,EdgeColors,n,false,1)
   lp_time = round(time()-start,digits=2)
   C = rowmin(X)
   lp_c = C[:,2]
   lp_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,lp_c)
   lp_ratio = lp_mistakes/LPval
   lp_sat = 1 - lp_mistakes/M

   # Run the Majority Vote Algorithm
   maj_c = NaiveLabel(EdgeList,EdgeColors,n,k,1)
   maj_mistakes = EdgeCatClusObj(EdgeList,EdgeColors,maj_c)
   maj_ratio = maj_mistakes/LPval
   maj_sat = 1 - maj_mistakes/M

   # Save the output
   LP_stats[:,ii] = [lp_sat; lp_ratio; lp_mistakes; lp_time]
   Maj_stats[:,ii] = [maj_sat; maj_ratio; maj_mistakes]

   # Now do a cuisine-by-cuisine analysis

   # We'll save all the output data
   lp_cuisine = zeros(9,20)
   maj_cuisine = zeros(9,20)
   LP_ing_list = zeros(n,20)
   Maj_ing_list = zeros(n,20)
   Shared_ing_list = zeros(n,20)
   Cuisine_Ings = zeros(20)
   Cuisine_Recipes = zeros(20)

   # How many "active" (i.e., partaking in some recipe of this cuisine) and
   # "used" (i.e., contributing to the edge satisfaction score) ingredients
   # are there in the overlap between LP and Majority, for each cuisine?
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

   for cuis = 1:20

       open(output_txt, "a") do f
               write(f, "\n\n-----------------------------------------------------------------------\n")
               write(f,Cuisines_1[cuis]*"\n I = ingredients; A = Active; U = Used\n")
               write(f, "-------------------------------------------------------------------------------------------------\n")
               write(f,"\t\t (I,A,U) \t # Recipes \t # Waste \t Ing Diff \t Active Diff \t Used Diff \t Recipes Diff \n")
       end
       println("---------------------------------------------------------------------------------\n")
       println("\n"*Cuisines_1[cuis]*"\n I = ingredients; A = Active; U = Used ")
       println("-----------------------------------------------------------------------------------------------")
       println(" \t\t (I,A,U) \t # Recipes \t # Waste \t Ing Diff \t Active Diff \t Used Diff \t Recipes Diff")

       # Record the number of ingredients and recipes total for this cuisine
       cuis_deg = Cdeg[cuis,:]
       cuis_ingredients = findall(!iszero,cuis_deg)
       cuis_recipes = findall(x->x==cuis, EdgeColors)
       ings = length(cuis_ingredients)
       recs = length(cuis_recipes)
       Cuisine_Ings[cuis] = ings                # Total number of ingredients
       Cuisine_Recipes[cuis] = recs             # Total number of recipes

       # Get indices for the ingredients each algorithm placed in this cuisine
       lp_cuis_inds = findall(x->x==cuis,lp_c)
       maj_cuis_inds = findall(x->x==cuis,maj_c)
       shared_cuis_inds = intersect(lp_cuis_inds,maj_cuis_inds)
       shared_total = length(shared_cuis_inds)

       # Save those for later
       LP_ing_list[lp_cuis_inds,cuis] .= 1
       Maj_ing_list[maj_cuis_inds,cuis] .= 1
       Shared_ing_list[shared_cuis_inds,cuis] .= 1

       # Count the number of ingredients each method associates with this cuisine
       lp_ing_num = length(lp_cuis_inds)
       maj_ing_num = length(maj_cuis_inds)
       shared_ing_num = length(shared_cuis_inds)
       All_shared[cuis] = shared_ing_num

       lp_diff = lp_ing_num - shared_ing_num
       maj_diff = maj_ing_num - shared_ing_num

       # Active ingredient = ingredient that is in SOME recipe of this cuisine
       LP_Active = Set{Int64}()
       Maj_Active = Set{Int64}()
       Shared_Active = Set{Int64}()

       # Used ingredient = ingredient that is in a recipe that is "satisfied"
       # by this method. I.e., if we remove the ingredient, the edge satisfaction
       # score of the method goes down.
       LP_Used = Set{Int64}()
       Maj_Used = Set{Int64}()
       Shared_Used = Set{Int64}()

       # Count how many recipes each algorithm can make with their ingredients.
       # This will require iterating through the list of recipes of this cuisine.
       lp_rec = 0
       maj_rec = 0
       share_rec = 0
       for i = 1:length(cuis_recipes)

           # Get i-th recipe from this cuisine. Check total number of ingredients
           r_num = cuis_recipes[i]
           recipe = EdgeList[r_num]
           ing_count = length(recipe)

           # How many of these ingredients does the lp have?
           lp_has = intersect(lp_cuis_inds,recipe)
           for ingredient in lp_has
               push!(LP_Active,ingredient)
           end

           # How many of these ingredients does maj-vote have?
           maj_has = intersect(maj_cuis_inds,recipe)
           for ingredient in maj_has
               push!(Maj_Active,ingredient)
           end

           # How many of these ingredients are shared by both?
           shared_has = intersect(shared_cuis_inds,recipe)
           for ingredient in shared_has
               push!(Shared_Active,ingredient)
           end

           # If Method X put all the recipes ingredients in the right cluster,
           # update the number of satisfied recipes for Method X, and update
           # the list of USED ingredients.
           if length(lp_has) == ing_count       # LP method has all recipes
               lp_rec += 1                      # Increase recipe satisfaction count
               for ingred in lp_has             # Put those ingredients in the "Used" category
                   push!(LP_Used,ingred)
               end
           end
           if length(maj_has) == ing_count      # Do the same for Majority Vot
               maj_rec += 1
               for ingred in maj_has
                   push!(Maj_Used,ingred)
               end
           end
           if length(shared_has) == ing_count   # Do the same for the intersection
               share_rec +=1
               for ingred in shared_has
                   push!(Shared_Used,ingred)
               end
           end
       end

       # Get the set of ingredients BEYOND the intersection that each method has.
       # Do this for both "active" and "used" ingredients.
       LP_over_Shared_active = setdiff(LP_Active,Shared_Active)
       Maj_over_Shared_active = setdiff(Maj_Active,Shared_Active)

       LP_over_Shared_used = setdiff(LP_Used,Shared_Used)
       Maj_over_Shared_used = setdiff(Maj_Used,Shared_Used)

       lp_active = length(unique(LP_Active))
       maj_active = length(unique(Maj_Active))
       shared_active = length(unique(Shared_Active))

       lp_used = length(unique(LP_Used))
       maj_used = length(unique(Maj_Used))
       shared_used = length(unique(Shared_Used))

       Shared_Active_Counts[cuis] = shared_active
       Shared_Used_Counts[cuis] = shared_used

       lp_over_shared_active = length(LP_over_Shared_active)
       maj_over_shared_active = length(Maj_over_Shared_active)
       lp_over_shared_used = length(LP_over_Shared_used)
       maj_over_shared_used = length(Maj_over_Shared_used)

       # Count how many more recipes beyond the intersection that Maj and LP can make
       count_lp = lp_rec-share_rec
       count_maj = maj_rec-share_rec

       # Count total number of wasted ingredients
       lpwaste = lp_ing_num-lp_used
       majwaste = maj_ing_num-maj_used
       sharedwaste = shared_ing_num-shared_used
       LP_wasted[cuis] = lpwaste
       Maj_wasted[cuis] = majwaste
       Shared_wasted[cuis] = sharedwaste

       # Store it all
       lp_cuisine[:,cuis] = [lp_ing_num; lp_active; lp_used; lp_rec; lpwaste; lp_diff; lp_over_shared_active; lp_over_shared_used; count_lp]
       maj_cuisine[:,cuis] = [maj_ing_num; maj_active; maj_used; maj_rec; majwaste; maj_diff; maj_over_shared_active; maj_over_shared_used; count_maj]

       # Print it all
       open(output_txt, "a") do f
           write(f,"Every:\t\t ($ings,-,-,-) \t\t $recs\n")
           write(f,"LPalg: \t\t ($lp_ing_num, $lp_active, $lp_used) \t\t $lp_rec \t\t $lpwaste \t\t $lp_diff \t\t $lp_over_shared_active \t\t $lp_over_shared_used \t\t $count_lp \n")
           write(f,"Major:\t\t ($maj_ing_num,$maj_active,$maj_used) \t\t $maj_rec \t\t $majwaste \t\t $maj_diff \t\t $maj_over_shared_active \t\t $maj_over_shared_used \t\t $count_maj \n")
           write(f,"Share:\t\t ($shared_total, $shared_active,$shared_used) \t\t $share_rec \t\t $sharedwaste\n")
       end
       println("Every:\t\t ($ings,-,-,-) \t\t $recs")
       println("LPalg: \t\t ($lp_ing_num,$lp_active,$lp_used) \t\t $lp_rec \t\t $lpwaste  \t\t $lp_diff \t\t $lp_over_shared_active \t\t $lp_over_shared_used \t\t $count_lp")
       println("Major:\t\t ($maj_ing_num,$maj_active,$maj_used) \t\t $maj_rec \t\t $majwaste \t\t $maj_diff \t\t $maj_over_shared_active \t\t $maj_over_shared_used \t\t $count_maj")
       println("Share:\t\t ($shared_total,$shared_active,$shared_used) \t\t $share_rec \t\t $sharedwaste")

       if lp_over_shared_active < printatmost
           # Then print them in the output text file
           lpos = collect(LP_over_Shared_active)
           open(output_txt, "a") do f
               write(f,"The LP includes $lp_over_shared_active active ingredients that Majority does not:\n")
           end
           for jj = 1:lp_over_shared_active
               nodeid = lpos[jj]
               ingre = Ingredients[nodeid]
               open(output_txt, "a") do f
                   write(f,ingre*"\n")
               end
           end
       end
       # Do the same for used ingredients
       if lp_over_shared_used < printatmost
           lpos = collect(LP_over_Shared_used)
           open(output_txt, "a") do f
               write(f,"The LP includes $lp_over_shared_used used ingredients that Majority does not.\n")
           end
           for jj = 1:lp_over_shared_used
               nodeid = lpos[jj]
               ingre = Ingredients[nodeid]
               open(output_txt, "a") do f
                   write(f,ingre*"\n")
               end
           end
       end

   end
   Total_shared[ii] = sum(All_shared)
   lp_total_wasted[ii] = sum(LP_wasted)
   maj_total_wasted[ii] = sum(Maj_wasted)

   matwrite(output_mat, Dict("EdgeList"=>EdgeList,
   "EdgeColors" => EdgeColors, "n" => n, "ingredient_subset" => ingredient_subset,
   "lp_c" => lp_c, "LPval" => LPval, "maj_c" => maj_c, "lp_cuisine" => lp_cuisine,
   "maj_cuisine" => maj_cuisine, "LP_ing_list" => LP_ing_list,
   "Maj_ing_list" => Maj_ing_list, "Shared_ing_list" => Shared_ing_list,
   "Shared_Active_Counts" => Shared_Active_Counts,
   "Shared_Used_Counts" => Shared_Used_Counts))


end

## Save the overall stats
matwrite(output_final,Dict("LP_stats" => LP_stats, "Maj_stats" => Maj_stats,
"lp_wasted" => lp_total_wasted, "maj_wasted" => maj_total_wasted,
"Problem_stats" => Problem_stats, "Total_shared" => Total_shared))
