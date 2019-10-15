## Load the original dataset
using MAT
using JLD
using Plots
data = load("../data/JLD_Files/Cooking.jld")
EdgeColors_1 = data["EdgeColors"]
EdgeList_1 = data["EdgeList"]
Cuisine2Label_1 = data["Cuisine2Label"]
Cuisines_1 = data["Cuisines"]
Ingredient2num_1 = data["Ingredient2num"]
Ingredients_1 = data["Ingredients"]
n_1 = length(Ingredients_1)

output_final = "Output/All_stats.mat"

## Load the data for plotting

mat = matread(output_final)
LP_stats = mat["LP_stats"][:,2:end]
Maj_stats = mat["Maj_stats"][:,2:end]
Problem_stats = mat["Problem_stats"][:,2:end]
Total_shared = mat["Total_shared"][2:end]
lp_wasted = mat["lp_wasted"][2:end]
maj_wasted = mat["maj_wasted"][2:end]
ns = Problem_stats[1,:]
lp_nodes = 1 .-lp_wasted./ns
maj_nodes = 1 .-maj_wasted./ns
shared_ratio = Total_shared ./ Problem_stats[1,:]
lp_ratios = LP_stats[2,:]
lp_sats = LP_stats[1,:]
maj_ratios = Maj_stats[2,:]
maj_sats = Maj_stats[1,:]


## Plot the fraction of ingredients that LP-round and Majority Vote agree on
xs = 10:10:200
s1 = 300
s2 = 250
x_label = "\\beta = Ingredient Degree Threshold"
y_label = "Shared Ingredients Fraction"
l_place = :best
ms = 5
lw = 2
title = ""
type1 = :circle
plot(xs,shared_ratio, title = title,
    labels = "LP-round",
    grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw,markerstrokewidth = 0,
    color = :purple, markershape = type1, markersize = ms)

savefig("Plots/SharedRatio.pdf")

## Edge Satisfaction

x_label = "\\beta = Ingredient Degree Threshold"
y_label = "Edge Satisfaction"
l_place = :best
ms = 5
lw = 2
title = ""
type1 = :circle
plot(xs,lp_sats, title = title,
    labels = "LP-round",
    grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw,markerstrokewidth = 0,
    color = :blue, markershape = type1, markersize = ms)

plot!(xs,maj_sats,
    labels = "M.Vote",linewidth = lw,
    markerstrokewidth = 0,color = :red, markershape = type1, markersize = ms)

savefig("Plots/EdgeSat.pdf")

## Wasted ingredients plot
title = ""
x_label = "\\beta = Ingredient Degree Threshold"
y_label = "# Unused Ingredients"
l_place = :bottomright
plot(xs,lp_wasted, title = title,
    labels = "LP-round",
    grid = false, size = (s1,s2),
    xlabel = x_label, ylabel = y_label, legend = l_place,
    linewidth = lw, markerstrokewidth = 0,
    color = :blue, markershape = :circle, markersize = ms)

plot!(xs,maj_wasted,
    labels = "M. Vote",linewidth = lw, markerstrokewidth = 0,
    color = :red, markershape = :circle, markersize = ms)

savefig("Plots/Unused_ingredients.pdf")

## Specifying a value of beta, we can summarize differences between LP-round
# and majority vote in a table of results
println("\n\n")
beta = 20
output_mat = "Output/Beta_$beta"*"_stats.mat"
mat = matread(output_mat)
maj_cuisine = round.(Int64,mat["maj_cuisine"])
lp_cuisine = round.(Int64,mat["lp_cuisine"])
p = sortperm(lp_cuisine[1,:])
LP_Used = lp_cuisine[3,p]
Maj_Used = maj_cuisine[3,p]

sdlp = 0
tlp = 0
sdmaj = 0
tmaj = 0
for cuis = 20:-1:1
    global sdlp, tlp, sdmaj, tmaj
   lp_ing, lp_active, lp_used, lp_rec, lp_waste, lp_diff, lp_over_shared_active, lp_over_shared_used, lp_count = lp_cuisine[:,p[cuis]]
   maj_ing, maj_active, maj_used, maj_rec, maj_waste, maj_diff, maj_over_shared_active, maj_over_shared_used, maj_count = maj_cuisine[:,p[cuis]]
   ratlp = round(lp_diff/lp_count,digits=2)
   ratmaj = round(maj_diff/maj_count,digits = 2)
   sdlp+=lp_diff; tlp+=lp_count; sdmaj+= maj_diff; tmaj+=maj_count
   println(Cuisines_1[p[cuis]]*" & $lp_diff & $lp_count & $ratlp & $maj_diff &  $maj_count & $ratmaj \\\\")
end

println("Total & $sdlp & $tlp& - & $sdmaj &  $tmaj & -\\\\")

## Bar Chart, showing unused ingredients in each cluster
using StatsPlots
beta = 10
numcuisines = 20    # Print the 'numcuisines' smallest cluster sizes


output_mat = "Output/Beta_$beta"*"_stats.mat"
mat = matread(output_mat)
maj_cuisine = round.(Int64,mat["maj_cuisine"])
lp_cuisine = round.(Int64,mat["lp_cuisine"])
p = sortperm(lp_cuisine[1,:])
lps = lp_cuisine[1,p][1:numcuisines]
majs = maj_cuisine[1,p][1:numcuisines]
LP_Used = lp_cuisine[3,p][1:numcuisines]
Maj_Used = maj_cuisine[3,p][1:numcuisines]

bw = 0.5
mn = [lps majs]
sx = repeat(["LP", "MV"], inner = numcuisines)
nam = [Cuisines_1[p[1:numcuisines]]; Cuisines_1[p[1:numcuisines]]]
StatsPlots.groupedbar(nam, mn, group = sx, ylabel = " Ingredients",
bar_width = bw, lw = 0, c = [:blue :red], markerstrokewidth = 1.5,
framestyle = :box, grid = false, yticks = 0:100:1000,xrotation = 60)

mn = [lps-LP_Used majs-Maj_Used]
sx = repeat(["LP Unused","Maj Unused"], inner = 18)
StatsPlots.groupedbar!(nam, mn, group = sx, ylabel = "# Ingredients",legend = :topleft,
bar_width = bw, lw = 0, c = [:cyan :pink], markerstrokewidth = 1.5)

savefig("Plots/BarChart_Unused_Beta_$beta.pdf")
## Showing shared ingredients

beta = 10
last = numcuisines

output_mat = "Output/Beta_$beta"*"_stats.mat"
mat = matread(output_mat)
maj_cuisine = round.(Int64,mat["maj_cuisine"])
lp_cuisine = round.(Int64,mat["lp_cuisine"])
shared = mat["Shared_ing_list"]
shared_num = sum(shared,dims =1)

shared = shared_num[p][1:last]
p = sortperm(lp_cuisine[1,:])
lps = lp_cuisine[1,p][1:last]
majs = maj_cuisine[1,p][1:last]
LP_Used = lp_cuisine[3,p][1:last]
Maj_Used = maj_cuisine[3,p][1:last]
bw = .5
s1 = 500
s2 = 400
mn = [lps majs]
sx = repeat(["LP", "MV"], inner = last)
nam = [Cuisines_1[p[1:last]]; Cuisines_1[p[1:last]]]
StatsPlots.groupedbar(nam, mn, group = sx, ylabel = " Ingredients",size = (s1,s2),
bar_width = bw, lw = 0, c = [:blue :red], markerstrokewidth = 0,
framestyle = :box, grid = false, yticks = 0:100:1000,xrotation = 50)

mn = [shared; shared]
sx = repeat(["Shared"], inner = last)
StatsPlots.groupedbar!(nam, mn, group = sx, ylabel = "# Ingredients",legend = :topleft,
bar_width = bw, lw = 0, c = [:gray],markerstrokewidth = 0)

savefig("Plots/BarChart_Shared_Beta_$beta.pdf")
