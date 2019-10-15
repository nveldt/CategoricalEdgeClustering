# print results for running all algorithms on the real world datasets
using MAT

println("")
datasets = ["Brain","MAG-10","Cooking","DAWN","Walmart-Trips"]

mat = matread("Output/dataset_stats.mat")
dataset_stats = round.(Int64,mat["dataset_stats"])

for j = 1:length(datasets)
    dataset = datasets[j]
    mat = matread("Output/"*dataset*"_CCC_results.mat")
    cbapp = round(mat["cb_ratio"], digits = 2)
    lcbapp = round(mat["lcb_ratio"], digits = 2)
    cbsat = round(mat["cb_edgesat"], digits = 2)
    lcbsat = round(mat["lcb_edgesat"], digits = 2)
    lcbtime = round(mat["lcb_run"],digits = 1)
    cbtime = round(mat["cb_run"],digits = 1)


    matiso = matread("Output/"*dataset*"_results_isocut.mat")
    isosat = round(matiso["iso_edgesat"], digits = 2)
    isoapp = round(matiso["iso_ratio"], digits = 2)
    isotime = round(matiso["iso_run"],digits = 1)
    majsat = round(matiso["maj_edgesat"], digits = 2)
    majapp = round(matiso["maj_ratio"], digits = 2)
    majtime = round(matiso["maj_run"],digits = 1)

    matlp = matread("Output/"*dataset*"_results.mat")
    lpsat = round(matlp["lp_edgesat"], digits = 2)
    lpapp = round(matlp["lp_ratio"], digits = 2)
    lptime = round(matlp["LPtime"],digits = 1)

    n = dataset_stats[j,1]
    M = dataset_stats[j,2]
    emax = dataset_stats[j,3]
    k = dataset_stats[j,4]
    println(dataset*" & $n & $M & $emax & $k &$lpapp & $majapp & $isoapp & $cbapp & $lcbapp & $lpsat & $majsat & $isosat & $cbsat & $lcbsat & $lptime & $majtime & $isotime & $cbtime & $lcbtime\\\\")
end
