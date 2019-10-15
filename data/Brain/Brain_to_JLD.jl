using MAT
mat1 = matread("original_data/Coactivation_matrix.mat")
mat2 = matread("original_data/GroupAverage_rsfMRI_matrix.mat")

A1 = mat1["Coactivation_matrix"]
A2 = mat2["GroupAverage_rsfMRI"]
n = size(A1,1)

EdgeColors = Vector{Int64}()
EdgeList = Vector{Vector{Int64}}()
for i = 1:n
    for j = i+1:n
        if A1[i,j] > 0 && A2[i,j] == 0.0
            push!(EdgeList,[i;j])
            push!(EdgeColors,1)
        elseif A2[i,j] > 0 && A1[i,j] == 0.0
            push!(EdgeList,[i;j])
            push!(EdgeColors,2)
        end
    end
end
save("../JLD_Files/Brain.jld","EdgeList", EdgeList, "EdgeColors", EdgeColors, "n", n)
