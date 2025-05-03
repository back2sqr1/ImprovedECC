using MAT
using LinearAlgebra
using MAT
using SparseArrays
using Plots
include("helpers.jl")


function test_incidence2elist()
    # Test incidence2elist
    H = sparse([1 1 0 0; 0 1 1 0; 1 0 1 0; 0 0 1 1])
    H = convert(SparseArrays.SparseMatrixCSC{Float64,Int64}, H)
    Hyperedges = incidence2elist(H)
    return Hyperedges == [[1, 2], [2, 3], [1, 3], [3, 4]]
end

function test_elist2incidence()
    # Test elist2incidence
    Hyperedges = [[1, 2], [2, 3], [1, 3], [3, 4]]
    H = elist2incidence(Hyperedges, size(Hyperedges)[1])
    x = sparse([1 1 0 0; 0 1 1 0; 1 0 1 0; 0 0 1 1])
    x = convert(SparseArrays.SparseMatrixCSC{Float64,Int64}, x)
    return H == x
end

if test_incidence2elist() && test_elist2incidence()
    println("All tests passed!")
else
    println("Tests failed!")
end


