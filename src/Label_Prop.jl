using MAT
using LinearAlgebra
using MAT
using SparseArrays
using Plots
using StatsBase
using JSON
using JLD
include("helpers.jl")

cooking = JSON.parsefile("include\\CategoricalEdgeClustering-master\\data\\Cooking\\train.json")


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

incidence_form = elist2incidence(EdgeList, length(Ingredients))

"""
Incidence matrix form:  H[e,u] = 1  iff node u is in hyperedge e

l x u matrix where l is the number of hyperedges and u is the number of nodes
"""

