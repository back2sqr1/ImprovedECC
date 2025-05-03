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
# get some of the dataset
for entry in cooking
    for i in entry["ingredients"]
        push!(Ingredients,i)
    end
    push!(Cuisines,entry["cuisine"])
end

c = length(Cuisines)

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

Num2Ing = Dict()
Num2Cu = Dict()
for i = 1:length(Ingredients)
    Num2Ing[i] = Ingredients[i]
end
for i = 1:length(Cuisines)
    Num2Cu[i] = Cuisines[i]
end
Num2Cu
Num2Ing


incidence_form = elist2incidence(EdgeList, length(Ingredients))

"""
Incidence matrix form:  H[e,u] = 1  iff node u is in hyperedge e

l x u matrix where l is the number of hyperedges and u is the number of nodes
"""

l, u = size(incidence_form)

bipartite_graph = zeros(Int64, l + u, l + u)
for i = 1:l
    for j = 1:u
        if incidence_form[i, j] == 1
            bipartite_graph[i, l + j] = 1
            bipartite_graph[l + j, i] = 1
        end
    end
end

alp = 0.99

EdgeColors
y = zeros(Int64, l + u, c)
for i = 1:l
    y[i, EdgeColors[i]] = 1
end 
B = sparse(bipartite_graph)                 # (l+u) × (l+u) adjacency

# row-normalise B  →  transition matrix P 
rowsums = vec(sum(B, dims = 2))
invrows  = 1.0 ./ rowsums
invrows[isinf.(invrows)] .= 0.0             # protect isolated vertices
Dinv     = spdiagm(0 => invrows)            # D⁻¹
P        = Dinv * B                         # each row now sums to 1

# linear system
alp   = 0.25
RHS = (1 - alp) * y
LHS = I - alp * P                             # NOT (I-α)*P

# solve linear system
x = LHS \ RHS
x
# find nodes with similar labels and get their names from CuNum and IngNum


highest = zeros(Int64, l + u)
for i = 1:(l + u)
    h = 0
    for j = 1:(c)
        if x[i, j] > h
            h = x[i, j]
            highest[i] = j
        end
    end

    j = highest[i]
    println("Node $i has label $j")
    if i <= l
        println("Node $i is a recipe with cuisine $(Num2Cu[j])")
    else
        println("Node $i is a ingredient $(Num2Ing[i - l]) is an ingredient with cuisine $(Num2Cu[j])")
    end
        
    
end

# find nodes with similar labels and get their names from Num2Cu and Num2Ing
sig = 1
for i = 1:100
    for j = 1:100
        similarity = exp(norm(x[i, :] - x[j, :]))^2 / (2 * sig)
        if similarity > 4.5
            println("Node $i and Node $j are similar with similarity $similarity")
            if i <= l
                println("Node $i is a recipe with cuisine $(Num2Cu[highest[i]])")
                for k = 1:length(EdgeList[i])
                    println("Node $i is connected to ingredient $(Num2Ing[EdgeList[i][k]])")
                end
            else
                println("Node $i is an ingredient $(Num2Ing[i - l]) is an ingredient with cuisine $(Num2Cu[highest[i]])")
            end
            if j <= l
                println("Node $j is a recipe with cuisine $(Num2Cu[highest[j]])")
                for k = 1:length(EdgeList[j])
                    println("Node $j is connected to ingredient $(Num2Ing[EdgeList[j][k]])")
                end
            else
                println("Node $j is an ingredient $(Num2Ing[j - l]) is an ingredient with cuisine $(Num2Cu[highest[j]])")
            end
        end
    end
end