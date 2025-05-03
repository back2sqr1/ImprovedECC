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


# randomly take a subset of the incidence matrix for testing

l, u = size(incidence_form)
l = 100
u = 100
bipartite_graph = zeros(Int64, l + u, l + u)
for i = 1:l
    for j = 1:u
        if incidence_form[i, j] == 1
            bipartite_graph[i, l + j] = 1
            bipartite_graph[l + j, i] = 1
        end
    end
end

# Plot the bipartite graph
using GraphPlot
using Graphs
using Colors


g = SimpleGraph(l + u)
for i = 1:l
    for j = 1:u
        if incidence_form[i, j] == 1
            add_edge!(g, i, l + j)
        end
    end
end

# Plot the bipartite graph delete nodes with no edges 

n = nv(g)                        # number of vertices
for v in n:-1:1                  # v = n, n-1, … , 1
    if outdegree(g, v) == 0
        rem_vertex!(g, v)        # remove vertex v
    end
end

function find_bipartite_sets(g::AbstractGraph)
    colors = fill(0, nv(g))       # 0 = unvisited, 1 / 2 = partitions
    left, right = Int[], Int[]

    for v in vertices(g)                     # handle disconnected graphs
        colors[v] == 0 || continue
        queue = [v]; colors[v] = 1; push!(left, v)

        while !isempty(queue)
            u = popfirst!(queue)
            for w in neighbors(g, u)
                if colors[w] == 0
                    colors[w] = 3 - colors[u]
                    push!(colors[w] == 1 ? left : right, w)
                    push!(queue, w)
                elseif colors[w] == colors[u]
                    error("Graph is not bipartite.")
                end
            end
        end
    end
    return (left, right)
end

left, right = find_bipartite_sets(g)
function bipartite_xy(left, right; xgap = 1.0, ygap = 3.0)
    xs = vcat(fill(0.0, length(left)), fill(xgap, length(right)))
    ys = vcat(range(0, step = -ygap, length = length(left)),
              range(0, step = -ygap, length = length(right)))
    return xs, ys
end
colours = fill(colorant"blue", nv(g))    # default = blue
colours[right] .= colorant"red"          # overwrite right–side vertices

xs, ys = bipartite_xy(left, right; xgap = 1)   # tweak gaps as you like
gplot(g, xs, ys, nodefillc=colours)
gplot(g, nodefillc=colours)

