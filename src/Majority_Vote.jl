using MAT
using LinearAlgebra
using MAT
using SparseArrays
using Plots
using StatsBase
include("helpers.jl")

function majority_vote(node_list_form :: Vector{Vector{Int64}}, colors :: Vector{Int64})
    # I want to color the nodes such that hyperedges tend to contain nodes that all have the same color as the hyperedge
    # The color of my node is the mode of the colors of the hyperedges it is in
    n = length(node_list_form)
    color_majority = zeros(Int64, n)
    for v = 1:n
        list_of_edges = node_list_form[v]
        color_array = colors[list_of_edges] # colors of the hyperedges that contain this node
        color_majority[v] = StatsBase.mode(color_array)
    end
end

M = matread("trivago-dataset/Trivago_Clickout_EdgeLabels.mat")

# Printing Binary Incidence Matrix
incidence_form = M["H"]
edge_list_form = incidence2elist(incidence_form)
node_list_form = incidence2elist(incidence_form, true)
incidence_form = elist2incidence(edge_list_form, size(H)[1])

# what are node labels vs edge labels? Labels labelNames?


# write out connection between the two from EC tto SSL
# one hot encoding for multiple colors

# pick maximum for possibility of color
# construct bipartitie graph


# left nodes, right sides positive or negative whether red or not (color or not)
# Zhou Bousquet 2004
