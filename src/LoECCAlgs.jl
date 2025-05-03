using SparseArrays
using LinearAlgebra
using StatsBase
using JuMP
using Gurobi
using Clustering
using Base.Sort
using DataStructures

gurobi_env = Gurobi.Env()

"""
LoECCCanonicalLP
Solves the Locally Budgeted Overlapping Edge-Colored Clustering LP relaxation

Input:
    EdgeList = the edge list of a LO-ECC instance
    EdgeColors = the edge colors
    n = number of vertices in the instance
    b = local (per vertex) budget for color assignments
    optimalflag = True iff you want a solution to the integer program. Defaults to false.
    outputflag = 1 for verbose output from Gurobi. Defaults to 0.

Output:
    LPval = The LP relaxation objective score. This is a lower bound on OPT (equal if you set optimalflag=true)
    X = The distance labels from solving the objective
    runtime = Solver runtime
"""
function LoECCCanonicalLP(EdgeList::Vector{Vector{Int64}}, EdgeColors::Array{Int64, 1}, n::Int64, b::Int64, optimalflag::Bool=false,outputflag::Int64=0)

    k = maximum(EdgeColors)
    M = length(EdgeList)

    # m = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=outputflag, gurobi_env))
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => outputflag))

    # Variables for nodes and edges
    @variable(m, y[1:M])

    @objective(m, Min, sum(y[i] for i=1:M))

    if optimalflag
        @variable(m, x[1:n,1:k],Bin)
    else
        @variable(m, x[1:n,1:k])
        @constraint(m,x .<= ones(n,k))
        @constraint(m,x .>= zeros(n,k))
        @constraint(m,y .<= ones(M))
        @constraint(m,y .>= zeros(M))
    end

    for i = 1:n
        @constraint(m, sum(x[i,j] for j = 1:k) >= k-b)
    end

    for e = 1:M
        color = EdgeColors[e]
        edge = EdgeList[e]

        # For every node in the edge, there's a constraint for the
        # node-color variable
        for v = edge
            @constraint(m, y[e] >= x[v,color])
        end

    end
    start = time()
    JuMP.optimize!(m)
    runtime = time()-start

    # Return clustering and objective value
    X = JuMP.value.(x)
    LPval= JuMP.objective_value(m)

    return LPval, X, runtime
end

"""
LoECCBicriteriaRound

A rounding scheme for the LP defined in LoECCBicriteriaLP. Produces a
(1/(1-epsilon), epsilon) approximation, where the first factor is on edge-penalties
and the second factor is on local clustering budget.

Input:
    EdgeList = the edges
    EdgeColors = the edge colors
    X = the distance labels from solving the LP relaxation
    LPval = the objective score from solving the LP relaxation
    b = the local budget for cluster assignments

Output:
    c = The cluster assignments
    RoundScore = The LO-ECC objective score of the rounded clustering c.
    RoundRatio = The approximation ratio of RoundScore with respect to LPval.
    RoundBudgetScore = The maximum number of colors assigned to any vertex by c.
    RoundBudgetRatio = The approximation score of RoundBudgetScore with respect to b.

"""
function LoECCBicriteriaRound(EdgeList::Union{Array{Int64,2}, Vector{Vector{Int64}}}, EdgeColors::Array{Int64, 1}, X::Array{Float64,2},LPval::Float64, b::Int64, epsilon::Float64)
    c = Vector{Vector{Int64}}()
    n = size(X, 1)
    maxcolors = 0
    for i = 1:n
        colors = Vector{Int64}()
        for j = 1:maximum(EdgeColors)
            if X[i, j] < epsilon
                push!(colors, j)
            end
        end
        push!(c, colors)
        maxcolors = maximum([maxcolors, length(colors)])
    end
    RoundScore = OverlappingEdgeCatClusObj(EdgeList, EdgeColors, c)
    RoundRatio = 1
    if LPval != 0
        RoundRatio = RoundScore/LPval
    end
    RoundBudgetScore = maxcolors
    RoundBudgetRatio = maxcolors/b
    return c, RoundScore, RoundRatio, RoundBudgetScore, RoundBudgetRatio
end

"""
OverlappingEdgeCatClusObj

Returns the number of mistakes made by an overlapping clustering in an instance of
a Categorical Edge Clustering problem.
"""
function OverlappingEdgeCatClusObj(EdgeList::Union{Array{Int64,2}, Vector{Vector{Int64}}}, EdgeColors::Array{Int64,1}, c::Vector{Vector{Int64}})
    n = length(c)
    mistakes = 0
    for i = 1:size(EdgeList,1)
        if size(EdgeList,2) == 2
            edge = EdgeList[i,:]
        else
            edge = EdgeList[i]
        end
        for v in edge
            if !(EdgeColors[i] in c[v])
                mistakes += 1
                break
            end
        end
    end
    return mistakes
end


"""
GreedyLocal
Returns a clustering according to the LO-ECC greedy algorithm. Each node is
assigned its b most frequent colors.

Input:
    EdgeList = the edges
    EdgeColors = the edge colors
    n = number of nodes
    k = number of colors
    b = the local budget for cluster assignments

Outputs:
    greedy_c = the cluster assignments
"""

function GreedyLocal(EdgeList::Vector{Vector{Int64}},EdgeColors::Vector{Int64},n::Int64,k::Int64,b::Int64)

    b = min(b, k)   # simplification

    # maintain a priority queue of color degrees for each node
    ColorDegree = Vector{PriorityQueue{Int64, Int64}}()

    # initialize priority queues
    for i = 1:n
        push!(ColorDegree, PriorityQueue{Int64, Int64}(Base.Order.Reverse))
        for j = 1:k
            ColorDegree[i][j] = 0
        end
    end

    # populate color degrees
    for t = 1:length(EdgeList)
        edge = EdgeList[t]
        color = EdgeColors[t]
        for node in edge
            ColorDegree[node][color] += 1
        end
    end

    # form the clustering
    greedy_c = Vector{Vector{Int64}}()
    for i = 1:n
        push!(greedy_c, Vector{Int64}())
        for j = 1:b
            push!(greedy_c[i], dequeue!(ColorDegree[i]))
        end
    end
    return greedy_c
end