#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

struct Node
    uid::Integer
    idx::Integer
end

struct Cost
    first::Float64
    second::Float64
    delta::Float64
end

struct Arc
    start_node::Node
    end_node::Node
    cost::Cost
end

@enum NeighbourhoodType INCLUSION EXCLUSION SYM_DIFF NEIGHBOURHOOD_NOT_SET

function stringToNeighbourhoodType(str::AbstractString)::NeighbourhoodType
    if str == "INC"
        return INCLUSION
    elseif str == "EXC"
        return EXCLUSION
    elseif str == "SYM_DIFF"
        return SYM_DIFF
    else
        return NEIGHBOURHOOD_NOT_SET
    end
end

struct Path
    arcs::Vector{Bool}  # path characteristic vector
end

function createEmptyPath(arcs_cardinality::Integer)::Path
    return Path([0 for _ in 1:arcs_cardinality])
end

struct Graph
    nodes::Vector{Node}
    arcs::Vector{Arc}
end

@enum AspTreeOp SERIES PARALLEL NONE

mutable struct AspTreeNode
    s::Node
    t::Node
    left::Integer
    right::Integer
    arc_idx::Integer  # the index of an arc if the node is a leaf
    operation::AspTreeOp
    is_leaf::Bool
end

struct AspTree
    root_idx::Integer
    nodes::Vector{AspTreeNode}
end

struct AspComponent
    s::Node
    t::Node
end

struct RrspInstance
    graph::Graph
    s_idx::Integer
    t_idx::Integer
    neighbourhood::NeighbourhoodType
    k::Integer
    gamma::Float64
end

mutable struct RrspSolution
    first_stage_path::Path
    second_stage_path::Path
    value::Float64
end

function createEmptyRrspSolution(arcs_cardinality::Integer)::RrspSolution
    return RrspSolution(createEmptyPath(arcs_cardinality), createEmptyPath(arcs_cardinality), Inf)
end

mutable struct AspNodeData
    opt_first_stage_path::Path
    opt_second_stage_path::Array{Path}
    opt_solution_paths::Array{RrspSolution}
end
