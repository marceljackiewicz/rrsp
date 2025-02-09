#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

import DataStructures

import Cbc
import JuMP

struct Node
    uid::Integer
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

@enum NeighbourhoodType INCLUSION EXCLUSION SYM_DIFF

struct Path
    arcs::Vector{Bool}  # path characteristic vector
end

struct Graph
    nodes::Vector{Node}
    arcs::Vector{Arc}
end

@enum AspTreeOp SERIES PARALLEL NONE

struct AspTreeNode
    s::Node
    t::Node
    left::Integer
    right::Integer
    arcs::Array{Arc}
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
    k::Integer
    gamma::Float64
end

struct RrspSolution
    first_stage_path::Path
    second_stage_path::Path
    value::Float64
end
