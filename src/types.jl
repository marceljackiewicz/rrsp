#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz
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

struct Graph
    nodes::Vector{Node}
    arcs::Vector{Arc}
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

"""
    struct RrspPath
        arcs::Vector{Bool}
    end

The structure for encapsulating a path of the problem instance graph.
Paths are represented by characteristic vectors of the arc set of the instance graph
as this representation allows for easier constraints formulation in MIP modeling -- the principal aspect of the solver.
"""
struct RrspPath
    arcs::Vector{Bool}  # path characteristic vector
end

function createEmptyPath(arcs_cardinality::Integer)::RrspPath
    return RrspPath([0 for _ in 1:arcs_cardinality])
end

"""
    struct RrspInstance
        graph::Graph
        s_idx::Integer
        t_idx::Integer
        neighbourhood::NeighbourhoodType
        k::Integer
        gamma::Float64
    end

The structure containing serialized RRSP problem instance.
It can be used as an opaque type, since the solving functions only input parameter is the instance itself.
"""
struct RrspInstance
    graph::Graph
    s_idx::Integer
    t_idx::Integer
    neighbourhood::NeighbourhoodType
    k::Integer
    gamma::Float64
end

"""
    struct RrspSolution
        first_stage_path::RrspPath
        second_stage_path::RrspPath
        value::Float64
    end

The structure containing optimal pair of paths and the optimal value of the solution to the RRSP problems.

Every solver API function returns `RrspSolution` structure object.
When solution is infeasible, the `RrspSolution` structure is returned *empty* -- the paths structures
are initialized to ``\\mathbf{0}`` vectors and `value` is `Inf`.

For convenience, the same structure is used when returning the solution for RRSP subproblems,
which concern only one optimal path. In this case, as to which path of the set is given by function docstring.
"""
mutable struct RrspSolution
    first_stage_path::RrspPath
    second_stage_path::RrspPath
    value::Float64
end

function createEmptyRrspSolution(arcs_cardinality::Integer)::RrspSolution
    return RrspSolution(createEmptyPath(arcs_cardinality), createEmptyPath(arcs_cardinality), Inf)
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

mutable struct AspNodeData
    opt_first_stage_path::RrspPath
    opt_second_stage_path::Array{RrspPath}
    opt_solution_paths::Array{RrspSolution}
end
