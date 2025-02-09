#= Robust Recoverable Shortest Path Solver
#
#  Module implementing incremental version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

# TODO
# - add model for incremental
#   [x] move "path constraints" to a function (see deterministic.jl)
#   [x] add constraints for different neighbourhoods
#   [-] add anti-cyclic constraints for model for general digraphs

#= Adds constraints to @a mdl for the path @a y to be in the @a neighbourhood_type neighbourhood, of the path @a x, of the size @a k.
#  
# The added constraints are anonymous, since there might be many variables modeling paths in a model
# and JuMP doesn't support dynamic allocated strings as parameter names.
=#
function addNeighbourhoodConstraints(
    model::JuMP.Model, x::Vector{JuMP.VariableRef}, y::Vector{JuMP.VariableRef},
    neighbourhood_type::NeighbourhoodType, k::Integer, arcs_num::Integer
)
    # z, to model neighbourhood containment, represents the xy product
    z::Vector{JuMP.VariableRef} = JuMP.@variable(model, [i in 1:arcs_num], lower_bound = 0.0)

    # constraints for second stage different arcs
    if neighbourhood_type == INCLUSION
        JuMP.@constraint(model, sum(y[i] - z[i] for i in 1:arcs_num) <= k)
    elseif neighbourhood_type == EXCLUSION
        JuMP.@constraint(model, sum(x[i] - z[i] for i in 1:arcs_num) <= k)
    elseif neighbourhood_type == SYM_DIFF
        JuMP.@constraint(model, sum(y[i] + x[i] - 2*z[i] for i in 1:arcs_num) <= k)
    end

    JuMP.@constraint(model, [i in 1:arcs_num], z[i] <= x[i])  # z doesn't have an arc if x doesn't have it
    JuMP.@constraint(model, [i in 1:arcs_num], z[i] <= y[i])  # z doesn't have an arc if y doesn't have it
end

#= Adds constraints to @a mdl for decision variables in @a x to create an s-t path in the graph @a g.
#  
# The added constraints are anonymous, since there might be many variables modeling paths in a model
# and JuMP doesn't support dynamic allocated strings as parameter names.
=#
function addPathConstraints(model::JuMP.Model, x::Vector{JuMP.VariableRef}, g::Graph, s_idx::Integer, t_idx::Integer)
    function sumOutFlowAtNodeIdx(node_idx::Integer)::JuMP.AffExpr
        out_arcs_indices = [i for i in 1:length(g.arcs) if g.arcs[i].start_node == g.nodes[node_idx]]
        return length(out_arcs_indices) > 0 ? sum(x[i] for i in out_arcs_indices) : 0
    end

    function sumInFlowAtNodeIdx(node_idx::Integer)::JuMP.AffExpr
        in_arcs_indices = [i for i in 1:length(g.arcs) if g.arcs[i].end_node == g.nodes[node_idx]]
        return length(in_arcs_indices) > 0 ? sum(x[i] for i in in_arcs_indices) : 0
    end

    function allInternalNodesIndices()::Vector{Integer}
        return [i for i in 1:length(g.nodes) if i != s_idx && i != t_idx]
   end

    # At any node (besides @a s and @a t) the flow is conserved
    JuMP.@constraint(
        model, [i in allInternalNodesIndices()], sumInFlowAtNodeIdx(i) == sumOutFlowAtNodeIdx(i)
    )  # flow conservation

    # Flow starts (flow in + 1 = flow out) and ends (flow in - 1 = flow out) in the specified nodes
    JuMP.@constraint(model, sumOutFlowAtNodeIdx(s_idx) - sumInFlowAtNodeIdx(s_idx) ==  1)  # flow_source
    JuMP.@constraint(model, sumOutFlowAtNodeIdx(t_idx) - sumInFlowAtNodeIdx(t_idx) == -1)  # flow_sink
end
