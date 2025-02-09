#= Robust Recoverable Shortest Path Solver
#
# Module implementing deterministic version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

import JuMP

include("types.jl")

#= Adds constraints to @a mdl for decision variables in @a x to create an s-t path in the graph @a g.
#  
# The added constraints are anonymous, since there might be many variables modeling paths in a model
# and JuMP doesn't support dynamic allocated strings as parameter names.
=#
function addPathConstraints(mdl::JuMP.Model, x::Vector{JuMP.VariableRef}, g::Graph, s_idx::Integer, t_idx::Integer)
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
        mdl, [i in allInternalNodesIndices()], sumInFlowAtNodeIdx(i) == sumOutFlowAtNodeIdx(i)
    )  # flow conservation

    # Flow starts (flow in + 1 = flow out) and ends (flow in - 1 = flow out) in the specified nodes
    JuMP.@constraint(mdl, sumOutFlowAtNodeIdx(s_idx) - sumInFlowAtNodeIdx(s_idx) ==  1)  # flow_source
    JuMP.@constraint(mdl, sumOutFlowAtNodeIdx(t_idx) - sumInFlowAtNodeIdx(t_idx) == -1)  # flow_sink
end
