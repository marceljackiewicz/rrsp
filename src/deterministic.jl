""" Robust Recoverable Shortest Path Solver

Module implementing deterministic version of Shortest Path.

Authors: Marcel Jackiewicz
"""

#= Adds constraints to @a mdl for decision variables in @a x to create an s-t path in the graph @a g.
#  
# The added constraints are anonymous, since there might be many variables modeling paths in a model
# and JuMP doesn't support dynamic allocated strings as parameter names.
=#
function addPathConstraints(mdl::JuMP.Model, x::Vector{JuMP.VariableRef}, g::Graph, s_idx::Integer, t_idx::Integer)::Nothing
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

    return
end

"""
    solveDeterministicShortestPath(instance::RrspInstance)::RrspSolution

Returns a shortest ``s-t`` path in the graph `instance.graph` with respect to `Cost.first` costs,
where ``s`` and ``t`` are nodes given by `instance.s_idx`, `instance.t_idx` indices in `instance.graph.arcs` array.

The path is stored in `RrspSolution.first_stage_path` of the returned structure and its cost is assigned to `RrspSolution.value`.

The path is computed using LP model.
"""
function solveDeterministicShortestPath(instance::RrspInstance)::RrspSolution
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, Rrsp.optimizer)

    arc_num::Integer = length(instance.graph.arcs)

    # x[i] --- is the i-th arc taken?
    JuMP.@variable(model, x[i in 1:arc_num] >= 0)
    addPathConstraints(model, model[:x], instance.graph, instance.s_idx, instance.t_idx)

    JuMP.@objective(model, Min, sum(x[i]*instance.graph.arcs[i].cost.first for i in 1:arc_num))

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return createEmptyRrspSolution(arc_num)
    end

    return RrspSolution(RrspPath([JuMP.value(x[i]) > 0.5 for i in 1:arc_num]), RrspPath([]), JuMP.objective_value(model))
end
