#= Robust Recoverable Shortest Path Solver
#
#  Module implementing incremental version of Shortest Path.
#
#  Authors: Marcel Jackiewicz
=#

#= Adds constraints to @a mdl for the path @a y to be in the @a neighbourhood_type neighbourhood, of the path @a x, of the size @a k.
#  
# The added constraints are anonymous, since there might be many variables modeling paths in a model
# and JuMP doesn't support dynamic allocated strings as parameter names.
=#
function addNeighbourhoodConstraints(
    model::JuMP.Model, x::Union{Vector{Bool}, Vector{JuMP.VariableRef}}, y::Vector{JuMP.VariableRef},
    neighbourhood_type::NeighbourhoodType, k::Integer, arcs_num::Integer
)::Nothing
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

    return
end

#= Adds constraints to @a mdl for decision variables in @a x for them to create a simple path in the graph @a g.
#  It is assumed that the constraints for @a x to create a path were added using different function!
#  
# The added constraints are anonymous, since there might be many variables modeling paths in a model
# and JuMP doesn't support dynamic allocated strings as parameter names.
=#
function addAntiCyclicConstraints(mdl::JuMP.Model, x::Vector{JuMP.VariableRef}, g::Graph, s_idx::Integer, t_idx::Integer)::Nothing
    # Arcs entering s or leaving t can not be in a simple s-t path
    JuMP.@constraint(
        mdl,
        [i in 1:length(g.arcs); g.arcs[i].end_node == g.nodes[s_idx] || g.arcs[i].start_node == g.nodes[t_idx]],
        x[i] == 0
    )

    # p represents the place in the ordering of nodes along the path; if the path contains a cycle, the ordering contains a contradiction
    p::Vector{JuMP.VariableRef} = JuMP.@variable(mdl, [i in 1:length(g.nodes)], integer=true, lower_bound=1.0, upper_bound=length(g.nodes))
    big_M::Float64 = length(g.nodes)
    JuMP.@constraint(mdl, [i in 1:length(g.arcs)], p[g.arcs[i].start_node.idx] + big_M*x[i] + 1 <= p[g.arcs[i].end_node.idx] + big_M)

    return
end

"""
    solveIncrementalShortestPath(instance::RrspInstance, x::Path)::RrspSolution

Returns a shortest ``s-t`` path in the graph `instance.graph` with respect to `Cost.second` costs,
where ``s`` and ``t`` are nodes given by `instance.s_idx`, `instance.t_idx` indices in `instance.graph.arcs` array.

If the costs for incremental problem are different than second stage costs lower bounds, they should be set before calling the function.

The path is in the neighbourhood of the path `x`; the neighbourhood type is selected with `instance.neighbourhood`
and the neighbourhood size with `instance.k` values.

The path is stored in `RrspSolution.second_stage_path` of the returned structure and its cost is assigned to `RrspSolution.value`.
For integrality of the solution, the input parameter path `x` is stored in `RrspSolution.first_stage_path`.
However, the first stage cost doesn't affect the `value` of the solution.

The path is computed using compact MIP model.
"""
function solveIncrementalShortestPath(instance::RrspInstance, x::RrspPath)::RrspSolution
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, Rrsp.optimizer)

    arc_num::Integer = length(instance.graph.arcs)

    # y[i] --- is the i-th arc taken?
    JuMP.@variable(model, y[i in 1:arc_num], Bin)

    addPathConstraints(model, model[:y], instance.graph, instance.s_idx, instance.t_idx)
    addNeighbourhoodConstraints(model, x.arcs, model[:y], instance.neighbourhood, instance.k, arc_num)

    JuMP.@objective(
        model, Min, sum(y[i]*instance.graph.arcs[i].cost.second for i in 1:arc_num)  # assume the adversary set the cost
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return createEmptyRrspSolution(arc_num)
    end

    return RrspSolution(x, RrspPath([JuMP.value(y[i]) > 0.5 for i in 1:arc_num]), JuMP.objective_value(model))
end
