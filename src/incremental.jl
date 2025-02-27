#= Robust Recoverable Shortest Path Solver
#
#  Module implementing incremental version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
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

#= Returns a shortest s-t path in the graph @a g, where s and t are nodes given by @a s_idx, t_idx indices in @a g arcs array.
#
#  The path is computed using LP model.
=#
function solveIncrementalShortestPath(instance::RrspInstance, x::Path)::Path
    optimizer = Cbc.Optimizer
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    # y[i] --- is the i-th arc taken?
    JuMP.@variable(model, y[i in 1:length(instance.graph.arcs)] >= 0)

    addPathConstraints(model, model[:y], instance.graph, instance.s_idx, instance.t_idx)
    addNeighbourhoodConstraints(model, x.arcs, model[:y], instance.neighbourhood, instance.k, length(instance.graph.arcs))

    JuMP.@objective(
        model, Min, sum(y[i]*instance.graph.arcs[i].cost.second for i in 1:length(instance.graph.arcs))  # assume the adversary set the cost
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return Path([])
    end

    path::Path = Path([JuMP.value(y[i]) > 0.5 for i in 1:length(instance.graph.arcs)])
    return path
end
