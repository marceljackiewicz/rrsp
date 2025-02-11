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

    JuMP.@objective(model, Min, sum(y[i]*instance.graph.arcs[i].cost.second for i in 1:length(instance.graph.arcs)))  # assume the adversary set the cost

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return Path([])
    end

    path::Path = Path([JuMP.value(y[i]) > 0.5 for i in 1:length(instance.graph.arcs)])
    return path
end
