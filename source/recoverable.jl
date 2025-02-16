#= Robust Recoverable Shortest Path Solver
#
#  Module implementing recoverable version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

#= Returns first and second stage s-t paths and the optimal objective function value for Recoverable Shortest Path,
#  to which the RRSP problem with interval uncertainty reduces to.
#
#  The paths are computed using compact ILP model.
=#
function solveRecoverableShortestPath(instance::RrspInstance)::RrspSolution
    optimizer = Cbc.Optimizer
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    # x[i], y[i] --- is the i-th arc in the first and second stage path, respectively, taken?
    JuMP.@variable(model, x[i in 1:length(instance.graph.arcs)], Bin)
    JuMP.@variable(model, y[i in 1:length(instance.graph.arcs)], Bin)

    addPathConstraints(model, model[:x], instance.graph, instance.s_idx, instance.t_idx)
    addPathConstraints(model, model[:y], instance.graph, instance.s_idx, instance.t_idx)

    #= if the graph contains cycles, the anti-cyclic constraints must be used for the first stage path,
       as taking more arcs may prove beneficial in the second stage =#
    addAntiCyclicConstraints(model, model[:x], instance.graph, instance.s_idx, instance.t_idx)
    # the anti-cyclic constraints can be droped for the second stage path iff inclusion neighbourhood is used
    if instance.neighbourhood != INCLUSION
        addAntiCyclicConstraints(model, model[:y], instance.graph, instance.s_idx, instance.t_idx)
    end

    addNeighbourhoodConstraints(model, model[:x], model[:y], instance.neighbourhood, instance.k, length(instance.graph.arcs))

    JuMP.@objective(
        model,
        Min,
        sum(
            x[i]*instance.graph.arcs[i].cost.first + y[i]*(instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta)
            for i in 1:length(instance.graph.arcs)
        )
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return RrspSolution(Path([]), Path([]), Inf)
    end

    x_path::Path = Path([JuMP.value(x[i]) > 0.5 for i in 1:length(instance.graph.arcs)])
    y_path::Path = Path([JuMP.value(y[i]) > 0.5 for i in 1:length(instance.graph.arcs)])
    return RrspSolution(x_path, y_path, JuMP.objective_value(model))
end
