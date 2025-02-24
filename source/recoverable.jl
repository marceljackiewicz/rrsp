#= Robust Recoverable Shortest Path Solver
#
#  Module implementing recoverable version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

#= Returns first and second stage s-t paths and the optimal objective function value for Recoverable Shortest Path,
#  to which the RRSP problem with interval uncertainty reduces to.
#
#  The paths are computed using compact MIP model.
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


#= Returns approximation ratio of approximating RRSP in acyclic digraphs and discrete budget with RecSP.
#  Works iff second stage cost for each arc is strictly > 0!
=#
function getRecSpToRrspAcyclicDiscreteBudgetApproxRatio(instance::RrspInstance)::Float64
    ratio::Float64 = maximum(
        (instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta)/instance.graph.arcs[i].cost.second
        for i in 1:length(instance.graph.arcs)
    )

    return ratio
end

#= Returns approximation ratio of approximating RRSP in acyclic digraphs and continuous budget with RecSP.
#  Works iff second stage cost for each arc is strictly > 0!
=#
function getRecSpToRrspAcyclicContBudgetApproxRatio(instance::RrspInstance)::Float64
    ratio_alpha::Float64 = maximum(
        (instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta)/instance.graph.arcs[i].cost.second
        for i in 1:length(instance.graph.arcs)
    )

    total_deviation_possible::Float64 = sum(instance.graph.arcs[i].cost.delta for i in 1:length(instance.graph.arcs))
    ratio_beta::Float64 = total_deviation_possible >= instance.gamma ? total_deviation_possible/instance.gamma : Inf

    # TODO: add computing X^ from incremental and F(X^) from adversarial
    ratio_gamma::Float64 = Inf

    return min(ratio_alpha, ratio_beta, ratio_gamma)
end
