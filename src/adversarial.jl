#= Robust Recoverable Shortest Path Solver
#
#  Module implementing adversarial version of Shortest Path.
#
#  Authors: Marcel Jackiewicz
=#

"""
    solveAdversarialShortestPathDAGContBudget(instance::RrspInstance, x::Path)::RrspSolution

    Returns optimal second stage shortest ``s-t`` path in the graph `instance.graph` in the neighbourhood of the path `x`
    for Adversarial Shortest Path in DAGs with continuous budgeted uncertainty.
    The nodes ``s`` and ``t`` are given by `instance.s_idx`, `instance.t_idx` indices in `instance.graph.arcs` array.
    The maximum budget value is equal to `instance.gamma`.

    The path is computed using compact MIP model.

    Note: Use only when the instance graph is acyclic!
"""
function solveAdversarialShortestPathDAGContBudget(instance::RrspInstance, x::RrspPath, eval_paths_num::Integer=0)::RrspSolution
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, Rrsp.optimizer)

    second_stage_paths_num::Integer = eval_paths_num == 0 ? length(instance.graph.arcs) + 1 : eval_paths_num

    # second stage paths characteristic vectors, y[i, j] --- is the j-th arc taken in the i-th path?
    JuMP.@variable(model, y[i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)], Bin)

    for i in 1:second_stage_paths_num
        # constraints for arcs in y[i] to create a path
        addPathConstraints(model, model[:y][i, :], instance.graph, instance.s_idx, instance.t_idx)
        # neighbourhood constraints for y[i]
        addNeighbourhoodConstraints(model, x.arcs, model[:y][i, :], instance.neighbourhood, instance.k, length(instance.graph.arcs))
    end

    #= the sum of flows on each y[i] (evaluation) path is 1.0 =#
    JuMP.@variable(model, lambda[i in 1:second_stage_paths_num] >= 0)
    JuMP.@constraint(model, sum(lambda) == 1.0)

    #= v[i,j] linearizes (to represent) the product of y[i,j]*lambda[i] =#
    #TODO(mj): move linearizing product of two variables into a separate function
    JuMP.@variable(model, v[i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)] >= 0)
    big_M::Float64 = 1.0  # y is a binary variable and 0.0 <= lambda <= 1.0
    #= y = 0 -> v = 0 =#
    JuMP.@constraint(
        model,
        [i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)],
        v[i,j] <= big_M*y[i,j]
    )
    #= y = 1 -> v = lambda =#
    JuMP.@constraint(
        model,
        [i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)],
        v[i,j] <= lambda[i]
    )
    JuMP.@constraint(
        model,
        [i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)],
        v[i,j] >= lambda[i] - big_M*(1 - y[i,j])
    )

    #= some budget must be spent to send some flow with y[i] path over an arc =#
    JuMP.@variable(model, theta >= 0)
    JuMP.@constraint(model, [j in 1:length(instance.graph.arcs)], theta >= sum(v[i, j] for i in 1:second_stage_paths_num))

    JuMP.@objective(
        model,
        Min,
        sum(
            instance.graph.arcs[j].cost.second*sum(v[i,j] for i in 1:second_stage_paths_num)
            for j in 1:length(instance.graph.arcs)
        )
        + instance.gamma*theta
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return createEmptyRrspSolution(length(instance.graph.arcs))
    end

    println("lambda: ", JuMP.value.(lambda))
    println("theta: ", JuMP.value(theta))
    for i in 1:second_stage_paths_num
        println("y_$i: ", JuMP.value.(y[i, :]))
    end

    # we have to decide among second stage path candidates =#
    second_stage_cost::Float64 = Inf
    second_stage_path_idx::Integer = 0
    for i in 1:second_stage_paths_num
        cost::Float64 = evaluateSecondStagePathContBudget(
            instance,
            RrspPath([JuMP.value(y[i,j]) > 0.5 for j in 1:length(instance.graph.arcs)])
        )
        if cost < second_stage_cost
            second_stage_cost = cost
            second_stage_path_idx = i
        end
    end
    second_stage_path::RrspPath = RrspPath([JuMP.value(y[second_stage_path_idx, j]) > 0.5 for j in 1:length(instance.graph.arcs)])

    return RrspSolution(x, second_stage_path, second_stage_cost)
end
