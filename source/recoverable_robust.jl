#= Robust Recoverable Shortest Path Solver
#
#  Module implementing recoverable robust version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

function evaluateSecondStagePathContBudget(instance::RrspInstance, y::Path)::Float64
    optimizer = Cbc.Optimizer
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    JuMP.@variable(model, u[i in 1:length(y.arcs)] >= 0)
    JuMP.@constraint(model, [i in 1:length(y.arcs)], u[i] <= instance.graph.arcs[i].cost.delta)
    JuMP.@constraint(model, sum(u) <= instance.gamma)

    JuMP.@variable(model, t >= 0)
    JuMP.@constraint(model, t <= sum(y.arcs[i]*(instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta) for i in 1:length(y.arcs)))

    JuMP.@objective(model, Max, t)

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return Inf
    end

    return JuMP.value(t)
end

#= Returns first and second stage s-t paths and the optimal objective function value
#  for RRSP with continuous budget in general digraphs.
#  Use another function for DAGs!
#
#  The paths are computed using compact ILP model.
=#
function solveRrspContBudget(instance::RrspInstance)::RrspSolution
    optimizer = Cbc.Optimizer
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    second_stage_paths_num::Integer = length(instance.graph.arcs) + 1

    # x[i] --- is the i-th arc in the first stage path taken?
    JuMP.@variable(model, x[i in 1:length(instance.graph.arcs)], Bin)
    # second stage paths characteristic vectors, y[i, j] --- is the j-th arc taken in the i-th path?
    JuMP.@variable(model, y[i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)], Bin)

    addPathConstraints(model, model[:x], instance.graph, instance.s_idx, instance.t_idx)
    for i in 1:second_stage_paths_num
        # constraints for arcs in y[i] to create a path
        addPathConstraints(model, model[:y][i, :], instance.graph, instance.s_idx, instance.t_idx)
        # neighbourhood constraints for y[i]
        addNeighbourhoodConstraints(model, model[:x], model[:y][i, :], instance.neighbourhood, instance.k, length(instance.graph.arcs))
    end

    #= if the graph contains cycles, the anti-cyclic constraints must be used for the first stage path,
    as taking more arcs may prove beneficial in the second stage =#
    addAntiCyclicConstraints(model, model[:x], instance.graph, instance.s_idx, instance.t_idx)
    # the anti-cyclic constraints can be droped for the second stage path iff inclusion neighbourhood is used
    if instance.neighbourhood != INCLUSION
        for i in 1:second_stage_paths_num
            addAntiCyclicConstraints(model, model[:y], instance.graph, instance.s_idx, instance.t_idx)
        end
    end

    #= the sum of flows on each y[i] (evaluation) path is 1.0 =#
    JuMP.@variable(model, lambda[i in 1:second_stage_paths_num] >= 0)
    JuMP.@constraint(model, sum(lambda) == 1.0)

    #= some budget must be spent to send some flow with y[i] path over an arc =#
    JuMP.@variable(model, gamma[j in 1:length(instance.graph.arcs)] >= 0)
    JuMP.@variable(model, theta >= 0)

    #= v[i,j] linearizes (to represent) the product of y[i,j]*lambda[i] =#
    JuMP.@variable(model, v[i in 1:second_stage_paths_num, j in 1:length(instance.graph.arcs)] >= 0)
    big_M::Float64 = 10e6  # ???
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

    JuMP.@constraint(model, [j in 1:length(instance.graph.arcs)], gamma[j] + theta >= sum(v[i, j] for i in 1:second_stage_paths_num))

    JuMP.@objective(
        model,
        Min,
        sum(
              instance.graph.arcs[j].cost.first*x[j]
            + instance.graph.arcs[j].cost.second*sum(v[i,j] for i in 1:second_stage_paths_num)
            + instance.graph.arcs[j].cost.delta*gamma[j]
            for j in 1:length(instance.graph.arcs)
        )
        + instance.gamma*theta
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        return RrspSolution(Path([]), Path([]), Inf)
    end

    first_stage_path::Path = Path([JuMP.value(x[j]) > 0.5 for j in 1:length(instance.graph.arcs)])
    first_stage_cost::Float64 = sum(instance.graph.arcs[j].cost.first*first_stage_path.arcs[j] for j in 1:length(instance.graph.arcs))
    # now that the first stage path is known we have to decide among second stage path candidates =#
    second_stage_cost::Float64 = Inf
    second_stage_path_idx::Integer = 0
    for i in 1:second_stage_paths_num
        cost::Float64 = evaluateSecondStagePathContBudget(
            instance,
            Path([JuMP.value(y[i,j]) > 0.5 for j in 1:length(instance.graph.arcs)])
        )
        if cost < second_stage_cost
            second_stage_cost = cost
            second_stage_path_idx = i
        end
    end
    second_stage_path::Path = Path([JuMP.value(y[second_stage_path_idx, j]) > 0.5 for j in 1:length(instance.graph.arcs)])

    return RrspSolution(first_stage_path, second_stage_path, first_stage_cost + second_stage_cost)
end
