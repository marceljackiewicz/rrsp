#= Robust Recoverable Shortest Path Solver
#
#  Module implementing recoverable robust version of Shortest Path.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

function evaluateSecondStagePathContBudget(instance::RrspInstance, y::RrspPath)::Float64
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, Rrsp.optimizer)

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

"""
    solveRrspContBudget(instance::RrspInstance)::RrspSolution

Returns optimal first and second stage shortest ``s-t`` path in the graph `instance.graph` for Recoverable Robust Shortest Path
with continuous budgeted uncertainty.
The nodes ``s`` and ``t`` are given by `instance.s_idx`, `instance.t_idx` indices in `instance.graph.arcs` array.
The maximum budget value is equal to `instance.gamma`.

If the graph is a DAG, then use solveRrspContBudgetDag function.

The paths are computed using compact MIP model.
"""
function solveRrspContBudget(instance::RrspInstance)::RrspSolution
    model::JuMP.Model = JuMP.Model()
    JuMP.set_optimizer(model, Rrsp.optimizer)

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
    big_M::Float64 = 10e6  # TODO: find minimal value for big_M
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
        return createEmptyRrspSolution(length(instance.graph.arcs))
    end

    first_stage_path::RrspPath = RrspPath([JuMP.value(x[j]) > 0.5 for j in 1:length(instance.graph.arcs)])
    first_stage_cost::Float64 = sum(instance.graph.arcs[j].cost.first*first_stage_path.arcs[j] for j in 1:length(instance.graph.arcs))
    # now that the first stage path is known we have to decide among second stage path candidates =#
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

    return RrspSolution(first_stage_path, second_stage_path, first_stage_cost + second_stage_cost)
end

function solveRrspContBudgetDagForTheta(instance::RrspInstance, t::Integer)::RrspSolution
    model = JuMP.Model()
    JuMP.set_optimizer(model, Rrsp.optimizer)

    theta::Float64 = t > 0 ? 1/t : 0.0
    num_of_snd_stage_paths = t > 0 ? t : 1

    # first stage path characteristic vector, x[i] --- is the i-th arc taken?
    JuMP.@variable(model, x[i in 1:length(instance.graph.arcs)], Bin)
    # second stage paths characteristic vectors, y[i, j] --- is the j-th arc taken in the i-th path?
    JuMP.@variable(model, y[i in 1:num_of_snd_stage_paths, j in 1:length(instance.graph.arcs)], Bin)
    # contiguous flow variables
    JuMP.@variable(model, 0 <= f_1[i in 1:length(instance.graph.arcs)] <= theta)
    JuMP.@variable(model, 0 <= f_2[i in 1:length(instance.graph.arcs)] <= 1 - theta)

    # constraint on second stage arc selection to correspond to optimal flow
    JuMP.@constraint(
        model,
        second_stage_opt_flow_to_selection[i in 1:length(instance.graph.arcs)],
        (theta == 0.0 ? 1.0 : theta)*sum(y[j, i] for j in 1:num_of_snd_stage_paths) == f_1[i] + f_2[i]
    )

    # constraints for arcs in x to create a path
    addPathConstraints(model, model[:x], instance.graph, instance.s_idx, instance.t_idx)

    # constraints for arcs in y[i] to create a path
    for i in 1:num_of_snd_stage_paths
        addPathConstraints(model, model[:y][i, :], instance.graph, instance.s_idx, instance.t_idx)
    end

    # arc inclusion neighbourhood constraints for y[i]
    for i in 1:num_of_snd_stage_paths
        addNeighbourhoodConstraints(model, model[:x], model[:y][i, :], INCLUSION, instance.k, length(instance.graph.arcs))
    end

    JuMP.@objective(
        model,
        Min,
          sum(x[i]*instance.graph.arcs[i].cost.first for i in 1:length(instance.graph.arcs))
        + sum(f_1[i]*instance.graph.arcs[i].cost.second for i in 1:length(instance.graph.arcs))
        + sum(f_2[i]*(instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta) for i in 1:length(instance.graph.arcs))
        + instance.gamma*theta
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        println("no primal solution", t)
        return createEmptyRrspSolution(length(instance.graph.arcs))
    end

    first_stage_path::RrspPath = RrspPath([JuMP.value(x[i]) > 0.5 for i in 1:length(instance.graph.arcs)])
    second_stage_path::RrspPath = argmin(
        p::RrspPath -> sum(
            p.arcs[i]*(
                  instance.graph.arcs[i].cost.second*JuMP.value(f_1[i])
                + (instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta)*JuMP.value(f_2[i]))
            for i in 1:length(instance.graph.arcs)),
        [RrspPath([JuMP.value(y[j, i]) > 0.5 for i in 1:length(instance.graph.arcs)]) for j in 1:num_of_snd_stage_paths]
    )

    return RrspSolution(first_stage_path, second_stage_path, JuMP.objective_value(model))
end


"""
    solveRrspContBudgetDag(instance::RrspInstance)::RrspSolution

Returns optimal first and second stage shortest ``s-t`` path in the graph `instance.graph` for Recoverable Robust Shortest Path
with continuous budgeted uncertainty when `instance.graph` is a directed acyclic graph.
The nodes ``s`` and ``t`` are given by `instance.s_idx`, `instance.t_idx` indices in `instance.graph.arcs` array.
The maximum budget value is equal to `instance.gamma`.

The paths are computed by solving a family of ``m + 1`` compact MIP models, where ``m`` is the number of arcs in the DAG.
Since the graph is acyclic, the model does not contain anti-cyclic constraints, unlike the model for general digraphs.
"""
function solveRrspContBudgetDag(instance::RrspInstance)::RrspSolution
    best_sol::RrspSolution = createEmptyRrspSolution(length(instance.graph.arcs))
    t::Integer = length(instance.graph.arcs) + 1
    for theta::Integer in 0:t
        fixed_theta_sol::RrspSolution = solveRrspContBudgetDagForTheta(instance, theta)
        if fixed_theta_sol.value < best_sol.value
            best_sol = fixed_theta_sol
        end
    end

    return best_sol
end

function solveRrspContBudgetAsp(instance::RrspInstance)::RrspSolution
    # TODO: not implemented
    tree::AspTree = getAspTreeDecomposition(instance.graph)

    function getLeftNode(tree::AspTree, idx::Integer)
        return tree.nodes[tree.nodes[i].left]
    end

    function getRightNode(tree::AspTree, idx::Integer)
        return tree.nodes[tree.nodes[i].right]
    end

    while length(tree.nodes) > 1
        comp_idx::Integer = 0
        for i in 1:length(tree.nodes)
            if !tree.nodes[i].is_leaf && getLeftNode(tree, i).is_leaf && getRightNode(tree, i).is_leaf
                comp_idx = i
                break
            end
        end

        
    end
end
