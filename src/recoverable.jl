""" Robust Recoverable Shortest Path Solver

Module implementing recoverable version of Shortest Path.

Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
"""

include("tree_decomposition.jl")

"""
    solveRecoverableShortestPath(instance::RrspInstance)::RrspSolution

Returns optimal first and second stage shortest ``s-t`` path in the graph `instance.graph` for Recoverable Shortest Path,
to which the RRSP problem with interval uncertainty reduces to.
The nodes ``s`` and ``t`` are given by `instance.s_idx`, `instance.t_idx` indices in `instance.graph.arcs` array.

The paths are computed using compact MIP model.
"""
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

"""
    solveRecoverableShortestPathInAsp(instance::RrspInstance)::RrspSolution

Returns optimal first and second stage shortest ``s-t`` path in the graph `instance.graph`
for Recoverable Shortest Path in Arc Series-Parallel (ASP) graph.

The `instance.graph` is assumed to be ASP. Hence, it's not validated further.
If this is not the case, then the behaviour and the return value of this function are undefined.

The paths are computed using combinatorial algoritm.

NOTE: The only supported neighbourhood is **inclusion** for now!
"""
function solveRecoverableShortestPathInAsp(instance::RrspInstance)::RrspSolution
    @assert instance.neighbourhood == INCLUSION ["solveRecoverableShortestPathInAsp; exclusion and sym diff not implemented!"]
    tree::AspTree = getAspTreeDecomposition(instance.graph)

    # reserve storage for all data for nodes composition
    asp_nodes_data::Vector{AspNodeData} = [
        AspNodeData(
            createEmptyPath(length(instance.graph.arcs)),
            [createEmptyPath(length(instance.graph.arcs)) for _ in 1:instance.k],
            [createEmptyRrspSolution(length(instance.graph.arcs)) for _ in 1:(1 + instance.k)]
        ) for _ in 1:length(tree.nodes)
    ]

    function getPathCardinality(p::Path)::Integer
        return sum(p.arcs)
    end

    function getPathCostFirstStage(p::Path)::Float64
        if getPathCardinality(p) == 0
            return Inf
        end

        return sum(instance.graph.arcs[i].cost.first for i in 1:length(p.arcs) if p.arcs[i])        
    end

    function getArcCostSecondStageUpperBound(cost::Cost)
        return cost.second + cost.delta
    end

    function getPathCostSecondStageUpperBound(p::Path)::Float64
        if getPathCardinality(p) == 0
            return Inf
        end

        return sum(getArcCostSecondStageUpperBound(instance.graph.arcs[i].cost) for i in 1:length(p.arcs) if p.arcs[i])
    end

    function getArcCostUpperBound(cost::Cost)
        return cost.first + getArcCostSecondStageUpperBound(cost)
    end

    # initialize data for leaves
    for node_idx in 1:length(tree.nodes)
        if tree.nodes[node_idx].is_leaf
            # set first stage path
            asp_nodes_data[node_idx].opt_first_stage_path.arcs[tree.nodes[node_idx].arc_idx] = 1
            # set second stage path
            asp_nodes_data[node_idx].opt_second_stage_path[1].arcs[tree.nodes[node_idx].arc_idx] = 1
            # set optimal pair
            asp_nodes_data[node_idx].opt_solution_paths[1].first_stage_path.arcs[tree.nodes[node_idx].arc_idx] = 1
            asp_nodes_data[node_idx].opt_solution_paths[1].second_stage_path.arcs[tree.nodes[node_idx].arc_idx] = 1
            asp_nodes_data[node_idx].opt_solution_paths[1].value = getArcCostUpperBound(instance.graph.arcs[tree.nodes[node_idx].arc_idx].cost)
        end
    end

    function isNodeValid(node::AspTreeNode)::Bool
        # either leaf and operation NONE or not a leaf and a valid operation
        return (node.is_leaf && node.operation == NONE) || (!node.is_leaf && node.operation != NONE)
    end

    function invalidNode(node::AspTreeNode)::Nothing
        node.is_leaf = false
        node.operation = NONE
        return
    end

    function hasNodeTwoLeaves(node::AspTreeNode)::Bool
        return isNodeValid(node) && !node.is_leaf && tree.nodes[node.left].is_leaf && tree.nodes[node.right].is_leaf
    end

    function mergeTwoPaths(p1::Path, p2::Path)::Path
        return Path([p1.arcs[idx] || p2.arcs[idx] for idx in 1:length(p1.arcs)])
    end

    function mergeTwoRrspSolutions(s1::RrspSolution, s2::RrspSolution)::RrspSolution
        return RrspSolution(
            mergeTwoPaths(s1.first_stage_path, s2.first_stage_path),
            mergeTwoPaths(s1.second_stage_path, s2.second_stage_path),
            s1.value + s2.value
        )
    end

    function performParallelComposition(node::AspNodeData, left::AspNodeData, right::AspNodeData)::Nothing
        # opt first stage path
        node.opt_first_stage_path = argmin(p::Path -> getPathCostFirstStage(p), [left.opt_first_stage_path, right.opt_first_stage_path])

        # opt second stage paths
        for l in 1:instance.k
            node.opt_second_stage_path[l] = argmin(p::Path -> getPathCostSecondStageUpperBound(p), [left.opt_second_stage_path[l], right.opt_second_stage_path[l]])
        end

        # opt feasible pairs
        node.opt_solution_paths[1] = argmin(s::RrspSolution -> s.value, [left.opt_solution_paths[1], right.opt_solution_paths[1]])
        for l in 2:(instance.k + 1)
            node.opt_solution_paths[l] = argmin(s::RrspSolution -> s.value, [
                left.opt_solution_paths[l],
                right.opt_solution_paths[l],
                RrspSolution(
                    left.opt_first_stage_path,
                    right.opt_second_stage_path[l - 1],
                    getPathCostFirstStage(left.opt_first_stage_path) + getPathCostSecondStageUpperBound(right.opt_second_stage_path[l - 1])
                ),
                RrspSolution(
                    right.opt_first_stage_path,
                    left.opt_second_stage_path[l - 1],
                    getPathCostFirstStage(right.opt_first_stage_path) + getPathCostSecondStageUpperBound(left.opt_second_stage_path[l - 1])
                ),
            ])
        end
    end

    function performSeriesComposition(node::AspNodeData, left::AspNodeData, right::AspNodeData)
        # opt first stage path
        node.opt_first_stage_path = mergeTwoPaths(left.opt_first_stage_path, right.opt_first_stage_path)

        # opt second stage paths
        node.opt_second_stage_path[1] = createEmptyPath(length(instance.graph.arcs))
        for l in 2:instance.k
            node.opt_second_stage_path[l] = argmin(
                p::Path -> getPathCostSecondStageUpperBound(p),
                [mergeTwoPaths(left.opt_second_stage_path[j], right.opt_second_stage_path[l - j]) for j in 1:(l - 1)]
            )
        end

        # opt feasible pairs
        node.opt_solution_paths[1] = mergeTwoRrspSolutions(left.opt_solution_paths[1], right.opt_solution_paths[1])
        for l in 2:(1 + instance.k)
            node.opt_solution_paths[l] = argmin(
                s::RrspSolution -> s.value,
                [mergeTwoRrspSolutions(left.opt_solution_paths[j], right.opt_solution_paths[l - j + 1]) for j in 1:l]
            )
        end
    end

    # perform leaves composition until root is the only tree node
    leaves_counter = length(instance.graph.arcs)
    while leaves_counter > 1
        for node_idx in 1:length(tree.nodes)
            # perform composition iff the node has two leaves as children
            if hasNodeTwoLeaves(tree.nodes[node_idx])
                if tree.nodes[node_idx].operation == PARALLEL
                    performParallelComposition(asp_nodes_data[node_idx], asp_nodes_data[tree.nodes[node_idx].left], asp_nodes_data[tree.nodes[node_idx].right])
                elseif tree.nodes[node_idx].operation == SERIES
                    performSeriesComposition(asp_nodes_data[node_idx], asp_nodes_data[tree.nodes[node_idx].left], asp_nodes_data[tree.nodes[node_idx].right])
                end

                # remove leaves from the tree
                invalidNode(tree.nodes[tree.nodes[node_idx].left])  # just invalid the nodes
                invalidNode(tree.nodes[tree.nodes[node_idx].right])
                
                # the node becomes a leaf then
                tree.nodes[node_idx].is_leaf = true
                tree.nodes[node_idx].operation = NONE

                leaves_counter -= 1
            end
        end
    end

    return argmin(s::RrspSolution -> s.value, asp_nodes_data[tree.root_idx].opt_solution_paths)
end
