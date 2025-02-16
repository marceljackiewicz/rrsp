#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

module rrsp

import DataStructures

import Cbc
import JuMP

include("types.jl")

include("deterministic.jl")
include("incremental.jl")
include("recoverable.jl")
include("recoverable_robust.jl")

# TODO: add serialized instance validation
function parseInstanceFromFile(file_name::String)::RrspInstance
    nodes::Vector{Node} = []
    uid_to_node_id::Dict{Integer, Integer} = Dict{Integer, Integer}()
    last_node_id::Integer = 0

    function deserializeArc(serialized_arc::String)::Arc
        n1, n2, fc, sc, dc = split(serialized_arc)
        uid_n1::Integer = parse(Int64, n1)
        uid_n2::Integer = parse(Int64, n2)
        fst_stg_cst::Float64 = parse(Float64, fc)
        snd_stg_cst::Float64 = parse(Float64, sc)
        delta::Float64 = parse(Float64, dc)

        if (!haskey(uid_to_node_id, uid_n1))
            push!(nodes, Node(uid_n1, length(nodes) + 1))
            last_node_id += 1
            uid_to_node_id[uid_n1] = last_node_id
        end
        start_node::Node = nodes[uid_to_node_id[uid_n1]]

        if (!haskey(uid_to_node_id, uid_n2))
            push!(nodes, Node(uid_n2, length(nodes) + 1))
            last_node_id += 1
            uid_to_node_id[uid_n2] = last_node_id
        end
        end_node::Node = nodes[uid_to_node_id[uid_n2]]

        return Arc(start_node, end_node, Cost(fst_stg_cst, snd_stg_cst, delta))
    end

    arcs::Vector{Arc} = []
    s_idx::Integer = -1
    t_idx::Integer = -1
    neighbourhood::NeighbourhoodType = NEIGHBOURHOOD_NOT_SET
    k::Integer = -1
    gamma::Float64 = 0.0
    open(file_name) do file
        params = split(readline(file))
        s_idx = parse(Int64, params[1])
        t_idx = parse(Int64, params[2])
        neighbourhood = stringToNeighbourhoodType(params[3])
        k = parse(Int64, params[4])
        gamma = parse(Float64, params[5])
        for line in eachline(file)
            push!(arcs, deserializeArc(line))
        end
    end

    return RrspInstance(Graph(nodes, arcs), s_idx, t_idx, neighbourhood, k, gamma)
end

function solveRrspContBudgetDagForTheta(instance::RrspInstance, t::Integer)::RrspSolution
    optimizer = Cbc.Optimizer
    model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

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
        return RrspSolution(Path([]), Path([]), Inf)
    end

    first_stage_path::Path = Path([JuMP.value(x[i]) > 0.5 for i in 1:length(instance.graph.arcs)])
    second_stage_path::Path = argmin(
        p::Path -> sum(
            p.arcs[i]*(
                  instance.graph.arcs[i].cost.second*JuMP.value(f_1[i])
                + (instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta)*JuMP.value(f_2[i]))
            for i in 1:length(instance.graph.arcs)),
        [Path([JuMP.value(y[j, i]) > 0.5 for i in 1:length(instance.graph.arcs)]) for j in 1:num_of_snd_stage_paths]
    )

    # println("solution ", t, " ", RrspSolution(first_stage_path, second_stage_path, JuMP.objective_value(model)))
    return RrspSolution(first_stage_path, second_stage_path, JuMP.objective_value(model))
end

function solveRrspContBudgetDag(instance::RrspInstance)::RrspSolution
    best_sol::RrspSolution = RrspSolution(Path([]), Path([]), Inf)
    t::Integer = length(instance.graph.arcs) + 1
    for theta::Integer in 0:t
        fixed_theta_sol::RrspSolution = solveRrspContBudgetDagForTheta(instance, theta)
        if fixed_theta_sol.value < best_sol.value
            best_sol = fixed_theta_sol
        end
    end

    return best_sol
end

function getAspTreeDecomposition(g::Graph)::AspTree
    function areParallelComponents(a::AspComponent, b::AspComponent)::Bool
        return (a.s == b.s) && (a.t == b.t)
    end

    function findParallelComponent(components::Vector{Array{AspComponent}})::Array{Integer}
        for idx1 in 1:(length(components) - 1)  # the last arc can't have a pair
            if length(components[idx1]) == 0
                continue
            end

            for idx2 in (idx1 + 1):length(components)
                if length(components[idx2]) == 0
                    continue
                end
    
                if areParallelComponents(components[idx1][1], components[idx2][1])
                    return [idx1, idx2]
                end
            end
        end

        return [0, 0]
    end

    function areSeriesComponents(a::AspComponent, b::AspComponent)::Bool
        return a.t == b.s
    end

    function findSeriesComponent(components::Vector{Array{AspComponent}})::Array{Integer}
        arcs_in::Dict{Node, Vector{Integer}} = Dict{Node, Vector{Integer}}()
        arcs_out::Dict{Node, Vector{Integer}} = Dict{Node, Vector{Integer}}()

        for idx in 1:length(components)
            if length(components[idx]) > 0
                haskey(arcs_in, components[idx][1].t) ? push!(arcs_in[components[idx][1].t], idx) : arcs_in[components[idx][1].t] = [idx]
                haskey(arcs_out, components[idx][1].s) ? push!(arcs_out[components[idx][1].s], idx) : arcs_out[components[idx][1].s] = [idx]
            end
        end

        for item in arcs_in
            node::Node = item.first
            if haskey(arcs_out, node) && length(arcs_in[node]) == length(arcs_out[node]) == 1
                return [arcs_in[node][1], arcs_out[node][1]]
            end
        end

        return [0, 0]
    end

    # create an ASP component for each arc in g
    components::Vector{Array{AspComponent}} = [[AspComponent(arc.start_node, arc.end_node)] for arc in g.arcs]
    tree_nodes::Vector{AspTreeNode} = [AspTreeNode(arc.start_node, arc.end_node, 0, 0, [arc], NONE, true) for arc in g.arcs]

    components_removed = 0
    while components_removed < length(g.arcs) - 1

        # while there are two series components
        s_comps = findSeriesComponent(components)
        while s_comps != [0, 0]
            push!(components, [AspComponent(components[s_comps[1]][1].s, components[s_comps[2]][1].t)])
            push!(tree_nodes, AspTreeNode(components[s_comps[1]][1].s, components[s_comps[2]][1].t, s_comps[1], s_comps[2], [], SERIES, false))
            components_removed += 1
            components[s_comps[2]] = []
            components[s_comps[1]] = []
            s_comps = findSeriesComponent(components)
        end

        # while there are two components parallel components
        p_comps = findParallelComponent(components)
        while p_comps != [0, 0]
            push!(components, [AspComponent(components[p_comps[1]][1].s, components[p_comps[1]][1].t)])
            push!(tree_nodes, AspTreeNode(components[p_comps[1]][1].s, components[p_comps[1]][1].t, p_comps[1], p_comps[2], [], PARALLEL, false))
            components_removed += 1
            components[p_comps[2]] = []
            components[p_comps[1]] = []
            p_comps = findParallelComponent(components)
        end
    end

    return AspTree(length(tree_nodes), tree_nodes)
end

function solveRrspContBudgetAsp(instance::RrspInstance)::RrspSolution
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

end  # module rrsp
