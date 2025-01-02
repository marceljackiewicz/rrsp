#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

module rrsp

import DataStructures

import Cbc
import JuMP

struct Node
    uid::Integer
end

struct Cost
    first::Float64
    second::Float64
    delta::Float64
end

struct Arc
    start_node::Node
    end_node::Node
    cost::Cost
end

struct Path
    arcs::Vector{Arc}
end

struct Graph
    nodes::Vector{Node}
    arcs::Vector{Arc}
end

struct RrspInstance
    graph::Graph
    s::Node
    t::Node
    k::Integer
    gamma::Float64
end

struct RrspSolution
    first_stage_path::Path
    second_stage_path::Path
    value::Float64
end

function getShortestPath(g::Graph, s::Node, t::Node)::Path
    optimizer = Cbc.Optimizer
    model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    # x[i] --- is the i-th arc taken?
    JuMP.@variable(model, x[i in 1:length(g.arcs)] >= 0)

    JuMP.@constraint(
        model,
        flow_constraint[i in 1:length(g.nodes); g.nodes[i] != s && g.nodes[i] != t],
        sum(x[j] for j in 1:length(g.arcs) if g.arcs[j].start_node == g.nodes[i]) == sum(x[j] for j in 1:length(g.arcs) if g.arcs[j].end_node == g.nodes[i])
    )

    JuMP.@constraint(
        model,
        flow_source,
          sum(x[i] for i in 1:length(g.arcs) if g.arcs[i].start_node == s)
        - sum(x[i] for i in 1:length(g.arcs) if g.arcs[i].end_node == s)
        == 1
    )
    JuMP.@constraint(
        model,
        flow_sink,
          sum(x[i] for i in 1:length(g.arcs) if g.arcs[i].start_node == t)
        - sum(x[i] for i in 1:length(g.arcs) if g.arcs[i].end_node == t)
        == -1
    )

    JuMP.@objective(model, Min, sum(x[i]*g.arcs[i].cost.first for i in 1:length(g.arcs)))

    # print(model)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        println("no primal solution")
        return Path([])
    end

    println("arcs: ", JuMP.value.(model[:x]))
    path::Path = Path([g.arcs[i] for i in 1:length(g.arcs) if JuMP.value(x[i]) > 0.5])
    return path
end

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
            push!(nodes, Node(uid_n1))
            last_node_id += 1
            uid_to_node_id[uid_n1] = last_node_id
        end
        start_node::Node = nodes[uid_to_node_id[uid_n1]]

        if (!haskey(uid_to_node_id, uid_n2))
            push!(nodes, Node(uid_n2))
            last_node_id += 1
            uid_to_node_id[uid_n2] = last_node_id
        end
        end_node::Node = nodes[uid_to_node_id[uid_n2]]

        return Arc(start_node, end_node, Cost(fst_stg_cst, snd_stg_cst, delta))
    end

    arcs::Vector{Arc} = []
    s::Node = Node(-1)
    t::Node = Node(-1)
    k::Integer = -1
    gamma::Float64 = 0.0
    open(file_name) do file
        params = split(readline(file))
        s = Node(parse(Int64, params[1]))
        t = Node(parse(Int64, params[2]))
        k = parse(Int64, params[3])
        gamma = parse(Float64, params[4])
        for line in eachline(file)
            push!(arcs, deserializeArc(line))
        end
    end

    return RrspInstance(Graph(nodes, arcs), s, t, k, gamma)
end

function solveRrspContBudgedDagForTheta(instance::RrspInstance, t::Integer)::RrspSolution
    # For now, inclusion neighbourhood only!
    optimizer = Cbc.Optimizer
    model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    theta::Float64 = 1/t

    # first stage path characteristic vector, x[i] --- is the i-th arc taken?
    JuMP.@variable(model, x[i in 1:length(instance.graph.arcs)], Bin)
    # second stage paths characteristic vectors, y[i, j] --- is the j-th arc taken in the i-th path?
    JuMP.@variable(model, y[i in 1:t, j in 1:length(instance.graph.arcs)], Bin)
    # second stage utility variables to model neighbourhood containment
    JuMP.@variable(model, z[i in 1:t, j in 1:length(instance.graph.arcs)] >= 0)
    # contiguous flow variables
    JuMP.@variable(model, 0 <= f_1[i in 1:length(instance.graph.arcs)] <= theta)
    JuMP.@variable(model, 0 <= f_2[i in 1:length(instance.graph.arcs)] <= 1 - theta)

    # constraint on second stage arc selection to correspond to optimal flow
    JuMP.@constraint(
        model,
        second_stage_opt_flow_to_selection[i in 1:length(instance.graph.arcs)],
        theta*sum(y[j, i] for j in 1:t) == f_1[i] + f_2[i]
    )

    # constraints for arcs in x to create a path
    JuMP.@constraint(
        model,
        first_stage_flow_constraint[i in 1:length(instance.graph.nodes); instance.graph.nodes[i] != instance.s && instance.graph.nodes[i] != instance.t],
        sum(x[j] for j in 1:length(instance.graph.arcs) if instance.graph.arcs[j].start_node == instance.graph.nodes[i]) == sum(x[j] for j in 1:length(instance.graph.arcs) if instance.graph.arcs[j].end_node == instance.graph.nodes[i])
    )
    JuMP.@constraint(
        model,
        first_stage_flow_source,
          sum(x[i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].start_node == instance.s)
        - sum(x[i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].end_node == instance.s)
        == 1
    )
    JuMP.@constraint(
        model,
        first_stage_flow_sink,
          sum(x[i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].start_node == instance.t)
        - sum(x[i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].end_node == instance.t)
        == -1
    )

    # constraints for arcs in y to create a path
    JuMP.@constraint(
        model,
        second_stage_flow_constraint[i in 1:length(instance.graph.nodes), j in 1:t; instance.graph.nodes[i] != instance.s && instance.graph.nodes[i] != instance.t],
        sum(y[j,k] for k in 1:length(instance.graph.arcs) if instance.graph.arcs[k].start_node == instance.graph.nodes[i]) == sum(y[j, k] for k in 1:length(instance.graph.arcs) if instance.graph.arcs[k].end_node == instance.graph.nodes[i])
    )
    JuMP.@constraint(
        model,
        second_stage_flow_source[j in 1:t],
          sum(y[j, i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].start_node == instance.s)
        - sum(y[j, i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].end_node == instance.s)
        == 1
    )
    JuMP.@constraint(
        model,
        second_stage_flow_sink[j in 1:t],
          sum(y[j, i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].start_node == instance.t)
        - sum(y[j, i] for i in 1:length(instance.graph.arcs) if instance.graph.arcs[i].end_node == instance.t)
        == -1
    )

    # arc inclusion neighbourhood constraints
    JuMP.@constraint(
        model,
        second_stage_different_arcs_bound[j in 1:t],
          sum(y[j, i] - z[j,i] for i in 1:length(instance.graph.arcs) ) <= instance.k
    )
    JuMP.@constraint(
        model,
        second_stage_z_vs_x[i in 1:t, j in 1:length(instance.graph.arcs)],
        z[i, j] <= x[j]
    )
    JuMP.@constraint(
        model,
        second_stage_z_vs_y[i in 1:t, j in 1:length(instance.graph.arcs)],
        z[i, j] <= y[i, j]
    )

    JuMP.@objective(
        model,
        Min,
          sum(x[i]*instance.graph.arcs[i].cost.first for i in 1:length(instance.graph.arcs))
        + sum(f_1[i]*instance.graph.arcs[i].cost.second for i in 1:length(instance.graph.arcs))
        + sum(f_2[i]*(instance.graph.arcs[i].cost.second + instance.graph.arcs[i].cost.delta) for i in 1:length(instance.graph.arcs))
        + instance.gamma*theta
    )

    # print(model)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    if (!JuMP.has_values(model))
        println("no primal solution")
        return RrspSolution(Path([]), Path([]), Inf)
    end

    println("arcs: ", JuMP.value.(model[:x]))
    first_stage_path::Path = Path([instance.graph.arcs[i] for i in 1:length(instance.graph.arcs) if JuMP.value(x[i]) > 0.5])
    # second_stage_path::Path = argmin(
    #     p::Path -> sum(instance.graph.arcs[i].cost.second*JuMP.value(f_1[i])),
    #     [Path([instance.graph.arcs[i] for i in 1:length(instance.graph.arcs) if JuMP.value(y[j, i]) > 0.5]) for j in 1:t]
    # ) # change Path to be a char. vector to get the proper solution
    second_stage_path::Path = Path([instance.graph.arcs[i] for i in 1:length(instance.graph.arcs) if JuMP.value(y[1, i]) > 0.5])

    return RrspSolution(first_stage_path, second_stage_path, JuMP.objective_value(model))
end

function getRrspContBudgetDag(instance::RrspInstance)::RrspSolution
    best_sol::RrspSolution = RrspSolution(Path([]), Path([]), Inf)
    t::Integer = length(instance.graph.arcs) + 1
    for theta::Integer in 1:t  # TODO 0::t, but refactor first
        fixed_theta_sol::RrspSolution = solveRrspContBudgedDagForTheta(instance, theta)
        if fixed_theta_sol.value < best_sol.value
            best_sol = fixed_theta_sol
        end
    end

    println("-----------")
    println(best_sol)
    return best_sol
end

end  # module rrsp
