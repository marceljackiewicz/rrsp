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
end

struct RrspSolution
    first_stage_path::Path
    second_stage_path::Path
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
    open(file_name) do file
        s, t = map(x -> Node(parse(Int64, x)), split(readline(file)))
        for line in eachline(file)
            push!(arcs, deserializeArc(line))
        end
    end

    return RrspInstance(Graph(nodes, arcs), s, t)
end

end  # module rrsp
