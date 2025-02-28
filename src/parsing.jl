#= Robust Recoverable Shortest Path Solver
#
#  Module implementing parsing the problem instances.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

"""
    parseInstanceFromFile(file_name::String)::RrspInstance

Return serialized RRSP instance given in `file_name`.
The input file should conform to the following format:

    s t N k G
    <arc>
    ...
    <arc>

where
- `s` -- path start node id,
- `t` -- path end node id,
- `N` -- neighbourhood type; one of `INC`, `EXC`, `SYM_DIFF`,
- `k` -- neighbourhood size parameter,
- `G` -- budget,
and `<arc>` is a line serializing one arc:

    n1 n2 c1 c2 d

where
- `n1` -- arc start node id,
- `n2` -- arc end node id,
- `c1` -- first stage cost,
- `c2` -- second stage cost lower bound,
- `d` -- second stage cost maximum deviation.

Example:

    1 3 INC 1 100.0
    1 2 0.0 50.0 10.0
    1 2 0.0 0.0 100.0
    2 3 0.0 20.0 30.0
    2 3 0.0 40.0 0.0
"""
function parseInstanceFromFile(file_name::String)::RrspInstance
    # TODO: add serialized instance validation
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
        neighbourhood = stringToNeighbourhoodType(params[3])
        k = parse(Int64, params[4])
        gamma = parse(Float64, params[5])
        for line in eachline(file)
            push!(arcs, deserializeArc(line))
        end
        # uid of s and t nodes can only be mapped to indices after the arcs are parsed
        s_idx = uid_to_node_id[parse(Int64, params[1])]
        t_idx = uid_to_node_id[parse(Int64, params[2])]
    end

    return RrspInstance(Graph(nodes, arcs), s_idx, t_idx, neighbourhood, k, gamma)
end
