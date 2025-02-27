#= Robust Recoverable Shortest Path Solver
#
#  Module implementing binary tree decomposition of arc series-parallel graphs.
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

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
    tree_nodes::Vector{AspTreeNode} = [AspTreeNode(g.arcs[arc_idx].start_node, g.arcs[arc_idx].end_node, 0, 0, arc_idx, NONE, true) for arc_idx in 1:length(g.arcs)]

    components_removed = 0
    while components_removed < length(g.arcs) - 1

        # while there are two series components
        s_comps = findSeriesComponent(components)
        while s_comps != [0, 0]
            push!(components, [AspComponent(components[s_comps[1]][1].s, components[s_comps[2]][1].t)])
            push!(tree_nodes, AspTreeNode(components[s_comps[1]][1].s, components[s_comps[2]][1].t, s_comps[1], s_comps[2], 0, SERIES, false))
            components_removed += 1
            components[s_comps[2]] = []
            components[s_comps[1]] = []
            s_comps = findSeriesComponent(components)
        end

        # while there are two components parallel components
        p_comps = findParallelComponent(components)
        while p_comps != [0, 0]
            push!(components, [AspComponent(components[p_comps[1]][1].s, components[p_comps[1]][1].t)])
            push!(tree_nodes, AspTreeNode(components[p_comps[1]][1].s, components[p_comps[1]][1].t, p_comps[1], p_comps[2], 0, PARALLEL, false))
            components_removed += 1
            components[p_comps[2]] = []
            components[p_comps[1]] = []
            p_comps = findParallelComponent(components)
        end
    end

    return AspTree(length(tree_nodes), tree_nodes)
end
