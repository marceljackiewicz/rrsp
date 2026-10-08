"""
    layer_vertex_ids(H, W, layer) -> UnitRange{Int}

Vertex ids in layer `layer` of a wide layered digraph with `H` arcs per
``s``–``t`` path and intermediate width `W`. Layer `0` is `{1}` (the source) and layer `H` is
the sink. Intermediate layers `1, …, H-1` have `W` vertices.
"""
function layer_vertex_ids(H::Integer, W::Integer, layer::Integer)::UnitRange{Int}
    hh = Int(H)
    ww = Int(W)
    ll = Int(layer)
    (0 <= ll <= hh) || throw(ArgumentError("layer outside 0:H"))
    ww >= 1 || throw(ArgumentError("width < 1"))
    if ll == 0
        return 1:1
    elseif ll == hh
        n = 2 + (hh - 1) * ww
        return n:n
    end
    start = 2 + (ll - 1) * ww
    return start:(start + ww - 1)
end

"""
    gen_wide_layered(H, W) -> Graph

Layered digraph with `H` arcs per ``s``–``t`` path: layer 0 is the source, layer `H` is
the sink, and each intermediate layer has `W` vertices. Consecutive layers
are joined by a complete bipartite digraph. Every ``s``–``t`` path has `H` arcs.
"""
function gen_wide_layered(H::Integer, W::Integer)::Graph
    hh = Int(H)
    ww = Int(W)
    hh >= 1 || throw(ArgumentError("H < 1"))
    ww >= 1 || throw(ArgumentError("width < 1"))
    n = 2 + (hh - 1) * ww
    tail = Int32[]
    head = Int32[]
    for i in 0:(hh - 1)
        for u in layer_vertex_ids(hh, ww, i)
            for v in layer_vertex_ids(hh, ww, i + 1)
                push!(tail, Int32(u))
                push!(head, Int32(v))
            end
        end
    end
    return graph_new(n, tail, head)
end

"""
    is_st_layered(g, s, t) -> Bool

`true` if `g` is a DAG and every simple ``s``–``t`` path has the same number
of arcs (the layered case of the neighborhood-equivalence lemma).
"""
function is_st_layered(g::Graph, s::Integer, t::Integer)::Bool
    # All s–t paths have the same number of arcs iff the fewest arcs on an s–t
    # path equals the most arcs on one. The longest path is only computable in
    # polynomial time on a DAG, hence the acyclicity requirement (checked on
    # the whole graph, not only on the part reachable from `s`).
    is_dag(g) || return false
    ss = _vertex(g, s)
    tt = _vertex(g, t)
    # The only s–t path is the empty one.
    ss == tt && return true
    lo, hi = _st_arc_count_bounds(g, ss, tt)
    # `lo` is infinite when `t` is unreachable: there is no s–t path at all.
    return isfinite(lo) && lo == hi
end

# (minimum, maximum) number of arcs over the s–t paths of a DAG, as Float64;
# (Inf, Inf) if `t` is not reachable from `s`.
function _st_arc_count_bounds(g::Graph, s::Int32, t::Int32)::Tuple{Float64,Float64}
    n = Int(g.n)
    # Minimum: a shortest path with every arc of cost 1.
    w1 = ones(Float64, Int(g.m))
    st, _p, zmin = _shortest_path(g, s, t, w1)
    st != ST_OK && return (Inf, Inf)
    # Defensive: without a topological order the maximum is not computed
    # (never reached from `is_st_layered`, which has checked `is_dag`).
    order = _topo_order(g)
    order === nothing && return (zmin, Inf)
    # Maximum: longest-path DP in topological order. `mx[v]` is the largest
    # number of arcs of an s–v path, or -Inf while `v` is not known to be
    # reachable from `s`.
    mx = fill(-Inf, n)
    mx[Int(s)] = 0.0
    @inbounds for i in 1:n
        v = order[i]
        dv = mx[Int(v)]
        # Not reachable from `s`: nothing to relax.
        dv == -Inf && continue
        # All predecessors of `v` come earlier in the order, so `mx[v]` is final
        # here; extend every s–v path by one arc.
        for a in outgoing(g, v)
            j = Int(g.head[a])
            nd = dv + 1.0
            nd > mx[j] && (mx[j] = nd)
        end
    end
    return (zmin, mx[Int(t)])
end
