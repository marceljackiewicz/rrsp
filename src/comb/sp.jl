function _check_weights(w::AbstractVector{<:Real}, m::Int)
    length(w) == m || throw(ArgumentError("weight length mismatch"))
    @inbounds for a in 1:m
        Float64(w[a]) < 0 && throw(ArgumentError("weight of arc $a is negative"))
    end
    return nothing
end

function _copy_weights(w::AbstractVector{<:Real}, m::Int)::Vector{Float64}
    _check_weights(w, m)
    out = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        out[a] = Float64(w[a])
    end
    return out
end

function _topo_order(g::Graph)::Union{Nothing,Vector{Int32}}
    n = Int(g.n)
    @inbounds for a in 1:Int(g.m)
        if g.tail[a] == g.head[a]
            return nothing
        end
    end
    indeg = Vector{Int32}(undef, n)
    @inbounds for v in 1:n
        indeg[v] = in_degree(g, Int32(v))
    end
    order = Vector{Int32}(undef, n)
    queue = Vector{Int32}(undef, n)
    qh = 1
    qt = 0
    @inbounds for v in Int32(1):g.n
        if indeg[v] == 0
            qt += 1
            queue[qt] = v
        end
    end
    seen = 0
    while qh <= qt
        v = queue[qh]
        qh += 1
        seen += 1
        order[seen] = v
        @inbounds for a in outgoing(g, v)
            w = g.head[a]
            indeg[w] -= Int32(1)
            if indeg[w] == 0
                qt += 1
                queue[qt] = w
            end
        end
    end
    seen == n || return nothing
    return order
end

function _path_from_pred(
    g::Graph,
    s::Int32,
    t::Int32,
    pred::AbstractVector{Int32},
    dist::AbstractVector{Float64},
)::Tuple{Status,Path,Float64}
    m = Int(g.m)
    if s == t
        return (ST_OK, path_empty(m), 0.0)
    end
    dt = dist[Int(t)]
    if dt == Inf || pred[Int(t)] == Int32(0)
        return (ST_INFEASIBLE, path_empty(m), Inf)
    end
    len = 0
    v = t
    while v != s
        a = pred[Int(v)]
        if a == Int32(0)
            return (ST_INFEASIBLE, path_empty(m), Inf)
        end
        len += 1
        len > m && return (ST_ERROR, path_empty(m), Inf)
        v = g.tail[a]
    end
    seq = Vector{Int32}(undef, len)
    v = t
    @inbounds for i in len:-1:1
        a = pred[Int(v)]
        seq[i] = a
        v = g.tail[a]
    end
    return (ST_OK, path_from_seq(g, s, t, seq), dt)
end

# Binary min-heap of (distance, vertex), ordered lexicographically so that ties
# in distance are extracted by increasing vertex index.
@inline _heap_less(d1::Float64, v1::Int32, d2::Float64, v2::Int32) = d1 < d2 || (d1 == d2 && v1 < v2)

function _heap_push!(hd::Vector{Float64}, hv::Vector{Int32}, d::Float64, v::Int32)
    push!(hd, d)
    push!(hv, v)
    i = length(hd)
    @inbounds while i > 1
        p = i >> 1
        _heap_less(hd[i], hv[i], hd[p], hv[p]) || break
        hd[i], hd[p] = hd[p], hd[i]
        hv[i], hv[p] = hv[p], hv[i]
        i = p
    end
    return nothing
end

function _heap_pop!(hd::Vector{Float64}, hv::Vector{Int32})::Tuple{Float64,Int32}
    d0, v0 = hd[1], hv[1]
    n = length(hd)
    hd[1], hv[1] = hd[n], hv[n]
    pop!(hd)
    pop!(hv)
    n -= 1
    i = 1
    @inbounds while true
        l = 2i
        l > n && break
        c = (l < n && _heap_less(hd[l + 1], hv[l + 1], hd[l], hv[l])) ? l + 1 : l
        _heap_less(hd[c], hv[c], hd[i], hv[i]) || break
        hd[i], hd[c] = hd[c], hd[i]
        hv[i], hv[c] = hv[c], hv[i]
        i = c
    end
    return d0, v0
end

"""
    _shortest_path_dijkstra(g, s, t, w) -> (status, path, z)

Nonnegative Dijkstra on the forward star of `g`, with a binary heap and lazy
deletion (``O((n + m) \\log n)``). Each vertex is settled at most once, in order of
distance and then vertex index. Ties keep the first predecessor written (`<`,
not `≤`).
"""
function _shortest_path_dijkstra(
    g::Graph,
    s::Int32,
    t::Int32,
    w::Vector{Float64},
)::Tuple{Status,Path,Float64}
    n = Int(g.n)
    dist = fill(Inf, n)
    pred = zeros(Int32, n)
    done = falses(n)
    dist[Int(s)] = 0.0
    hd = Float64[0.0]
    hv = Int32[s]
    while !isempty(hd)
        du, u = _heap_pop!(hd, hv)
        @inbounds begin
            done[u] && continue
            done[u] = true
            u == t && break
            for a in outgoing(g, u)
                vv = g.head[a]
                nd = du + w[a]
                if nd < dist[vv]
                    dist[vv] = nd
                    pred[vv] = a
                    _heap_push!(hd, hv, nd, vv)
                end
            end
        end
    end
    return _path_from_pred(g, s, t, pred, dist)
end

function _shortest_path_on_order(
    g::Graph,
    s::Int32,
    t::Int32,
    w::Vector{Float64},
    order::Vector{Int32},
)::Tuple{Status,Path,Float64}
    n = Int(g.n)
    dist = Vector{Float64}(undef, n)
    pred = Vector{Int32}(undef, n)
    @inbounds for v in 1:n
        dist[v] = Inf
        pred[v] = Int32(0)
    end
    dist[Int(s)] = 0.0
    @inbounds for i in 1:n
        u = order[i]
        du = dist[Int(u)]
        du == Inf && continue
        for a in outgoing(g, u)
            vv = Int(g.head[a])
            nd = du + w[a]
            if nd < dist[vv]
                dist[vv] = nd
                pred[vv] = a
            end
        end
    end
    return _path_from_pred(g, s, t, pred, dist)
end

"""
    _shortest_path_dag(g, s, t, w) -> (status, path, z)

Single-source shortest path on a DAG, in topological order.

# Throws
- `ArgumentError`: if `g` is not a DAG.
"""
function _shortest_path_dag(
    g::Graph,
    s::Int32,
    t::Int32,
    w::Vector{Float64},
)::Tuple{Status,Path,Float64}
    order = _topo_order(g)
    order === nothing && throw(ArgumentError("graph is not a DAG"))
    return _shortest_path_on_order(g, s, t, w, order)
end

function _shortest_path(
    g::Graph,
    s::Int32,
    t::Int32,
    w::Vector{Float64},
)::Tuple{Status,Path,Float64}
    m = Int(g.m)
    if s == t
        return (ST_OK, path_empty(m), 0.0)
    end
    order = _topo_order(g)
    if order === nothing
        return _shortest_path_dijkstra(g, s, t, w)
    end
    return _shortest_path_on_order(g, s, t, w, order)
end

function _finish_sp(
    net::Network,
    status::Status,
    path::Path,
    z::Float64,
    t0::UInt64;
    rob::Bool = false,
    method::Method = METHOD_COMB,
    gap::Float64 = 0.0,
)::Solution
    m = Int(net.graph.m)
    dt = (time_ns() - t0) / 1e9
    if status != ST_OK
        return solution_empty(m; status = status, method = method, time_sec = dt)
    end
    z1 = path_cost(path, net.costs.C)
    z2 = rob ? z : 0.0
    return Solution(path, path_empty(m), z, z1, z2, ST_OK, method, dt, gap)
end

"""
    solve_sp(net, w, solver) -> Solution

Shortest ``s``–``t`` path in `net` under arc weights `w`.

The path is returned in `first`. Combinatorial methods require nonnegative
weights: Dijkstra on the forward star, or a topological scan when the graph
is a DAG. `METHOD_AUTO` and `METHOD_COMB` both use this algorithm.
`METHOD_MIP` uses a compact path formulation when `solver.optimizer` is set,
and returns `ST_NOT_IMPL` if it is `nothing`.

If ``s = t``, the result is the trivial (empty) path of cost `0`.

# Throws
- `ArgumentError`: if `w` has length other than `m`, or some weight is negative.
"""
function solve_sp(net::Network, w::AbstractVector{<:Real}, solver::Solver)::Solution
    t0 = time_ns()
    m = Int(net.graph.m)
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        ww = _copy_weights(w, m)
        return _solve_sp_mip(net, ww, solver, t0; rob = false)
    end
    ww = _copy_weights(w, m)
    status, path, z = _shortest_path(net.graph, net.s, net.t, ww)
    return _finish_sp(net, status, path, z, t0; rob = false)
end
