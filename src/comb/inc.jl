function _inc_vertices(g::Graph, s::Int32, x::Path)::Vector{Int32}
    k = length(x.seq)
    vx = Vector{Int32}(undef, k + 1)
    vx[1] = s
    @inbounds for i in 1:k
        vx[i + 1] = g.head[x.seq[i]]
    end
    return vx
end

function _inc_from_dyn(
    net::Network,
    x::Path,
    w::Vector{Float64},
    dtail::Vector{Int32},
    dhead::Vector{Int32},
    dtime::Vector{Int32},
    dcost::Vector{Float64},
    dyseq::Vector{Vector{Int32}},
    R::Int,
    t0::UInt64,
)::Solution
    g = net.graph
    m = Int(g.m)
    st, dseq, _z = _csp(Int(g.n), net.s, net.t, dtail, dhead, dtime, dcost, R)
    dt = (time_ns() - t0) / 1e9
    if st != ST_OK
        return solution_empty(m; status = st, method = METHOD_COMB, time_sec = dt)
    end
    yseq = _concat_seqs(dseq, dyseq)
    yseq = _simplify_arc_walk(g, net.s, net.t, yseq)
    y = try
        path_from_seq(g, net.s, net.t, yseq)
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = (time_ns() - t0) / 1e9)
    end
    return _solution_inc(net, x, y, w, t0, METHOD_COMB, 0.0)
end

function _inc_inclusion(net::Network, x::Path, w::Vector{Float64}, k::Int, t0::UInt64)::Solution
    g = net.graph
    m = Int(g.m)
    times = Vector{Int32}(undef, m)
    @inbounds for a in 1:m
        times[a] = x.chi[a] == 0x00 ? Int32(1) : Int32(0)
    end
    R = min(k, m)
    st, seq, _z = _csp(Int(g.n), net.s, net.t, g.tail, g.head, times, w, R)
    dt = (time_ns() - t0) / 1e9
    if st != ST_OK
        return solution_empty(m; status = st, method = METHOD_COMB, time_sec = dt)
    end
    seq = _simplify_arc_walk(g, net.s, net.t, seq)
    y = try
        path_from_seq(g, net.s, net.t, seq)
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = (time_ns() - t0) / 1e9)
    end
    return _solution_inc(net, x, y, w, t0, METHOD_COMB, 0.0)
end

function _inc_excl_dag(net::Network, x::Path, w::Vector{Float64}, k::Int, t0::UInt64)::Solution
    g = net.graph
    n = Int(g.n)
    order = _topo_order(g)
    order === nothing && return _not_impl(n > 0 ? Int(g.m) : 0, METHOD_COMB, t0)
    vx = _inc_vertices(g, net.s, x)
    nxs = length(x.seq)
    dtail = Int32[]
    dhead = Int32[]
    dtime = Int32[]
    dcost = Float64[]
    dyseq = Vector{Vector{Int32}}()
    @inbounds for i in 1:nxs
        a = x.seq[i]
        push!(dtail, g.tail[a])
        push!(dhead, g.head[a])
        push!(dtime, Int32(0))
        push!(dcost, w[a])
        push!(dyseq, Int32[a])
    end
    skip = x.chi
    dist = Vector{Float64}(undef, n)
    pred = Vector{Int32}(undef, n)
    nvx = nxs + 1
    @inbounds for i in 1:nvx
        _sssp_fill!(dist, pred, g, vx[i], w, order, skip)
        for j in (i + 1):nvx
            p = j - i
            p <= k || continue
            tgt = vx[j]
            st, path, zsp = _path_from_pred(g, vx[i], tgt, pred, dist)
            st == ST_OK || continue
            push!(dtail, vx[i])
            push!(dhead, tgt)
            push!(dtime, Int32(p))
            push!(dcost, zsp)
            push!(dyseq, copy(path.seq))
        end
    end
    return _inc_from_dyn(net, x, w, dtail, dhead, dtime, dcost, dyseq, k, t0)
end

function _inc_sym_dag(net::Network, x::Path, w::Vector{Float64}, k::Int, t0::UInt64)::Solution
    g = net.graph
    n = Int(g.n)
    m = Int(g.m)
    order = _topo_order(g)
    order === nothing && return _not_impl(m, METHOD_COMB, t0)
    vx = _inc_vertices(g, net.s, x)
    nxs = length(x.seq)
    dtail = Int32[]
    dhead = Int32[]
    dtime = Int32[]
    dcost = Float64[]
    dyseq = Vector{Vector{Int32}}()
    @inbounds for i in 1:nxs
        a = x.seq[i]
        push!(dtail, g.tail[a])
        push!(dhead, g.head[a])
        push!(dtime, Int32(0))
        push!(dcost, w[a])
        push!(dyseq, Int32[a])
    end
    L = min(k, m)
    skip = x.chi
    n_arc_dist = fill(Inf, L + 1, n)
    n_arc_pred = zeros(Int32, L + 1, n)
    nvx = nxs + 1
    @inbounds for i in 1:nvx
        _n_arc_fill!(n_arc_dist, n_arc_pred, g, vx[i], w, L, skip)
        for j in (i + 1):nvx
            p = j - i
            p < k || continue
            tgt = vx[j]
            for l in 1:(k - p)
                best, br = _n_arc_best(n_arc_dist, tgt, l)
                best == Inf && continue
                br <= 0 && continue
                yseq = _n_arc_seq(g, vx[i], tgt, n_arc_pred, br)
                isempty(yseq) && continue
                push!(dtail, vx[i])
                push!(dhead, tgt)
                push!(dtime, Int32(p + l))
                push!(dcost, best)
                push!(dyseq, yseq)
            end
        end
    end
    return _inc_from_dyn(net, x, w, dtail, dhead, dtime, dcost, dyseq, k, t0)
end

function _solve_inc_comb(
    net::Network,
    x::Path,
    w::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    t0::UInt64,
)::Solution
    if net.s == net.t || k == 0
        return _solution_inc(net, x, x, w, t0, METHOD_COMB, 0.0)
    end
    k = _clamp_k(net.graph, k)
    tree = _asp_tree_cached(net.graph, net.s, net.t)
    if tree !== nothing
        return _solve_inc_asp(net, tree, x, w, nb, k, t0)
    end
    if nb == NB_INCLUSION
        return _inc_inclusion(net, x, w, k, t0)
    elseif nb == NB_EXCLUSION
        return _inc_excl_dag(net, x, w, k, t0)
    else
        return _inc_sym_dag(net, x, w, k, t0)
    end
end
