function _push_dyn!(
    dtail::Vector{Int32},
    dhead::Vector{Int32},
    dtime::Vector{Int32},
    dcost::Vector{Float64},
    dxseq::Vector{Vector{Int32}},
    dyseq::Vector{Vector{Int32}},
    u::Int32,
    v::Int32,
    tau::Int,
    c::Float64,
    xs::Vector{Int32},
    ys::Vector{Int32},
)
    push!(dtail, u)
    push!(dhead, v)
    push!(dtime, Int32(tau))
    push!(dcost, c)
    push!(dxseq, xs)
    push!(dyseq, ys)
    return nothing
end

function _solve_rec_dag(
    net::Network,
    c2::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    t0::UInt64,
)::Solution
    g = net.graph
    n = Int(g.n)
    m = Int(g.m)
    C = net.costs.C
    order = _topo_order(g)
    if order === nothing
        return _not_impl(m, METHOD_COMB, t0)
    end
    idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        idx[Int(order[i])] = Int32(i)
    end
    dtail = Int32[]
    dhead = Int32[]
    dtime = Int32[]
    dcost = Float64[]
    dxseq = Vector{Vector{Int32}}()
    dyseq = Vector{Vector{Int32}}()
    @inbounds for a in 1:m
        _push_dyn!(
            dtail, dhead, dtime, dcost, dxseq, dyseq,
            g.tail[a], g.head[a], 0, C[a] + c2[a], Int32[Int32(a)], Int32[Int32(a)],
        )
    end
    if k > 0
        L = min(k, max(m, 1))
        need_spC = nb == NB_INCLUSION
        need_sp2 = nb == NB_EXCLUSION
        need_n_arc_C = nb == NB_EXCLUSION || nb == NB_SYMDIFF
        need_n_arc_2 = nb == NB_INCLUSION || nb == NB_SYMDIFF
        # Allocate only the tables the neighborhood reads.
        nsc = need_spC ? n : 0
        ns2 = need_sp2 ? n : 0
        spC_dist = fill(Inf, nsc, nsc)
        spC_pred = zeros(Int32, nsc, nsc)
        sp2_dist = fill(Inf, ns2, ns2)
        sp2_pred = zeros(Int32, ns2, ns2)
        if need_spC || need_sp2
            dbuf = Vector{Float64}(undef, n)
            pbuf = Vector{Int32}(undef, n)
            @inbounds for src in Int32(1):Int32(n)
                if need_spC
                    _sssp_fill!(dbuf, pbuf, g, src, C, order, nothing)
                    for v in 1:n
                        spC_dist[Int(src), v] = dbuf[v]
                        spC_pred[Int(src), v] = pbuf[v]
                    end
                end
                if need_sp2
                    _sssp_fill!(dbuf, pbuf, g, src, c2, order, nothing)
                    for v in 1:n
                        sp2_dist[Int(src), v] = dbuf[v]
                        sp2_pred[Int(src), v] = pbuf[v]
                    end
                end
            end
        end
        nac = need_n_arc_C ? n : 0
        na2 = need_n_arc_2 ? n : 0
        n_arc_C_dist = fill(Inf, nac, need_n_arc_C ? L + 1 : 0, nac)
        n_arc_C_pred = zeros(Int32, nac, need_n_arc_C ? L + 1 : 0, nac)
        n_arc_2_dist = fill(Inf, na2, need_n_arc_2 ? L + 1 : 0, na2)
        n_arc_2_pred = zeros(Int32, na2, need_n_arc_2 ? L + 1 : 0, na2)
        if need_n_arc_C || need_n_arc_2
            @inbounds for src in Int32(1):Int32(n)
                si = Int(src)
                if need_n_arc_C
                    _n_arc_fill!(view(n_arc_C_dist, si, :, :), view(n_arc_C_pred, si, :, :), g, src, C, L, nothing)
                end
                if need_n_arc_2
                    _n_arc_fill!(view(n_arc_2_dist, si, :, :), view(n_arc_2_pred, si, :, :), g, src, c2, L, nothing)
                end
            end
        end
        @inbounds for u in Int32(1):Int32(n)
            ui = Int(u)
            for v in Int32(1):Int32(n)
                idx[ui] < idx[Int(v)] || continue
                if nb == NB_INCLUSION
                    st, xpath, zx = _path_from_pred(g, u, v, view(spC_pred, ui, :), view(spC_dist, ui, :))
                    st == ST_OK || continue
                    xs = copy(xpath.seq)
                    n_arc_dist = view(n_arc_2_dist, ui, :, :)
                    n_arc_pred = view(n_arc_2_pred, ui, :, :)
                    for l in 1:k
                        ll = min(l, L)
                        best, br = _n_arc_best(n_arc_dist, v, ll)
                        best == Inf && continue
                        br <= 0 && continue
                        ys = _n_arc_seq(g, u, v, n_arc_pred, br)
                        isempty(ys) && continue
                        _push_dyn!(dtail, dhead, dtime, dcost, dxseq, dyseq, u, v, l, zx + best, xs, ys)
                    end
                elseif nb == NB_EXCLUSION
                    st, ypath, zy = _path_from_pred(g, u, v, view(sp2_pred, ui, :), view(sp2_dist, ui, :))
                    st == ST_OK || continue
                    ys = copy(ypath.seq)
                    n_arc_dist = view(n_arc_C_dist, ui, :, :)
                    n_arc_pred = view(n_arc_C_pred, ui, :, :)
                    for l in 1:k
                        ll = min(l, L)
                        best, br = _n_arc_best(n_arc_dist, v, ll)
                        best == Inf && continue
                        br <= 0 && continue
                        xs = _n_arc_seq(g, u, v, n_arc_pred, br)
                        isempty(xs) && continue
                        _push_dyn!(dtail, dhead, dtime, dcost, dxseq, dyseq, u, v, l, best + zy, xs, ys)
                    end
                else
                    dist_C = view(n_arc_C_dist, ui, :, :)
                    pred_C = view(n_arc_C_pred, ui, :, :)
                    dist_2 = view(n_arc_2_dist, ui, :, :)
                    pred_2 = view(n_arc_2_pred, ui, :, :)
                    for l in 2:k
                        best = Inf
                        bxs = Int32[]
                        bys = Int32[]
                        for p in 1:(l - 1)
                            pc = min(p, L)
                            py = min(l - p, L)
                            cx, rx = _n_arc_best(dist_C, v, pc)
                            cy, ry = _n_arc_best(dist_2, v, py)
                            (cx == Inf || cy == Inf) && continue
                            tot = cx + cy
                            if tot < best
                                xs = _n_arc_seq(g, u, v, pred_C, rx)
                                ys = _n_arc_seq(g, u, v, pred_2, ry)
                                (isempty(xs) || isempty(ys)) && continue
                                best = tot
                                bxs = xs
                                bys = ys
                            end
                        end
                        best == Inf && continue
                        _push_dyn!(dtail, dhead, dtime, dcost, dxseq, dyseq, u, v, l, best, bxs, bys)
                    end
                end
            end
        end
    end
    st, dseq, _z = _csp(n, net.s, net.t, dtail, dhead, dtime, dcost, k)
    dt = (time_ns() - t0) / 1e9
    if st != ST_OK
        return solution_empty(m; status = st, method = METHOD_COMB, time_sec = dt)
    end
    xseq = _concat_seqs(dseq, dxseq)
    yseq = _concat_seqs(dseq, dyseq)
    xseq = _simplify_arc_walk(g, net.s, net.t, xseq)
    yseq = _simplify_arc_walk(g, net.s, net.t, yseq)
    xp = try
        path_from_seq(g, net.s, net.t, xseq)
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = dt)
    end
    yp = try
        path_from_seq(g, net.s, net.t, yseq)
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = dt)
    end
    return _solution_rec(net, xp, yp, c2, t0, METHOD_COMB, 0.0)
end
