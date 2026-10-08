function _rec_asp_finish(
    net::Network,
    c2::Vector{Float64},
    xp::_AspPair,
    t0::UInt64,
)::Solution
    g = net.graph
    m = Int(g.m)
    dt = (time_ns() - t0) / 1e9
    x = try
        path_from_seq(g, net.s, net.t, _rope_flatten(xp.xseq))
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = dt)
    end
    y = try
        path_from_seq(g, net.s, net.t, _rope_flatten(xp.yseq))
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = dt)
    end
    return _solution_rec(net, x, y, c2, t0, METHOD_COMB, 0.0)
end

function _n_arc_add!(
    tab::Matrix{_AspPath},
    i::Int,
    lc::Int,
    rc::Int,
    nL::Int,
    nR::Int,
    L::Int,
    ni::Int,
)
    cap = min(L, ni)
    tab[i, 1] = _ASP_EMPTY
    @inbounds for l in 2:cap
        best = _ASP_EMPTY
        for j in 1:(l - 1)
            j > min(L, nL) && continue
            (l - j) > min(L, nR) && continue
            cand = _asp_cat(tab[lc, j], tab[rc, l - j])
            best = _asp_better(best, cand)
        end
        tab[i, l] = best
    end
    return
end

function _solve_rec_asp(
    net::Network,
    tree::AspTree,
    c2::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    t0::UInt64,
)::Solution
    nn = Int(tree.n_nodes)
    L = k
    C = net.costs.C
    n_arcs = zeros(Int, nn)
    pair = [_ASP_PAIR_EMPTY for _ in 1:nn, _ in 1:(L + 1)]
    need_spC = nb == NB_INCLUSION
    need_sp2 = nb == NB_EXCLUSION
    need_n_arc_C = nb == NB_EXCLUSION || nb == NB_SYMDIFF
    need_n_arc_2 = nb == NB_INCLUSION || nb == NB_SYMDIFF
    spC = [_ASP_EMPTY for _ in 1:nn]
    sp2 = [_ASP_EMPTY for _ in 1:nn]
    n_arc_C = [_ASP_EMPTY for _ in 1:nn, _ in 1:max(L, 1)]
    n_arc_2 = [_ASP_EMPTY for _ in 1:nn, _ in 1:max(L, 1)]
    @inbounds for i in 1:nn
        nd = tree.nodes[i]
        if nd.op == ASP_LEAF
            n_arcs[i] = 1
            a = Int(nd.arc)
            xa = _rope_leaf(a)
            pC = _AspPath(C[a], xa)
            p2 = _AspPath(c2[a], xa)
            pair[i, 1] = _AspPair(C[a], c2[a], xa, xa)
            if need_spC
                spC[i] = pC
            end
            if need_sp2
                sp2[i] = p2
            end
            if need_n_arc_C && L >= 1
                n_arc_C[i, 1] = pC
            end
            if need_n_arc_2 && L >= 1
                n_arc_2[i, 1] = p2
            end
        else
            lc = Int(nd.left)
            rc = Int(nd.right)
            nL = n_arcs[lc]
            nR = n_arcs[rc]
            n_arcs[i] = nL + nR
            cap = min(L, n_arcs[i])
            if nd.op == ASP_SERIES
                if need_spC
                    spC[i] = _asp_cat(spC[lc], spC[rc])
                end
                if need_sp2
                    sp2[i] = _asp_cat(sp2[lc], sp2[rc])
                end
                if need_n_arc_C
                    _n_arc_add!(n_arc_C, i, lc, rc, nL, nR, L, n_arcs[i])
                end
                if need_n_arc_2
                    _n_arc_add!(n_arc_2, i, lc, rc, nL, nR, L, n_arcs[i])
                end
                for l in 0:cap
                    best = _ASP_PAIR_EMPTY
                    for j in 0:l
                        j > min(L, nL) && continue
                        (l - j) > min(L, nR) && continue
                        cand = _pair_cat(pair[lc, j + 1], pair[rc, (l - j) + 1])
                        best = _pair_better(best, cand)
                    end
                    pair[i, l + 1] = best
                end
            else
                if need_spC
                    spC[i] = _asp_better(spC[lc], spC[rc])
                end
                if need_sp2
                    sp2[i] = _asp_better(sp2[lc], sp2[rc])
                end
                if need_n_arc_C
                    for l in 1:cap
                        n_arc_C[i, l] = _asp_better(n_arc_C[lc, l], n_arc_C[rc, l])
                    end
                end
                if need_n_arc_2
                    for l in 1:cap
                        n_arc_2[i, l] = _asp_better(n_arc_2[lc, l], n_arc_2[rc, l])
                    end
                end
                pair[i, 1] = _pair_better(pair[lc, 1], pair[rc, 1])
                for l in 1:cap
                    best = _pair_better(pair[lc, l + 1], pair[rc, l + 1])
                    if nb == NB_INCLUSION
                        best = _pair_better(best, _pair_from(spC[lc], n_arc_2[rc, l]))
                        best = _pair_better(best, _pair_from(spC[rc], n_arc_2[lc, l]))
                    elseif nb == NB_EXCLUSION
                        best = _pair_better(best, _pair_from(n_arc_C[lc, l], sp2[rc]))
                        best = _pair_better(best, _pair_from(n_arc_C[rc, l], sp2[lc]))
                    else
                        for p in 1:(l - 1)
                            if p <= min(L, nL) && (l - p) <= min(L, nR)
                                best = _pair_better(best, _pair_from(n_arc_C[lc, p], n_arc_2[rc, l - p]))
                            end
                            if p <= min(L, nR) && (l - p) <= min(L, nL)
                                best = _pair_better(best, _pair_from(n_arc_C[rc, p], n_arc_2[lc, l - p]))
                            end
                        end
                    end
                    pair[i, l + 1] = best
                end
            end
        end
    end
    root = Int(tree.root)
    best = _ASP_PAIR_EMPTY
    @inbounds for l in 0:L
        best = _pair_better(best, pair[root, l + 1])
    end
    _pair_tot(best) == Inf && return solution_empty(
        Int(net.graph.m);
        status = ST_INFEASIBLE,
        method = METHOD_COMB,
        time_sec = (time_ns() - t0) / 1e9,
    )
    return _rec_asp_finish(net, c2, best, t0)
end
