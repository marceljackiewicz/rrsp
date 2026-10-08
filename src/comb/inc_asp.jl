function _inc_asp_finish(net::Network, x::Path, w::Vector{Float64}, yseq::Vector{Int32}, t0::UInt64)::Solution
    g = net.graph
    m = Int(g.m)
    y = try
        path_from_seq(g, net.s, net.t, yseq)
    catch
        return solution_empty(m; status = ST_ERROR, method = METHOD_COMB, time_sec = (time_ns() - t0) / 1e9)
    end
    return _solution_inc(net, x, y, w, t0, METHOD_COMB, 0.0)
end

function _inc_asp_best(slot::Matrix{_AspPath}, root::Int, k::Int)::_AspPath
    best = _ASP_EMPTY
    @inbounds for l in 0:k
        best = _asp_better(best, slot[root, l + 1])
    end
    return best
end

function _inc_asp_incl(
    net::Network,
    tree::AspTree,
    x::Path,
    w::Vector{Float64},
    k::Int,
    t0::UInt64,
)::Solution
    nn = Int(tree.n_nodes)
    L = k
    n_arcs = zeros(Int, nn)
    slot = [_ASP_EMPTY for _ in 1:nn, _ in 1:(L + 1)]
    @inbounds for i in 1:nn
        nd = tree.nodes[i]
        if nd.op == ASP_LEAF
            n_arcs[i] = 1
            a = Int(nd.arc)
            p = _AspPath(w[a], _rope_leaf(a))
            if x.chi[a] != 0x00
                slot[i, 1] = p
            elseif L >= 1
                slot[i, 2] = p
            end
        else
            lc = Int(nd.left)
            rc = Int(nd.right)
            nL = n_arcs[lc]
            nR = n_arcs[rc]
            n_arcs[i] = nL + nR
            cap = min(L, n_arcs[i])
            if nd.op == ASP_SERIES
                for l in 0:cap
                    best = _ASP_EMPTY
                    for j in 0:l
                        j > min(L, nL) && continue
                        (l - j) > min(L, nR) && continue
                        cand = _asp_cat(slot[lc, j + 1], slot[rc, (l - j) + 1])
                        best = _asp_better(best, cand)
                    end
                    slot[i, l + 1] = best
                end
            else
                for l in 0:cap
                    slot[i, l + 1] = _asp_better(slot[lc, l + 1], slot[rc, l + 1])
                end
            end
        end
    end
    best = _inc_asp_best(slot, Int(tree.root), L)
    best.cost == Inf && return solution_empty(Int(net.graph.m); status = ST_INFEASIBLE, method = METHOD_COMB, time_sec = (time_ns() - t0) / 1e9)
    return _inc_asp_finish(net, x, w, _rope_flatten(best.seq), t0)
end

function _inc_asp_excl(
    net::Network,
    tree::AspTree,
    x::Path,
    w::Vector{Float64},
    k::Int,
    t0::UInt64,
)::Solution
    nn = Int(tree.n_nodes)
    L = k
    n_arcs = zeros(Int, nn)
    nx = zeros(Int, nn)
    slot = [_ASP_EMPTY for _ in 1:nn, _ in 1:(L + 1)]
    avoid = [_ASP_EMPTY for _ in 1:nn]
    @inbounds for i in 1:nn
        nd = tree.nodes[i]
        if nd.op == ASP_LEAF
            n_arcs[i] = 1
            a = Int(nd.arc)
            p = _AspPath(w[a], _rope_leaf(a))
            if x.chi[a] != 0x00
                nx[i] = 1
                slot[i, 1] = p
            else
                avoid[i] = p
            end
        else
            lc = Int(nd.left)
            rc = Int(nd.right)
            nL = n_arcs[lc]
            nR = n_arcs[rc]
            n_arcs[i] = nL + nR
            nx[i] = nx[lc] + nx[rc]
            cap = min(L, n_arcs[i])
            if nd.op == ASP_SERIES
                avoid[i] = _asp_cat(avoid[lc], avoid[rc])
                if nx[lc] > 0 && nx[rc] > 0
                    for l in 0:cap
                        best = _ASP_EMPTY
                        for j in 0:l
                            j > min(L, nL) && continue
                            (l - j) > min(L, nR) && continue
                            cand = _asp_cat(slot[lc, j + 1], slot[rc, (l - j) + 1])
                            best = _asp_better(best, cand)
                        end
                        slot[i, l + 1] = best
                    end
                elseif nx[lc] > 0
                    for l in 0:cap
                        slot[i, l + 1] = _asp_cat(slot[lc, l + 1], avoid[rc])
                    end
                elseif nx[rc] > 0
                    for l in 0:cap
                        slot[i, l + 1] = _asp_cat(avoid[lc], slot[rc, l + 1])
                    end
                end
            else
                avoid[i] = _asp_better(avoid[lc], avoid[rc])
                if nx[lc] > 0 && nx[rc] == 0
                    for l in 0:cap
                        slot[i, l + 1] = slot[lc, l + 1]
                    end
                    kk = nx[lc]
                    if kk <= L
                        slot[i, kk + 1] = _asp_better(slot[i, kk + 1], avoid[rc])
                    end
                elseif nx[rc] > 0 && nx[lc] == 0
                    for l in 0:cap
                        slot[i, l + 1] = slot[rc, l + 1]
                    end
                    kk = nx[rc]
                    if kk <= L
                        slot[i, kk + 1] = _asp_better(slot[i, kk + 1], avoid[lc])
                    end
                elseif nx[lc] > 0 && nx[rc] > 0
                    for l in 0:cap
                        slot[i, l + 1] = _asp_better(slot[lc, l + 1], slot[rc, l + 1])
                    end
                end
            end
        end
    end
    best = _inc_asp_best(slot, Int(tree.root), L)
    best.cost == Inf && return solution_empty(Int(net.graph.m); status = ST_INFEASIBLE, method = METHOD_COMB, time_sec = (time_ns() - t0) / 1e9)
    return _inc_asp_finish(net, x, w, _rope_flatten(best.seq), t0)
end

function _inc_asp_sym(
    net::Network,
    tree::AspTree,
    x::Path,
    w::Vector{Float64},
    k::Int,
    t0::UInt64,
)::Solution
    nn = Int(tree.n_nodes)
    L = k
    n_arcs = zeros(Int, nn)
    nx = zeros(Int, nn)
    slot = [_ASP_EMPTY for _ in 1:nn, _ in 1:(L + 1)]
    @inbounds for i in 1:nn
        nd = tree.nodes[i]
        if nd.op == ASP_LEAF
            n_arcs[i] = 1
            a = Int(nd.arc)
            p = _AspPath(w[a], _rope_leaf(a))
            if x.chi[a] != 0x00
                nx[i] = 1
                slot[i, 1] = p
            elseif L >= 1
                slot[i, 2] = p
            end
        else
            lc = Int(nd.left)
            rc = Int(nd.right)
            nL = n_arcs[lc]
            nR = n_arcs[rc]
            n_arcs[i] = nL + nR
            nx[i] = nx[lc] + nx[rc]
            cap = min(L, n_arcs[i])
            if nd.op == ASP_SERIES
                for l in 0:cap
                    best = _ASP_EMPTY
                    for j in 0:l
                        j > min(L, nL) && continue
                        (l - j) > min(L, nR) && continue
                        cand = _asp_cat(slot[lc, j + 1], slot[rc, (l - j) + 1])
                        best = _asp_better(best, cand)
                    end
                    slot[i, l + 1] = best
                end
            else
                for l in 0:cap
                    best = _ASP_EMPTY
                    kk = l - nx[rc]
                    if kk >= 0 && kk <= L
                        best = _asp_better(best, slot[lc, kk + 1])
                    end
                    kk = l - nx[lc]
                    if kk >= 0 && kk <= L
                        best = _asp_better(best, slot[rc, kk + 1])
                    end
                    slot[i, l + 1] = best
                end
            end
        end
    end
    best = _inc_asp_best(slot, Int(tree.root), L)
    best.cost == Inf && return solution_empty(Int(net.graph.m); status = ST_INFEASIBLE, method = METHOD_COMB, time_sec = (time_ns() - t0) / 1e9)
    return _inc_asp_finish(net, x, w, _rope_flatten(best.seq), t0)
end

function _solve_inc_asp(
    net::Network,
    tree::AspTree,
    x::Path,
    w::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    t0::UInt64,
)::Solution
    if nb == NB_INCLUSION
        return _inc_asp_incl(net, tree, x, w, k, t0)
    elseif nb == NB_EXCLUSION
        return _inc_asp_excl(net, tree, x, w, k, t0)
    else
        return _inc_asp_sym(net, tree, x, w, k, t0)
    end
end
