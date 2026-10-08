# Constrained shortest path (integer traversal times, nonnegative costs).
# No zero-traversal-time cycle is assumed. Waiting (unused budget) is allowed.

function _sssp_fill!(
    dist::AbstractVector{Float64},
    pred::AbstractVector{Int32},
    g::Graph,
    src::Int32,
    w::Vector{Float64},
    order::Vector{Int32},
    skip::Union{Nothing,AbstractVector{UInt8}},
)
    n = Int(g.n)
    @inbounds for v in 1:n
        dist[v] = Inf
        pred[v] = Int32(0)
    end
    dist[Int(src)] = 0.0
    @inbounds for i in 1:n
        u = order[i]
        du = dist[Int(u)]
        du == Inf && continue
        for a in outgoing(g, u)
            skip !== nothing && skip[a] != 0x00 && continue
            vv = Int(g.head[a])
            nd = du + w[a]
            if nd < dist[vv]
                dist[vv] = nd
                pred[vv] = a
            end
        end
    end
    return nothing
end

# Layered shortest-path table by number of arcs: `dist[r + 1, v]` is the cheapest
# walk from `src` to `v` with exactly `r` arcs (`r = 0, …, L`), and `pred[r + 1, v]`
# is its last arc. Arcs flagged in `skip` are not used.
function _n_arc_fill!(
    dist::AbstractMatrix{Float64},
    pred::AbstractMatrix{Int32},
    g::Graph,
    src::Int32,
    w::Vector{Float64},
    L::Int,
    skip::Union{Nothing,AbstractVector{UInt8}},
)
    n = Int(g.n)
    fill!(dist, Inf)
    fill!(pred, Int32(0))
    dist[1, Int(src)] = 0.0
    L <= 0 && return nothing
    @inbounds for r in 0:(L - 1)
        col = r + 1
        ncol = r + 2
        for v in 1:n
            dv = dist[col, v]
            dv == Inf && continue
            for a in outgoing(g, Int32(v))
                skip !== nothing && skip[a] != 0x00 && continue
                vv = Int(g.head[a])
                nd = dv + w[a]
                if nd < dist[ncol, vv]
                    dist[ncol, vv] = nd
                    pred[ncol, vv] = a
                end
            end
        end
    end
    return nothing
end

# Cheapest entry of the table for `tgt` over all arc counts `0:L`: (cost, arc count).
function _n_arc_best(dist::AbstractMatrix{Float64}, tgt::Int32, L::Int)::Tuple{Float64,Int}
    best = Inf
    br = -1
    tv = Int(tgt)
    @inbounds for r in 0:L
        d = dist[r + 1, tv]
        if d < best
            best = d
            br = r
        end
    end
    return (best, br)
end

# Arcs of the walk with `br` arcs ending at `tgt`, read back through `pred`;
# empty if `br == 0` or the walk does not start at `src`.
function _n_arc_seq(
    g::Graph,
    src::Int32,
    tgt::Int32,
    pred::AbstractMatrix{Int32},
    br::Int,
)::Vector{Int32}
    br == 0 && return Int32[]
    seq = Vector{Int32}(undef, br)
    v = Int(tgt)
    @inbounds for i in br:-1:1
        a = pred[i + 1, v]
        a == Int32(0) && return Int32[]
        seq[i] = a
        v = Int(g.tail[a])
    end
    v == Int(src) || return Int32[]
    return seq
end

function _simplify_arc_walk(g::Graph, s::Int32, t::Int32, seq::Vector{Int32})::Vector{Int32}
    isempty(seq) && return seq
    n = Int(g.n)
    onstack = zeros(UInt8, n)
    stack_v = Vector{Int32}(undef, n + 1)
    stack_a = Vector{Int32}(undef, length(seq))
    sv = 1
    sa = 0
    stack_v[1] = s
    onstack[Int(s)] = 0x01
    @inbounds for k in 1:length(seq)
        a = seq[k]
        w = g.head[a]
        if onstack[Int(w)] != 0x00
            while stack_v[sv] != w
                onstack[Int(stack_v[sv])] = 0x00
                sv -= 1
                sa -= 1
                sv < 1 && return seq
            end
            continue
        end
        sa += 1
        stack_a[sa] = a
        sv += 1
        stack_v[sv] = w
        onstack[Int(w)] = 0x01
    end
    stack_v[sv] == t || return seq
    out = Vector{Int32}(undef, sa)
    @inbounds for i in 1:sa
        out[i] = stack_a[i]
    end
    return out
end

"""
    _csp(n, s, t, tails, heads, times, costs, R) -> (status, seq, z)

Minimum-cost `s`–`t` walk of total traversal time at most `R`. `seq` is a list of
arc indices into the parallel arrays `tails`, `heads`, `times`, `costs`.
"""
function _csp(
    n::Int,
    s::Int32,
    t::Int32,
    tails::Vector{Int32},
    heads::Vector{Int32},
    times::Vector{Int32},
    costs::Vector{Float64},
    R::Int,
)::Tuple{Status,Vector{Int32},Float64}
    md = length(tails)
    if s == t
        return (ST_OK, Int32[], 0.0)
    end
    R < 0 && return (ST_INFEASIBLE, Int32[], Inf)
    zero_as = Int32[]
    pos_as = Int32[]
    @inbounds for a in 1:md
        if times[a] == Int32(0)
            push!(zero_as, Int32(a))
        else
            push!(pos_as, Int32(a))
        end
    end
    nr = R + 1
    dist = fill(Inf, n, nr)
    pred_a = zeros(Int32, n, nr)
    pred_layer = fill(Int32(-1), n, nr)
    dist[Int(s), 1] = 0.0
    @inbounds for r in 0:R
        col = r + 1
        for _ in 1:n
            changed = false
            for a in zero_as
                u = Int(tails[a])
                v = Int(heads[a])
                du = dist[u, col]
                du == Inf && continue
                nd = du + costs[a]
                if nd < dist[v, col]
                    dist[v, col] = nd
                    pred_a[v, col] = a
                    pred_layer[v, col] = Int32(r)
                    changed = true
                end
            end
            changed || break
        end
        if r < R
            ncol = col + 1
            for v in 1:n
                if dist[v, col] < dist[v, ncol]
                    dist[v, ncol] = dist[v, col]
                    pred_a[v, ncol] = Int32(0)
                    pred_layer[v, ncol] = Int32(r)
                end
            end
        end
        for a in pos_as
            tau = Int(times[a])
            rr = r + tau
            rr > R && continue
            u = Int(tails[a])
            v = Int(heads[a])
            du = dist[u, col]
            du == Inf && continue
            nd = du + costs[a]
            c2 = rr + 1
            if nd < dist[v, c2]
                dist[v, c2] = nd
                pred_a[v, c2] = a
                pred_layer[v, c2] = Int32(r)
            end
        end
    end
    tt = Int(t)
    best = Inf
    br = -1
    @inbounds for r in 0:R
        d = dist[tt, r + 1]
        if d < best
            best = d
            br = r
        end
    end
    best == Inf && return (ST_INFEASIBLE, Int32[], Inf)
    seq_rev = Int32[]
    v = tt
    r = br
    steps = 0
    max_steps = n * nr + 8
    while true
        steps += 1
        steps > max_steps && return (ST_ERROR, Int32[], Inf)
        col = r + 1
        a = pred_a[v, col]
        pr = pred_layer[v, col]
        if pr < 0 && a == Int32(0)
            break
        end
        if a == Int32(0)
            r = Int(pr)
            continue
        end
        push!(seq_rev, a)
        v = Int(tails[a])
        r = Int(pr)
    end
    v == Int(s) || return (ST_ERROR, Int32[], Inf)
    reverse!(seq_rev)
    return (ST_OK, seq_rev, best)
end

function _concat_seqs(dseq::Vector{Int32}, parts::Vector{Vector{Int32}})::Vector{Int32}
    ntot = 0
    @inbounds for a in dseq
        ntot += length(parts[a])
    end
    out = Vector{Int32}(undef, ntot)
    k = 0
    @inbounds for a in dseq
        part = parts[a]
        for i in 1:length(part)
            k += 1
            out[k] = part[i]
        end
    end
    return out
end

# A simple s–t path has at most n-1 arcs, so every neighborhood (inclusion,
# exclusion, symmetric difference) is already saturated at k = 2(n-1). Larger
# values give the same problem; clamping bounds the size of the time-expanded
# networks.
_clamp_k(g::Graph, k::Integer)::Int = min(Int(k), 2 * max(Int(g.n) - 1, 1))
