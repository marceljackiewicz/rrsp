"""
    gen_random_dag(n, m; rng) -> Graph

Random DAG on vertices `1, …, n`: every arc has `tail < head`. If
`m ≥ n-1`, the path `1–2–…–n` is included so `1` reaches `n`. Extra arcs
are sampled uniformly among forward pairs, with replacement (parallel arcs
are allowed).
"""
function gen_random_dag(
    n::Integer,
    m::Integer;
    rng::AbstractRNG = Random.default_rng(),
)::Graph
    n < 1 && throw(ArgumentError("n < 1"))
    m < 0 && throw(ArgumentError("m < 0"))
    nn = Int(n)
    mm = Int(m)
    nn == 1 && mm > 0 && throw(ArgumentError("no acyclic arcs on one vertex"))
    tail = Vector{Int32}(undef, mm)
    head = Vector{Int32}(undef, mm)
    k = 0
    if nn >= 2
        npath = min(mm, nn - 1)
        for i in 1:npath
            k += 1
            tail[k] = Int32(i)
            head[k] = Int32(i + 1)
        end
    end
    while k < mm
        u = rand(rng, 1:(nn - 1))
        v = rand(rng, (u + 1):nn)
        k += 1
        tail[k] = Int32(u)
        head[k] = Int32(v)
    end
    return graph_new(nn, tail, head)
end

"""
    gen_layered_skips(H, W; p, rng) -> Graph

Wide layered digraph of [`gen_wide_layered`](@ref), plus each arc that
advances exactly two layers included independently with probability `p`.
Every arc spans one or two layers, so the shortest ``s``–``t`` path has
between ``ceil(H/2)`` and `H` arcs. `p = 0` reproduces the layered graph.
"""
function gen_layered_skips(
    H::Integer,
    W::Integer;
    p::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
)::Graph
    hh = Int(H)
    ww = Int(W)
    pp = Float64(p)
    (0.0 <= pp <= 1.0) || throw(ArgumentError("p outside [0, 1]"))
    base = gen_wide_layered(hh, ww)
    tail = Int32[base.tail;]
    head = Int32[base.head;]
    if pp > 0.0 && hh >= 2
        for i in 0:(hh - 2)
            for u in layer_vertex_ids(hh, ww, i)
                for v in layer_vertex_ids(hh, ww, i + 2)
                    if pp == 1.0 || rand(rng) < pp
                        push!(tail, Int32(u))
                        push!(head, Int32(v))
                    end
                end
            end
        end
    end
    return graph_new(Int(base.n), tail, head)
end
