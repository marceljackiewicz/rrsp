"""
    gen_repsel(widths) -> Graph

Representatives-selection ASP graph: a series of parallel bundles. Vertex
`i` is joined to `i+1` by `widths[i]` parallel arcs. Terminals are `1` and
`length(widths)+1`. With equal widths this is a chain of beads, as in the
`two_beads` example.
"""
function gen_repsel(widths::AbstractVector{<:Integer})::Graph
    isempty(widths) && throw(ArgumentError("widths is empty"))
    tail = Int32[]
    head = Int32[]
    for (i, w) in enumerate(widths)
        w >= 1 || throw(ArgumentError("bundle width < 1"))
        u = Int32(i)
        v = Int32(i + 1)
        for _ in 1:Int(w)
            push!(tail, u)
            push!(head, v)
        end
    end
    return graph_new(length(widths) + 1, tail, head)
end

"""
    gen_asp(n_arcs; rng) -> Graph

Random two-terminal arc-series-parallel graph with exactly `n_arcs` arcs,
built by a sequence of series subdivisions and parallel duplications from a
single `1 → 2` arc. Terminals remain `1` and `2`.
"""
function gen_asp(n_arcs::Integer; rng::AbstractRNG = Random.default_rng())::Graph
    n_arcs >= 1 || throw(ArgumentError("n_arcs < 1"))
    tail = Int32[1]
    head = Int32[2]
    n = 2
    while length(tail) < Int(n_arcs)
        i = rand(rng, 1:length(tail))
        if length(tail) == 1 || rand(rng, Bool)
            n += 1
            v = head[i]
            head[i] = Int32(n)
            push!(tail, Int32(n))
            push!(head, v)
        else
            push!(tail, tail[i])
            push!(head, head[i])
        end
    end
    return graph_new(n, tail, head)
end

"""
    gen_asp_detours(H) -> Graph

Series of `H` gadgets. Each gadget is a direct arc in parallel with a
path of three arcs. A shortest ``s``–``t`` path has `H` arcs and a longest
has `3H` arcs. Terminals are `1` and `3H+1`.
"""
function gen_asp_detours(H::Integer)::Graph
    hh = Int(H)
    hh >= 1 || throw(ArgumentError("H < 1"))
    n = 1 + 3 * hh
    tail = Int32[]
    head = Int32[]
    for g in 0:(hh - 1)
        u = Int32(1 + 3 * g)
        p = Int32(u + 1)
        q = Int32(u + 2)
        v = Int32(u + 3)
        push!(tail, u)
        push!(head, v)
        push!(tail, u)
        push!(head, p)
        push!(tail, p)
        push!(head, q)
        push!(tail, q)
        push!(head, v)
    end
    return graph_new(n, tail, head)
end
