"""
    enumerate_st_paths(g, s, t; limit=1_000_000) -> Vector{Path}

All simple ``s``–``t`` paths of `g`, in depth-first order (arcs of each vertex
in increasing index). For ``s = t`` the result is the single empty path.

The number of paths is exponential in general, so the function refuses to
return more than `limit` of them.

# Throws
- `ArgumentError`: if a vertex is out of range or there are more than `limit`
  paths.
"""
function enumerate_st_paths(g::Graph, s::Integer, t::Integer; limit::Integer = 1_000_000)::Vector{Path}
    ss = _vertex(g, s)
    tt = _vertex(g, t)
    m = Int(g.m)
    found = Path[]
    if ss == tt
        push!(found, path_empty(m))
        return found
    end
    seq = Int32[]
    on_path = falses(Int(g.n))
    on_path[ss] = true
    function dfs(v::Int32)
        if v == tt
            chi = zeros(UInt8, m)
            for a in seq
                chi[a] = 0x01
            end
            length(found) >= limit && throw(ArgumentError("more than $limit simple s–t paths"))
            push!(found, Path(chi, copy(seq)))
            return
        end
        for a in outgoing(g, v)
            w = g.head[a]
            on_path[w] && continue
            on_path[w] = true
            push!(seq, a)
            dfs(w)
            pop!(seq)
            on_path[w] = false
        end
        return
    end
    dfs(ss)
    return found
end
