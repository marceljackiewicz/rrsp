function _arc(g::Graph, a::Integer)::Int32
    ia = Int(a)
    (1 <= ia <= Int(g.m)) || throw(ArgumentError("arc id $a is outside 1:$(g.m)"))
    return Int32(ia)
end

"""
    path_empty(m) -> Path

A path with no arcs on a graph of `m` arcs: zero characteristic vector and
empty sequence.
"""
function path_empty(m::Integer)::Path
    m < 0 && throw(ArgumentError("m < 0"))
    return Path(zeros(UInt8, Int(m)), Int32[])
end

"""
    path_from_seq(g, s, t, seq) -> Path

Construct a simple ``s``–``t`` path from a sequence of arc indices.

The empty sequence is admissible if and only if `s == t` (the trivial path).
The sequence is copied.

# Throws
- `ArgumentError`: if the walk is not a simple ``s``–``t`` path in `g`.
"""
function path_from_seq(g::Graph, s::Integer, t::Integer, seq::AbstractVector{<:Integer})::Path
    ss = _vertex(g, s)
    tt = _vertex(g, t)
    m = Int(g.m)
    if isempty(seq)
        ss == tt || throw(ArgumentError("empty sequence does not join $ss to $tt"))
        return Path(zeros(UInt8, m), Int32[])
    end
    chi = zeros(UInt8, m)
    visited = falses(Int(g.n))
    visited[ss] = true
    v = ss
    seq32 = Vector{Int32}(undef, length(seq))
    @inbounds for k in 1:length(seq)
        a = _arc(g, seq[k])
        chi[a] == 0x00 || throw(ArgumentError("repeated arc $a; path not simple"))
        g.tail[a] == v || throw(ArgumentError("arc $a does not leave $v"))
        w = g.head[a]
        visited[w] && throw(ArgumentError("walk is not a simple path"))
        visited[w] = true
        chi[a] = 0x01
        seq32[k] = a
        v = w
    end
    v == tt || throw(ArgumentError("walk does not end at t"))
    return Path(chi, seq32)
end

"""
    path_from_chi(g, s, t, chi) -> Path

Reconstruct the unique simple ``s``–``t`` path whose characteristic vector is
`chi`, or throw if `chi` does not encode such a path.

The characteristic vector is copied.
"""
function path_from_chi(g::Graph, s::Integer, t::Integer, chi::AbstractVector{<:Integer})::Path
    ss = _vertex(g, s)
    tt = _vertex(g, t)
    m = Int(g.m)
    length(chi) == m || throw(ArgumentError("characteristic vector has length $(length(chi)), expected $m"))
    chi_u8 = Vector{UInt8}(undef, m)
    nsel = 0
    @inbounds for a in 1:m
        c = Int(chi[a])
        (c == 0 || c == 1) || throw(ArgumentError("characteristic vector entries must be 0 or 1"))
        chi_u8[a] = UInt8(c)
        nsel += c
    end
    if ss == tt
        nsel == 0 || throw(ArgumentError("a nonempty closed walk is not a simple s–t path"))
        return Path(chi_u8, Int32[])
    end
    seq = Vector{Int32}(undef, nsel)
    len = 0
    visited = falses(Int(g.n))
    visited[ss] = true
    v = ss
    used = 0
    while v != tt
        nxt = Int32(0)
        count = 0
        @inbounds for a in outgoing(g, v)
            if chi_u8[a] != 0x00
                count += 1
                nxt = a
            end
        end
        count == 1 || throw(ArgumentError("characteristic vector is not an s–t path"))
        w = g.head[nxt]
        visited[w] && throw(ArgumentError("characteristic vector is not a simple path"))
        visited[w] = true
        len += 1
        seq[len] = nxt
        used += 1
        v = w
        used > nsel && throw(ArgumentError("characteristic vector is not an s–t path"))
    end
    used == nsel || throw(ArgumentError("characteristic vector is not an s–t path"))
    return Path(chi_u8, seq)
end

"""
    path_cost(p, w) -> Float64

Return the inner product of `w` with the characteristic vector of `p`.

# Throws
- `ArgumentError`: if `w` and `p.chi` differ in length.
"""
function path_cost(p::Path, w::AbstractVector{<:Real})::Float64
    length(w) == length(p.chi) || throw(ArgumentError("weight length mismatch"))
    z = 0.0
    @inbounds for a in 1:length(p.chi)
        if p.chi[a] != 0x00
            z += Float64(w[a])
        end
    end
    return z
end

"""
    path_length(p) -> Int32

Number of arcs in `p.seq`.
"""
function path_length(p::Path)::Int32
    return Int32(length(p.seq))
end
