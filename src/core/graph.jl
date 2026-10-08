"""
    graph_new(n, tail, head; label=nothing) -> Graph

Construct a multidigraph on `n` vertices from endpoint arrays.

`tail` and `head` have length `m`. Arc `a` is the directed pair
`(tail[a], head[a])`, with both endpoints in `1, …, n`. The forward and reverse
stars are built so that, for every vertex `v`, the block
`out_arc[first_out[v]:first_out[v+1]-1]` contains the indices of all arcs with
tail `v`, in increasing order, and likewise for `first_in` and `in_arc`.

If `label` is omitted, `label[v] = v`. The arrays passed by the caller are
copied; the resulting [`Graph`](@ref) does not alias them.

# Arguments
- `n`: number of vertices.
- `tail`: tails of the arcs.
- `head`: heads of the arcs.
- `label`: optional external vertex identifiers, of length `n`.

# Throws
- `ArgumentError`: if `n < 1`, if `tail` and `head` differ in length, if an
  endpoint lies outside `1, …, n`, if `label` has the wrong length or repeats
  a value, or if `m` exceeds the `Int32` range.
"""
function graph_new(
    n::Integer,
    tail::AbstractVector{<:Integer},
    head::AbstractVector{<:Integer};
    label::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)::Graph
    n < 1 && throw(ArgumentError("n < 1"))
    n > typemax(Int32) && throw(ArgumentError("n exceeds Int32"))
    length(tail) == length(head) || throw(ArgumentError("tail and head have different lengths"))
    length(tail) > typemax(Int32) && throw(ArgumentError("m exceeds Int32"))
    n32 = Int32(n)
    m32 = Int32(length(tail))
    if label !== nothing && length(label) != Int(n32)
        throw(ArgumentError("label has length $(length(label)), expected $n32"))
    end
    if label !== nothing && length(unique(label)) != length(label)
        throw(ArgumentError("vertex labels are not unique"))
    end

    tail32 = Vector{Int32}(undef, Int(m32))
    head32 = Vector{Int32}(undef, Int(m32))
    out_deg = zeros(Int32, Int(n32))
    in_deg = zeros(Int32, Int(n32))
    @inbounds for a in 1:Int(m32)
        u = Int32(tail[a])
        v = Int32(head[a])
        (Int32(1) <= u <= n32 && Int32(1) <= v <= n32) ||
            throw(ArgumentError("arc $a has an endpoint outside 1:$n32"))
        tail32[a] = u
        head32[a] = v
        out_deg[u] += Int32(1)
        in_deg[v] += Int32(1)
    end

    first_out = Vector{Int32}(undef, Int(n32) + 1)
    first_in = Vector{Int32}(undef, Int(n32) + 1)
    first_out[1] = 1
    first_in[1] = 1
    @inbounds for v in 1:Int(n32)
        first_out[v + 1] = first_out[v] + out_deg[v]
        first_in[v + 1] = first_in[v] + in_deg[v]
    end

    out_arc = Vector{Int32}(undef, Int(m32))
    in_arc = Vector{Int32}(undef, Int(m32))
    next_out = Vector{Int32}(undef, Int(n32))
    next_in = Vector{Int32}(undef, Int(n32))
    @inbounds for v in 1:Int(n32)
        next_out[v] = first_out[v]
        next_in[v] = first_in[v]
    end
    @inbounds for a in Int32(1):m32
        u = tail32[a]
        v = head32[a]
        out_arc[next_out[u]] = a
        next_out[u] += Int32(1)
        in_arc[next_in[v]] = a
        next_in[v] += Int32(1)
    end

    labels = Vector{Int32}(undef, Int(n32))
    if label === nothing
        @inbounds for v in 1:Int(n32)
            labels[v] = Int32(v)
        end
    else
        length(label) == Int(n32) || throw(ArgumentError("label has the wrong length"))
        @inbounds for v in 1:Int(n32)
            labels[v] = Int32(label[v])
        end
    end
    return Graph(n32, m32, tail32, head32, first_out, out_arc, first_in, in_arc, labels)
end

function _vertex(g::Graph, v::Integer)::Int32
    iv = Int(v)
    (1 <= iv <= Int(g.n)) || throw(ArgumentError("vertex $v is outside 1:$(g.n)"))
    return Int32(iv)
end

"""
    out_degree(g, v) -> Int32

Number of arcs in the forward star of vertex `v`.
"""
function out_degree(g::Graph, v::Integer)
    vv = _vertex(g, v)
    return g.first_out[vv + Int32(1)] - g.first_out[vv]
end

"""
    in_degree(g, v) -> Int32

Number of arcs in the reverse star of vertex `v`.
"""
function in_degree(g::Graph, v::Integer)
    vv = _vertex(g, v)
    return g.first_in[vv + Int32(1)] - g.first_in[vv]
end

"""
    outgoing(g, v)

View of the arc indices in the forward star of vertex `v`.
"""
function outgoing(g::Graph, v::Integer)
    vv = _vertex(g, v)
    lo = Int(g.first_out[vv])
    hi = Int(g.first_out[vv + Int32(1)] - 1)
    return view(g.out_arc, lo:hi)
end

"""
    incoming(g, v)

View of the arc indices in the reverse star of vertex `v`.
"""
function incoming(g::Graph, v::Integer)
    vv = _vertex(g, v)
    lo = Int(g.first_in[vv])
    hi = Int(g.first_in[vv + Int32(1)] - 1)
    return view(g.in_arc, lo:hi)
end

"""
    is_dag(g) -> Bool

Return `true` if and only if `g` contains no directed cycle (including
self-loops).
"""
function is_dag(g::Graph)::Bool
    # Kahn's algorithm: repeatedly remove a vertex without incoming arcs.
    # The graph is acyclic iff every vertex can be removed this way. O(n + m).

    # A self-loop is a cycle; exit early (the peeling below would also fail
    # on it, since the loop keeps its vertex's in-degree above zero).
    @inbounds for a in 1:Int(g.m)
        if g.tail[a] == g.head[a]
            return false
        end
    end
    n = Int(g.n)
    # Number of arcs still entering each vertex (parallel arcs count separately).
    indeg = Vector{Int32}(undef, n)
    @inbounds for v in 1:n
        indeg[v] = in_degree(g, Int32(v))
    end
    # FIFO queue stored in an array of length n: `qh` is the index of the next
    # vertex to remove, `qt` the index of the last one inserted. Each vertex is
    # inserted at most once, so n slots suffice.
    queue = Vector{Int32}(undef, n)
    qh = 1
    qt = 0
    # Start with the sources, i.e. the vertices with no incoming arc.
    @inbounds for v in Int32(1):g.n
        if indeg[v] == 0
            qt += 1
            queue[qt] = v
        end
    end
    # `seen` counts the vertices removed so far.
    seen = 0
    while qh <= qt
        v = queue[qh]
        qh += 1
        seen += 1
        # Delete `v`: each outgoing arc no longer enters its head, and a head
        # left without incoming arcs becomes removable.
        @inbounds for a in outgoing(g, v)
            w = g.head[a]
            indeg[w] -= Int32(1)
            if indeg[w] == 0
                qt += 1
                queue[qt] = w
            end
        end
    end
    # Vertices on a directed cycle never reach in-degree 0, so they are never
    # removed: `seen < n` exactly when the graph has a cycle.
    return seen == n
end
