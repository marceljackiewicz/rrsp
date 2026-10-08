function _del_id!(v::Vector{Int32}, id::Int32)
    @inbounds for i in 1:length(v)
        if v[i] == id
            v[i] = v[end]
            pop!(v)
            return
        end
    end
    return
end

# Arc sequences of the ASP dynamic programs are persistent ropes: joining two
# sequences is O(1) and only the sequence at the root is flattened. Copying the
# sequences at every merge would cost O(m) per merge, O(m²) per table.
struct _Rope
    len::Int
    leaf::Int32
    left::Union{Nothing,_Rope}
    right::Union{Nothing,_Rope}
end

const _ROPE_EMPTY = _Rope(0, Int32(0), nothing, nothing)

_rope_leaf(a::Integer) = _Rope(1, Int32(a), nothing, nothing)

function _rope_cat(a::_Rope, b::_Rope)::_Rope
    a.len == 0 && return b
    b.len == 0 && return a
    return _Rope(a.len + b.len, Int32(0), a, b)
end

"""The arcs of a rope, in order (iterative, so deep ropes do not overflow the stack)."""
function _rope_flatten(r::_Rope)::Vector{Int32}
    out = Vector{Int32}(undef, r.len)
    r.len == 0 && return out
    stack = _Rope[r]
    i = 0
    while !isempty(stack)
        node = pop!(stack)
        if node.left === nothing
            i += 1
            out[i] = node.leaf
        else
            push!(stack, node.right::_Rope)
            push!(stack, node.left::_Rope)
        end
    end
    return out
end

struct _AspPath
    cost::Float64
    seq::_Rope
end

const _ASP_EMPTY = _AspPath(Inf, _ROPE_EMPTY)

function _asp_cat(a::_AspPath, b::_AspPath)::_AspPath
    (a.cost == Inf || b.cost == Inf) && return _ASP_EMPTY
    return _AspPath(a.cost + b.cost, _rope_cat(a.seq, b.seq))
end

function _asp_better(a::_AspPath, b::_AspPath)::_AspPath
    return a.cost <= b.cost ? a : b
end

struct _AspPair
    xc::Float64
    yc::Float64
    xseq::_Rope
    yseq::_Rope
end

const _ASP_PAIR_EMPTY = _AspPair(Inf, Inf, _ROPE_EMPTY, _ROPE_EMPTY)

function _pair_tot(p::_AspPair)::Float64
    return p.xc + p.yc
end

function _pair_better(a::_AspPair, b::_AspPair)::_AspPair
    return _pair_tot(a) <= _pair_tot(b) ? a : b
end

function _pair_cat(a::_AspPair, b::_AspPair)::_AspPair
    (_pair_tot(a) == Inf || _pair_tot(b) == Inf) && return _ASP_PAIR_EMPTY
    return _AspPair(a.xc + b.xc, a.yc + b.yc, _rope_cat(a.xseq, b.xseq), _rope_cat(a.yseq, b.yseq))
end

function _pair_from(xp::_AspPath, yp::_AspPath)::_AspPair
    (xp.cost == Inf || yp.cost == Inf) && return _ASP_PAIR_EMPTY
    return _AspPair(xp.cost, yp.cost, xp.seq, yp.seq)
end

"""
    asp_decompose(g, s, t) -> Union{AspTree, Nothing}

Binary decomposition tree of the two-terminal multidigraph `(g, s, t)` if it is
arc-series-parallel, and `nothing` otherwise. Does not mutate `g`.
"""
function asp_decompose(g::Graph, s::Integer, t::Integer)::Union{AspTree,Nothing}
    ss = _vertex(g, s)
    tt = _vertex(g, t)
    ss == tt && return nothing
    m = Int(g.m)
    m == 0 && return nothing
    # An arc-series-parallel two-terminal graph is acyclic. Rejecting cycles up
    # front also keeps the reduction loop below finite (a self-loop is a series
    # pair with itself and would be requeued forever).
    is_dag(g) || return nothing
    n = Int(g.n)
    nodes = Vector{AspNode}(undef, m)
    @inbounds for a in 1:m
        nodes[a] = AspNode(ASP_LEAF, Int32(0), Int32(0), Int32(a), g.tail[a], g.head[a])
    end
    sizehint!(nodes, 2 * m - 1)
    out_ids = [Int32[] for _ in 1:n]
    in_ids = [Int32[] for _ in 1:n]
    @inbounds for a in Int32(1):Int32(m)
        push!(out_ids[Int(g.tail[a])], a)
        push!(in_ids[Int(g.head[a])], a)
    end

    function tail_of(id::Int32)::Int32
        return nodes[id].s
    end
    function head_of(id::Int32)::Int32
        return nodes[id].t
    end

    function add_live!(id::Int32)
        u = Int(tail_of(id))
        v = Int(head_of(id))
        push!(out_ids[u], id)
        push!(in_ids[v], id)
        return
    end

    function remove_live!(id::Int32)
        u = Int(tail_of(id))
        v = Int(head_of(id))
        _del_id!(out_ids[u], id)
        _del_id!(in_ids[v], id)
        return
    end

    q = Int32[]
    inq = zeros(UInt8, n)
    function requeue!(z::Int)
        if z != Int(ss) && z != Int(tt) && inq[z] == 0x00
            push!(q, Int32(z))
            inq[z] = 0x01
        end
        return
    end

    function parallel_reduce!(a::Int32, b::Int32)::Int32
        u = tail_of(a)
        v = head_of(a)
        remove_live!(a)
        remove_live!(b)
        nid = Int32(length(nodes) + 1)
        push!(nodes, AspNode(ASP_PARALLEL, a, b, Int32(0), u, v))
        add_live!(nid)
        requeue!(Int(u))
        requeue!(Int(v))
        return nid
    end

    function series_reduce!(ain::Int32, aout::Int32)::Int32
        u = tail_of(ain)
        w = head_of(aout)
        remove_live!(ain)
        remove_live!(aout)
        nid = Int32(length(nodes) + 1)
        push!(nodes, AspNode(ASP_SERIES, ain, aout, Int32(0), u, w))
        add_live!(nid)
        requeue!(Int(u))
        requeue!(Int(w))
        return nid
    end

    function reduce_parallel_in!(v::Int)
        changed = true
        while changed
            changed = false
            ins = in_ids[v]
            nvins = length(ins)
            @inbounds for i in 1:nvins
                for j in (i + 1):nvins
                    if tail_of(ins[i]) == tail_of(ins[j])
                        parallel_reduce!(ins[i], ins[j])
                        changed = true
                        break
                    end
                end
                changed && break
            end
        end
        return
    end

    function reduce_parallel_out!(v::Int)
        changed = true
        while changed
            changed = false
            outs = out_ids[v]
            nouts = length(outs)
            @inbounds for i in 1:nouts
                for j in (i + 1):nouts
                    if head_of(outs[i]) == head_of(outs[j])
                        parallel_reduce!(outs[i], outs[j])
                        changed = true
                        break
                    end
                end
                changed && break
            end
        end
        return
    end

    @inbounds for v in Int32(1):Int32(n)
        if v != ss && v != tt
            push!(q, v)
            inq[Int(v)] = 0x01
        end
    end
    qh = 1
    while qh <= length(q)
        v = Int(q[qh])
        qh += 1
        inq[v] = 0x00
        reduce_parallel_in!(v)
        reduce_parallel_out!(v)
        if length(in_ids[v]) == 1 && length(out_ids[v]) == 1 && in_ids[v][1] != out_ids[v][1]
            series_reduce!(in_ids[v][1], out_ids[v][1])
        end
    end

    live = Int32[]
    @inbounds for v in 1:n
        for id in out_ids[v]
            push!(live, id)
        end
    end
    isempty(live) && return nothing
    @inbounds for id in live
        if tail_of(id) != ss || head_of(id) != tt
            return nothing
        end
    end
    while length(out_ids[Int(ss)]) >= 2
        a = out_ids[Int(ss)][1]
        b = out_ids[Int(ss)][2]
        parallel_reduce!(a, b)
    end
    length(out_ids[Int(ss)]) == 1 || return nothing
    root = out_ids[Int(ss)][1]
    head_of(root) == tt || return nothing
    nn = Int32(length(nodes))
    nn == Int32(2 * m - 1) || return nothing
    return AspTree(nodes, root, nn)
end

# A `Graph` is immutable once built, so the decomposition of a given
# `(graph, s, t)` can be reused. The cache is keyed by the identity of the
# graph's `tail` array and dropped when the graph is garbage collected. This
# matters when many solves share one digraph (path enumeration, sampled
# scenarios).
const _ASP_CACHE = WeakKeyDict{Vector{Int32},Dict{Tuple{Int32,Int32},Union{Nothing,AspTree}}}()
const _ASP_LOCK = ReentrantLock()

function _asp_tree_cached(g::Graph, s::Integer, t::Integer)::Union{AspTree,Nothing}
    key = (Int32(s), Int32(t))
    lock(_ASP_LOCK) do
        per_graph = get!(_ASP_CACHE, g.tail) do
            Dict{Tuple{Int32,Int32},Union{Nothing,AspTree}}()
        end
        haskey(per_graph, key) && return per_graph[key]
        tree = asp_decompose(g, s, t)
        per_graph[key] = tree
        return tree
    end
end
