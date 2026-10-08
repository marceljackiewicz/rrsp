"""
    costs_new(C, c_hat, d) -> Costs

Construct arc costs, copying the three arrays.

# Throws
- `ArgumentError`: if the arrays differ in length, or if some cost or
  deviation is negative or not finite.
"""
function costs_new(
    C::AbstractVector{<:Real},
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
)::Costs
    length(C) == length(c_hat) == length(d) ||
        throw(ArgumentError("Costs length mismatch"))
    m = length(C)
    Cc = Vector{Float64}(undef, m)
    ch = Vector{Float64}(undef, m)
    dd = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        da = Float64(d[a])
        ca = Float64(C[a])
        cha = Float64(c_hat[a])
        (isfinite(ca) && ca >= 0) || throw(ArgumentError("first-stage cost of arc $a is negative or not finite"))
        (isfinite(cha) && cha >= 0) || throw(ArgumentError("nominal cost of arc $a is negative or not finite"))
        (isfinite(da) && da >= 0) || throw(ArgumentError("deviation of arc $a is negative or not finite"))
        Cc[a] = ca
        ch[a] = cha
        dd[a] = da
    end
    return Costs(Cc, ch, dd)
end

"""
    network_new(graph, s, t, costs) -> Network

Construct a shortest-path instance. The cost arrays are copied.

# Throws
- `ArgumentError`: if `s` or `t` lies outside `1, …, n`, or if the cost
  arrays have length other than `m`.
"""
function network_new(graph::Graph, s::Integer, t::Integer, costs::Costs)::Network
    ss = _vertex(graph, s)
    tt = _vertex(graph, t)
    length(costs.C) == Int(graph.m) ||
        throw(ArgumentError("Costs.m != Graph.m"))
    copied = costs_new(costs.C, costs.c_hat, costs.d)
    return Network(graph, ss, tt, copied)
end
