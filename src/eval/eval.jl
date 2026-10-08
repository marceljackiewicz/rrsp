"""
    eval_nominal(net, x) -> Float64

Nominal second-stage cost ``\\hat{c}(x)`` of path `x` on `net`.

# Throws
- `ArgumentError`: if `x.chi` has length other than `m`.
"""
function eval_nominal(net::Network, x::Path)::Float64
    length(x.chi) == Int(net.graph.m) ||
        throw(ArgumentError("path does not match network"))
    return path_cost(x, net.costs.c_hat)
end

"""
    eval_max(net, x) -> Float64

Worst-case interval cost ``(\\hat{c}+d)(x)`` of path `x` on `net`.

# Throws
- `ArgumentError`: if `x.chi` has length other than `m`.
"""
function eval_max(net::Network, x::Path)::Float64
    length(x.chi) == Int(net.graph.m) ||
        throw(ArgumentError("path does not match network"))
    z = 0.0
    c_hat = net.costs.c_hat
    d = net.costs.d
    @inbounds for a in 1:length(x.chi)
        if x.chi[a] != 0x00
            z += c_hat[a] + d[a]
        end
    end
    return z
end
