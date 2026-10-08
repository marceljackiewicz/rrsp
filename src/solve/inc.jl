function _require_st_path(net::Network, x::Path)::Path
    m = Int(net.graph.m)
    length(x.chi) == m || throw(ArgumentError("path does not match network"))
    return path_from_chi(net.graph, net.s, net.t, x.chi)
end

"""
    solve_inc(net, params, x, w, solver) -> Solution

Incremental shortest path: a minimum-weight ``s``–``t`` path ``y`` in the
neighborhood of the committed path `x`.

The objective is ``z = w(y)``. The recovered path is returned in `second`;
`first` is a copy of `x`. `z_second` equals `z`.

- Inclusion is combinatorial on every digraph (constrained shortest path /
  time-expanded network, or the ASP decomposition-tree DP).
- Exclusion and symmetric difference are combinatorial on DAGs and on
  arc-series-parallel graphs, and mixed-integer on cyclic instances.
- ``k = 0`` returns ``y = x`` without a MIP.

# Throws
- `ArgumentError`: if `x` is not a simple ``s``–``t`` path, if `w` has length
  other than ``m``, or if some weight is negative.
"""
function solve_inc(
    net::Network,
    params::Params,
    x::Path,
    w::AbstractVector{<:Real},
    solver::Solver,
)::Solution
    t0 = time_ns()
    m = Int(net.graph.m)
    xp = _require_st_path(net, x)
    k = Int(params.k)
    nb = params.neighborhood
    g = net.graph
    dag = is_dag(g)
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        ww = _copy_weights(w, m)
        return _solve_inc_mip(net, xp, ww, nb, k, solver, t0)
    end
    if solver.method == METHOD_COMB
        if nb != NB_INCLUSION && !dag && k > 0
            return _not_impl(m, METHOD_COMB, t0)
        end
        ww = _copy_weights(w, m)
        return _solve_inc_comb(net, xp, ww, nb, k, t0)
    end
    # METHOD_AUTO
    if nb == NB_INCLUSION || dag || k == 0
        ww = _copy_weights(w, m)
        return _solve_inc_comb(net, xp, ww, nb, k, t0)
    end
    solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
    ww = _copy_weights(w, m)
    return _solve_inc_mip(net, xp, ww, nb, k, solver, t0)
end

"""
    eval_recovered(net, params, x, w, solver) -> Float64

Incremental cost of recovering from `x` under weights `w`: `solve_inc(...).z`.
"""
function eval_recovered(
    net::Network,
    params::Params,
    x::Path,
    w::AbstractVector{<:Real},
    solver::Solver,
)::Float64
    return solve_inc(net, params, x, w, solver).z
end
