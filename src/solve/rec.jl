function _second_stage_weights(net::Network, params::Params)::Union{Nothing,Vector{Float64}}
    m = Int(net.graph.m)
    if params.uncertainty == U_INTERVAL
        w = Vector{Float64}(undef, m)
        c_hat = net.costs.c_hat
        d = net.costs.d
        @inbounds for a in 1:m
            w[a] = c_hat[a] + d[a]
        end
        _check_weights(w, m)
        return w
    elseif params.uncertainty == U_NOMINAL
        return _copy_weights(net.costs.c_hat, m)
    else
        return nothing
    end
end

function _solve_rec_k0(net::Network, c2::Vector{Float64}, t0::UInt64, method::Method)::Solution
    m = Int(net.graph.m)
    w = Vector{Float64}(undef, m)
    C = net.costs.C
    @inbounds for a in 1:m
        w[a] = C[a] + c2[a]
    end
    _check_weights(w, m)
    status, path, _z = _shortest_path(net.graph, net.s, net.t, w)
    dt = (time_ns() - t0) / 1e9
    if status != ST_OK
        return solution_empty(m; status = status, method = method, time_sec = dt)
    end
    return _solution_rec(net, path, path, c2, t0, method, 0.0)
end

"""
    solve_rec(net, params, solver) -> Solution

Recoverable shortest path: ``\\min C(x) + c_2(y)`` over ``s``–``t`` paths with
``y`` in the neighborhood of ``x``.

Under interval uncertainty, ``c_2 = \\hat{c}+d``. Under nominal uncertainty,
``c_2 = \\hat{c}``. Continuous and discrete budgets return `ST_NOT_IMPL`.

The first-stage path is in `first`, the recovery in `second`. `z_first` is
``C(x)`` and `z_second` is ``c_2(y)``.

Combinatorial solution is used on DAGs and arc-series-parallel graphs
(decomposition-tree DP when an [`AspTree`](@ref) exists;
otherwise the DAG dynamic network) and when ``k = 0`` (shortest path under
``C+c_2``). Cyclic instances with ``k > 0`` use a mixed-integer formulation.

# Throws
- `ArgumentError`: if some first-stage or second-stage weight used by the
  combinatorial routine is negative.
"""
function solve_rec(net::Network, params::Params, solver::Solver)::Solution
    t0 = time_ns()
    m = Int(net.graph.m)
    c2 = _second_stage_weights(net, params)
    if c2 === nothing
        meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
        return _not_impl(m, meth, t0)
    end
    return _solve_rec_given_c2(net, c2, params.neighborhood, Int(params.k), solver, t0)
end

function _solve_rec_given_c2(
    net::Network,
    c2::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    dag = is_dag(net.graph)
    k = _clamp_k(net.graph, k)
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        return _solve_rec_mip(net, c2, nb, k, solver, t0)
    end
    if solver.method == METHOD_COMB
        if !dag && k > 0
            return _not_impl(m, METHOD_COMB, t0)
        end
        if k == 0 || net.s == net.t
            return _solve_rec_k0(net, c2, t0, METHOD_COMB)
        end
        tree = _asp_tree_cached(net.graph, net.s, net.t)
        if tree !== nothing
            return _solve_rec_asp(net, tree, c2, nb, k, t0)
        end
        return _solve_rec_dag(net, c2, nb, k, t0)
    end
    # METHOD_AUTO
    if dag || k == 0 || net.s == net.t
        if k == 0 || net.s == net.t
            return _solve_rec_k0(net, c2, t0, METHOD_COMB)
        end
        tree = _asp_tree_cached(net.graph, net.s, net.t)
        if tree !== nothing
            return _solve_rec_asp(net, tree, c2, nb, k, t0)
        end
        return _solve_rec_dag(net, c2, nb, k, t0)
    end
    solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
    return _solve_rec_mip(net, c2, nb, k, solver, t0)
end
