"""
    solve_rrsp(net, params, solver) -> Solution

Two-stage recoverable-robust shortest path: ``\\min_x C(x) + ADV(x)``.

- Interval and nominal uncertainty dispatch to [`solve_rec`](@ref).
- Continuous budget, ``k = 0``: two shortest paths, ``C(x)+\\hat{c}(x)+\\min(Γ,d(x))``.
- Continuous budget, ``Γ = 0``: recoverable problem under ``\\hat{c}``.
- Continuous budget, ``Γ`` at least ``\\sum d_a``: recoverable problem under ``\\hat{c}+d``.
- Continuous budget otherwise: Bold–Goerigk mixed-integer formulation, or the
  compact ASP MIP on inclusion when [`asp_decompose`](@ref) succeeds.
- Discrete budget, ``k = 0``: Bertsimas–Sim with first-stage cost in the
  objective. ``Δ = 0`` / ``Δ ≥ m`` reduce to recoverable under ``\\hat{c}`` /
  ``\\hat{c}+d``. Mid-budget discrete RRSP with ``k > 0`` returns
  `ST_NOT_IMPL` (no compact MIP); use [`approx_rrsp`](@ref) on DAGs.

`first` is the first-stage path and `second` a recovery. Combinatorial
solution of the general continuous-budget case returns `ST_NOT_IMPL`.

# Throws
- `ArgumentError`: if `gamma < 0` under a continuous budget, or if `delta < 0`
  under a discrete budget.
"""
function solve_rrsp(net::Network, params::Params, solver::Solver)::Solution
    U = params.uncertainty
    if U == U_INTERVAL || U == U_NOMINAL
        return solve_rec(net, params, solver)
    end
    t0 = time_ns()
    m = Int(net.graph.m)
    if U == U_DISC_BUDGET
        delta = _require_delta(params.delta)
        k = Int(params.k)
        if k == 0 || net.s == net.t
            return _solve_rrsp_k0_disc(net, delta, solver, t0)
        end
        if delta == 0
            return solve_rec(net, params_rrsp(U_NOMINAL, params.neighborhood, k), solver)
        end
        if delta >= m
            return solve_rec(net, params_rrsp(U_INTERVAL, params.neighborhood, k), solver)
        end
        uniq = _only_st_path(net)
        if uniq !== nothing
            return _solve_rrsp_k0_disc(net, delta, solver, t0)
        end
        meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
        return _not_impl(m, meth, t0)
    end
    gamma = _require_gamma(params.gamma)
    k = Int(params.k)
    if k == 0 || net.s == net.t
        return _solve_rrsp_k0_cont(net, gamma, solver, t0)
    end
    if gamma == 0.0
        return solve_rec(net, params_rrsp(U_NOMINAL, params.neighborhood, k), solver)
    end
    if _gamma_covers_all(net.costs.d, gamma)
        return solve_rec(net, params_rrsp(U_INTERVAL, params.neighborhood, k), solver)
    end
    uniq = _only_st_path(net)
    if uniq !== nothing
        return _solve_rrsp_k0_cont(net, gamma, solver, t0)
    end
    if solver.method == METHOD_COMB
        return _not_impl(m, METHOD_COMB, t0)
    end
    solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
    if params.neighborhood == NB_INCLUSION
        tree = _asp_tree_cached(net.graph, net.s, net.t)
        if tree !== nothing
            return _solve_rrsp_cont_asp(net, tree, params, solver, t0)
        end
    end
    return _solve_rrsp_cont_bg(net, params, solver, t0)
end
