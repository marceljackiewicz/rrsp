function _solution_adv(
    net::Network,
    x::Path,
    y::Path,
    z::Float64,
    t0::UInt64,
    method::Method,
    gap::Float64,
)::Solution
    dt = (time_ns() - t0) / 1e9
    z1 = path_cost(x, net.costs.C)
    return Solution(_copy_path(x), _copy_path(y), z, z1, z, ST_OK, method, dt, gap)
end

function _adv_recovery(
    net::Network,
    x::Path,
    w::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    t0::UInt64,
)::Path
    sol = _solve_inc_comb(net, x, w, nb, k, t0)
    sol.status == ST_OK && return sol.second
    return _copy_path(x)
end

function _solve_adv_dump(
    net::Network,
    x::Path,
    gamma::Float64,
    t0::UInt64,
)::Solution
    z = _dump_gamma(x, net.costs.c_hat, net.costs.d, gamma)
    return _solution_adv(net, x, x, z, t0, METHOD_COMB, 0.0)
end

function _solve_adv_interval(
    net::Network,
    params::Params,
    x::Path,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    w = Vector{Float64}(undef, m)
    c_hat = net.costs.c_hat
    d = net.costs.d
    @inbounds for a in 1:m
        w[a] = c_hat[a] + d[a]
    end
    _check_weights(w, m)
    sol = solve_inc(net, params, x, w, solver)
    sol.status != ST_OK && return sol
    return _solution_adv(net, x, sol.second, sol.z, t0, sol.method_used, sol.mip_gap)
end

function _solve_adv_nominal(
    net::Network,
    params::Params,
    x::Path,
    solver::Solver,
    t0::UInt64,
)::Solution
    w = _copy_weights(net.costs.c_hat, Int(net.graph.m))
    sol = solve_inc(net, params, x, w, solver)
    sol.status != ST_OK && return sol
    return _solution_adv(net, x, sol.second, sol.z, t0, sol.method_used, sol.mip_gap)
end

"""
    solve_adv(net, params, x, solver) -> Solution

Worst-case second-stage cost of a committed path `x`:
``\\max_{c \\in U} \\min_{y \\in N(x,k)} c(y)``.

- Interval uncertainty: incremental shortest path under ``\\hat{c}+d``.
- Nominal uncertainty: incremental shortest path under ``\\hat{c}``.
- Continuous budget, ``k = 0``: dump ``Γ`` onto `x` (greedy on ``d_a`` along `x`).
- Continuous budget, inclusion: Nasrabadi–Orlin LP on the time-expanded network.
- Continuous budget, exclusion / symmetric difference: LP on DAGs, mixed-integer
  on cyclic graphs.
- Discrete budget, ``k = 0``: dump the ``Δ`` largest deviations onto `x`.
- Discrete budget, ``Δ = 0`` or ``Δ ≥ m``: incremental under ``\\hat{c}`` or
  ``\\hat{c}+d``. Mid-budget discrete ADV with ``k > 0`` returns `ST_NOT_IMPL`.

`z` is the adversarial value. `first` is `x` and `second` is a recovery that
attains the inner minimum under a worst-case scenario.

# Throws
- `ArgumentError`: if `x` is not a simple ``s``–``t`` path, if `gamma < 0`
  under a continuous budget, or if `delta < 0` under a discrete budget.
"""
function solve_adv(net::Network, params::Params, x::Path, solver::Solver)::Solution
    t0 = time_ns()
    m = Int(net.graph.m)
    xp = _require_st_path(net, x)
    U = params.uncertainty
    if U == U_DISC_BUDGET
        delta = _require_delta(params.delta)
        k = Int(params.k)
        if k == 0 || net.s == net.t
            z = _dump_delta(xp, net.costs.c_hat, net.costs.d, delta)
            return _solution_adv(net, xp, xp, z, t0, METHOD_COMB, 0.0)
        end
        if delta == 0
            p2 = params_rrsp(U_NOMINAL, params.neighborhood, k)
            return _solve_adv_nominal(net, p2, xp, solver, t0)
        end
        if delta >= m
            p2 = params_rrsp(U_INTERVAL, params.neighborhood, k)
            return _solve_adv_interval(net, p2, xp, solver, t0)
        end
        meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
        return _not_impl(m, meth, t0)
    end
    if U == U_INTERVAL
        return _solve_adv_interval(net, params, xp, solver, t0)
    end
    if U == U_NOMINAL
        return _solve_adv_nominal(net, params, xp, solver, t0)
    end
    # U_CONT_BUDGET
    gamma = _require_gamma(params.gamma)
    k = Int(params.k)
    if k == 0 || net.s == net.t
        return _solve_adv_dump(net, xp, gamma, t0)
    end
    if _gamma_covers_all(net.costs.d, gamma)
        p2 = params_rrsp(U_INTERVAL, params.neighborhood, k)
        return _solve_adv_interval(net, p2, xp, solver, t0)
    end
    if gamma == 0.0
        p2 = params_rrsp(U_NOMINAL, params.neighborhood, k)
        return _solve_adv_nominal(net, p2, xp, solver, t0)
    end
    if solver.method == METHOD_COMB
        return _not_impl(m, METHOD_COMB, t0)
    end
    solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
    return _solve_adv_cont_mip(net, xp, params.neighborhood, k, gamma, solver, t0)
end

"""
    eval_worstcase(net, params, x, solver) -> Float64

Adversarial cost of a committed path: `solve_adv(...).z`.
"""
function eval_worstcase(
    net::Network,
    params::Params,
    x::Path,
    solver::Solver,
)::Float64
    return solve_adv(net, params, x, solver).z
end
