"""
    cost_structure_alpha(net) -> Union{Nothing,Float64}
    cost_structure_alpha(costs) -> Union{Nothing,Float64}

Largest ``α ∈ (0, 1]`` such that ``\\hat{c}_a \\ge α(\\hat{c}_a + d_a)`` for
every arc, or `nothing` if some arc has zero nominal second-stage cost and
positive deviation (the ``α``-condition fails).
"""
function cost_structure_alpha(c::Costs)::Union{Nothing,Float64}
    m = length(c.c_hat)
    a = Inf
    c_hat = c.c_hat
    d = c.d
    @inbounds for i in 1:m
        mx = c_hat[i] + d[i]
        mx <= 0.0 && continue
        c_hat[i] <= 0.0 && return nothing
        r = c_hat[i] / mx
        r < a && (a = r)
    end
    a == Inf && return 1.0
    return a
end

cost_structure_alpha(net::Network) = cost_structure_alpha(net.costs)

function _cmax_weights(net::Network)::Vector{Float64}
    m = Int(net.graph.m)
    w = Vector{Float64}(undef, m)
    c_hat = net.costs.c_hat
    d = net.costs.d
    @inbounds for a in 1:m
        w[a] = c_hat[a] + d[a]
    end
    _check_weights(w, m)
    return w
end

function _eval_bar_F(
    net::Network,
    params::Params,
    x::Path,
    solver::Solver,
    t0::UInt64,
)::Solution
    w = _cmax_weights(net)
    p_int = params_rrsp(U_INTERVAL, params.neighborhood, Int(params.k))
    inc = solve_inc(net, p_int, x, w, solver)
    dt = (time_ns() - t0) / 1e9
    inc.status != ST_OK && return inc
    z1 = path_cost(x, net.costs.C)
    z2 = inc.z
    return Solution(
        _copy_path(x),
        _copy_path(inc.second),
        z1 + z2,
        z1,
        z2,
        ST_OK,
        inc.method_used,
        dt,
        0.0,
    )
end

function _eval_heuristic_path(
    net::Network,
    params::Params,
    x::Path,
    solver::Solver,
    t0::UInt64,
)::Solution
    adv = solve_adv(net, params, x, solver)
    if adv.status == ST_OK
        z1 = path_cost(x, net.costs.C)
        dt = (time_ns() - t0) / 1e9
        return Solution(
            _copy_path(x),
            _copy_path(adv.second),
            z1 + adv.z,
            z1,
            adv.z,
            ST_OK,
            adv.method_used,
            dt,
            adv.mip_gap,
        )
    end
    adv.status == ST_NOT_IMPL && return _eval_bar_F(net, params, x, solver, t0)
    return adv
end

function _better_approx(a::Solution, b::Solution)::Solution
    a.status != ST_OK && return b
    b.status != ST_OK && return a
    return a.z <= b.z ? a : b
end

"""
    approx_rrsp(net, params, solver) -> Solution

Polynomial heuristic for budgeted recoverable-robust shortest path
(combinatorial Rec requires a DAG when ``k > 0``), from the single-scenario bounds of Hradovich, Kasperski, and Zieliński.

The first-stage path is an optimum of the nominal recoverable problem
``\\min \\hat{F}(X)``. Under a continuous budget a second candidate is the
optimum under the proportional scenario ``S'`` with costs
``\\hat{c} + κ \\min(d, Γ)``, where ``κ = Γ / D`` and
``D = \\sum_a \\min(d_a, Γ)`` (see [`cost_structure_kappa`](@ref)); it is tried
only if ``0 < Γ < \\sum_a d_a``. `z` is a polynomial-time upper
bound on the recoverable-robust optimum: the true adversarial cost of the
chosen path when that evaluation is implemented, otherwise the interval
evaluation ``\\overline{F}``.

If a cost-structure constant ``α`` exists ([`cost_structure_alpha`](@ref)),
the returned value satisfies ``z \\le \\mathrm{OPT}/α``. Under a continuous
budget the same construction is an ``m``-approximation, and a ``1/κ``
approximation for ``κ = Γ/D`` (with ``D`` as above) when ``0 < Γ < \\sum_a d_a``.

Interval and nominal uncertainty dispatch to [`solve_rec`](@ref) (exact).
Combinatorial Rec on a cyclic graph with ``k > 0`` returns `ST_NOT_IMPL`.
"""
function approx_rrsp(net::Network, params::Params, solver::Solver)::Solution
    U = params.uncertainty
    if U == U_INTERVAL || U == U_NOMINAL
        return solve_rec(net, params, solver)
    end
    t0 = time_ns()
    m = Int(net.graph.m)
    k = Int(params.k)
    nb = params.neighborhood
    if U == U_DISC_BUDGET
        _require_delta(params.delta)
    elseif U == U_CONT_BUDGET
        _require_gamma(params.gamma)
    else
        meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
        return _not_impl(m, meth, t0)
    end
    rec_hat = solve_rec(net, params_rrsp(U_NOMINAL, nb, k), solver)
    rec_hat.status != ST_OK && return rec_hat
    best = _eval_heuristic_path(net, params, rec_hat.first, solver, t0)
    if U == U_CONT_BUDGET
        gamma = params.gamma
        D = _sum_vec(net.costs.d)
        if gamma > 0.0 && D > 0.0 && gamma < D
            d = net.costs.d
            c_hat = net.costs.c_hat
            dcap = Vector{Float64}(undef, m)
            Dcap = 0.0
            @inbounds for a in 1:m
                da = d[a] < gamma ? d[a] : gamma
                dcap[a] = da
                Dcap += da
            end
            kappa = gamma / Dcap
            c2 = Vector{Float64}(undef, m)
            @inbounds for a in 1:m
                c2[a] = c_hat[a] + kappa * dcap[a]
            end
            rec_p = _solve_rec_given_c2(net, c2, nb, k, solver, t0)
            if rec_p.status == ST_OK
                ev = _eval_heuristic_path(net, params, rec_p.first, solver, t0)
                best = _better_approx(best, ev)
            elseif best.status != ST_OK
                return rec_p
            end
        end
    end
    return best
end

"""
    cost_structure_kappa(net, gamma) -> Union{Nothing,Float64}
    cost_structure_kappa(costs, gamma) -> Union{Nothing,Float64}

Continuous-budget fraction ``κ = Γ / D`` with ``D = ∑_a min(d_a, Γ)``.
`nothing` if ``Γ ≤ 0``, if there is no positive deviation mass, or if
``Γ ≥ \\sum_a d_a`` (the budget covers all deviations, so the proportional
scenario coincides with the interval scenario and no factor is needed).
"""
function cost_structure_kappa(c::Costs, gamma::Real)::Union{Nothing,Float64}
    g = Float64(gamma)
    g <= 0.0 && return nothing
    D = 0.0
    Draw = 0.0
    @inbounds for da in c.d
        Draw += da
        D += da < g ? da : g
    end
    D <= 0.0 && return nothing
    g >= Draw && return nothing
    return g / D
end

cost_structure_kappa(net::Network, gamma::Real) = cost_structure_kappa(net.costs, gamma)

"""
    approx_bound_factors(net, params, solver=nothing) -> NamedTuple

Theoretical approximation factors from the uncertainty chapter: ``1/α`` when
the cost-structure condition holds, and for a continuous budget also ``1/κ``
and ``1/(1-ν)`` (the last needs `solver` to evaluate ``F_c(X')``). The
unconditional continuous factor never exceeds ``m``.
"""
function approx_bound_factors(
    net::Network,
    params::Params,
    solver::Union{Nothing,Solver} = nothing,
)
    alpha = cost_structure_alpha(net)
    m = Float64(net.graph.m)
    kappa = nothing
    nu = nothing
    if params.uncertainty == U_CONT_BUDGET
        kappa = cost_structure_kappa(net, params.gamma)
        if solver !== nothing && kappa !== nothing && 0.0 < kappa < 1.0
            c2 = Vector{Float64}(undef, Int(net.graph.m))
            c_hat = net.costs.c_hat
            d = net.costs.d
            @inbounds for a in 1:Int(net.graph.m)
                da = d[a]
                dcap = da < params.gamma ? da : params.gamma
                c2[a] = c_hat[a] + kappa * dcap
            end
            rec_p = _solve_rec_given_c2(
                net, c2, params.neighborhood, Int(params.k), solver, time_ns(),
            )
            if rec_p.status == ST_OK
                fc = path_cost(rec_p.first, net.costs.C) +
                    eval_worstcase(net, params, rec_p.first, solver)
                if isfinite(fc) && fc > 0.0
                    nu = params.gamma / fc
                end
            end
        end
    end
    b_alpha = alpha === nothing ? NaN : 1.0 / alpha
    b_kappa = kappa === nothing || kappa <= 0.0 ? NaN : 1.0 / kappa
    b_nu = (nu === nothing || !(nu < 1.0) || nu < 0.0) ? NaN : 1.0 / (1.0 - nu)
    b_m = params.uncertainty == U_CONT_BUDGET ? m : NaN
    bound = Inf
    for b in (b_alpha, b_kappa, b_nu, b_m)
        isfinite(b) && b < bound && (bound = b)
    end
    return (
        alpha = alpha === nothing ? NaN : alpha,
        kappa = kappa === nothing ? NaN : kappa,
        nu = nu === nothing ? NaN : nu,
        bound_alpha = b_alpha,
        bound_kappa = b_kappa,
        bound_nu = b_nu,
        bound_m = b_m,
        bound = isfinite(bound) ? bound : NaN,
    )
end
