# Exact RRSP on instances with a modest number of s–t paths.
#
# The commitment is chosen from an explicit set of paths. For each candidate
# the adversarial problem is solved by cut generation: a master problem over
# the attack (an LP for a continuous budget, a MIP with binary attack variables
# for a discrete budget) and, as separation oracle, the incremental shortest
# path problem for the current attack.

const _ENUM_BUDGETS = (U_CONT_BUDGET, U_DISC_BUDGET)
const _ENUM_POOL = 16    # attacks remembered as pruning certificates

"""
    solve_adv_cuts(net, params, x, solver; max_iter=10_000, tol=1e-7) -> Solution

Adversarial value of the committed path `x` under a continuous or discrete
budget: the maximum over attacks of the cheapest recovery in the neighborhood
of `x`,
``\\max_{δ} \\min_{y \\in N(x, k)} (\\hat c + δ)(y)``.

Works for every neighborhood and on any digraph for which
[`solve_inc`](@ref) works with `solver`. The master problem needs
`solver.optimizer`. The result's `z` is the value of the best attack found,
`second` is its cheapest recovery. Compare [`solve_adv`](@ref), which uses a
compact formulation.

Returns `ST_NOT_IMPL` for other uncertainty sets or without an optimizer, and
`ST_ERROR` if `max_iter` cuts are not enough.
"""
function solve_adv_cuts(
    net::Network,
    params::Params,
    x::Path,
    solver::Solver;
    max_iter::Integer = 10_000,
    tol::Real = 1e-7,
)::Solution
    t0 = time_ns()
    m = Int(net.graph.m)
    r = _adv_cuts(net, params, x, solver, max_iter, tol, Inf, ())
    r.status == ST_OK || return solution_empty(m; status = r.status, method = METHOD_MIP, time_sec = (time_ns() - t0) / 1e9)
    return _solution_adv(net, r.xp, r.y, r.value, t0, METHOD_MIP, 0.0)
end

# Worker of `solve_adv_cuts`. `seeds` are attack vectors (additional cost on
# every arc) tried first: each gives a lower bound on the adversarial value, and
# if some lower bound reaches `cutoff` the search stops early with `pruned`.
# The same early exit applies to the lower bounds met during cut generation.
function _adv_cuts(net::Network, params::Params, x::Path, solver::Solver, max_iter, tol, cutoff::Float64, seeds)
    g = net.graph
    m = Int(g.m)
    U = params.uncertainty
    nothing_result(st) = (status = st, xp = x, y = x, value = NaN, attack = Float64[], pruned = false)
    U in _ENUM_BUDGETS || return nothing_result(ST_NOT_IMPL)
    xp = _require_st_path(net, x)
    c_hat = net.costs.c_hat
    d = net.costs.d
    solver.optimizer === nothing && return nothing_result(ST_NOT_IMPL)
    cont = U == U_CONT_BUDGET
    cont ? _require_gamma(params.gamma) : _require_delta(params.delta)
    seen = Set{Vector{UInt8}}([copy(xp.chi)])
    sep_params = params_rrsp(U_INTERVAL, params.neighborhood, Int(params.k))
    best = -Inf
    best_y = _copy_path(xp)
    best_attack = zeros(Float64, m)
    w = Vector{Float64}(undef, m)
    result(st, pruned) = (status = st, xp = xp, y = best_y, value = best, attack = best_attack, pruned = pruned)
    function separate!(attack::AbstractVector{Float64})
        @inbounds for a in 1:m
            w[a] = c_hat[a] + attack[a]
        end
        sep = solve_inc(net, sep_params, xp, w, solver)
        sep.status == ST_OK || return sep
        if sep.z > best
            best = sep.z
            best_y = _copy_path(sep.second)
            best_attack = Vector{Float64}(attack)
        end
        return sep
    end
    # Seed attacks first: cheap lower bounds that may settle the cutoff before
    # any model is built.
    cuts = Path[xp]
    for seed in seeds
        sep = separate!(seed)
        sep.status == ST_OK || return result(sep.status, false)
        best >= cutoff && return result(ST_OK, true)
        if !(sep.second.chi in seen)
            push!(seen, copy(sep.second.chi))
            push!(cuts, _copy_path(sep.second))
        end
    end
    model = _new_mip(solver)
    if cont
        att = JuMP.@variable(model, [1:m], lower_bound = 0.0)
        @inbounds for a in 1:m
            JuMP.set_upper_bound(att[a], d[a])
        end
        JuMP.@constraint(model, sum(att) <= params.gamma)
        gain = att
    else
        att = JuMP.@variable(model, [1:m], Bin)
        JuMP.@constraint(model, sum(att) <= Float64(params.delta))
        gain = [d[a] * att[a] for a in 1:m]
    end
    tau = JuMP.@variable(model)
    JuMP.@objective(model, Max, tau)
    function add_cut!(y::Path)
        JuMP.@constraint(model, tau <= sum(c_hat[a] + gain[a] for a in y.seq; init = 0.0))
        return nothing
    end
    foreach(add_cut!, cuts)
    attack = Vector{Float64}(undef, m)
    for _ in 1:max_iter
        st, _gap = _optimize_mip(model)
        st == ST_OK || return result(st, false)
        upper = Float64(JuMP.objective_value(model))
        @inbounds for a in 1:m
            v = Float64(JuMP.value(att[a]))
            attack[a] = cont ? clamp(v, 0.0, d[a]) : (v >= 0.5 ? d[a] : 0.0)
        end
        sep = separate!(attack)
        sep.status == ST_OK || return result(sep.status, false)
        best >= cutoff && return result(ST_OK, true)
        upper <= best + Float64(tol) * max(1.0, abs(best)) && return result(ST_OK, false)
        sep.second.chi in seen && return result(ST_OK, false)
        push!(seen, copy(sep.second.chi))
        add_cut!(sep.second)
    end
    return result(ST_ERROR, false)
end

_require_budget(params::Params) =
    params.uncertainty == U_CONT_BUDGET ? _require_gamma(params.gamma) : _require_delta(params.delta)

# The adversary's natural attack on a path: spend the budget on its arcs with
# the largest deviations. Exact when no recovery is possible, and a good
# certificate otherwise.
function _greedy_attack(net::Network, params::Params, x::Path)::Vector{Float64}
    d = net.costs.d
    attack = zeros(Float64, length(d))
    order = sort(collect(x.seq); by = a -> -d[a])
    if params.uncertainty == U_CONT_BUDGET
        left = Float64(params.gamma)
        for a in order
            left <= 0.0 && break
            v = min(d[a], left)
            attack[a] = v
            left -= v
        end
    else
        for a in Iterators.take(order, Int(params.delta))
            attack[a] = d[a]
        end
    end
    return attack
end

"""
    enum_lower_bounds(net, params, paths, solver) -> Vector{Float64}

For each candidate commitment `x` in `paths`, the lower bound
``C(x) + \\min_{y \\in N(x, k)} \\hat c(y)`` on its recoverable-robust value:
the adversary can only raise costs above the nominal ones.
"""
function enum_lower_bounds(net::Network, params::Params, paths::AbstractVector{Path}, solver::Solver)::Vector{Float64}
    p = params_rrsp(U_INTERVAL, params.neighborhood, Int(params.k))
    c_hat = net.costs.c_hat
    lb = Vector{Float64}(undef, length(paths))
    for (i, x) in enumerate(paths)
        sol = solve_inc(net, p, x, c_hat, solver)
        sol.status == ST_OK || throw(ErrorException("recovery bound failed: $(sol.status)"))
        lb[i] = path_cost(x, net.costs.C) + sol.z
    end
    return lb
end

"""
    solve_rrsp_enum(net, params, solver; paths=nothing, lower_bounds=nothing,
                    limit=1_000_000, tol=1e-6) -> Solution

Recoverable-robust shortest path under a continuous or discrete budget by
enumerating the commitment ``x`` over the ``s``–``t`` paths (all of them by
default; pass `paths` to restrict the candidates).

Candidates are visited in increasing order of their lower bound
([`enum_lower_bounds`](@ref), or the supplied `lower_bounds`, one entry per
path). The search stops when the next lower bound reaches the incumbent, or
when the incumbent reaches the value that a recovery neighborhood containing
every path would give, which is a lower bound for all of them. Each visited
candidate is evaluated with [`solve_adv_cuts`](@ref), or in closed form for
``k = 0``.

This covers the discrete budget at positive ``k``, for which no compact
formulation is known, and serves as an independent check of the compact
formulations on small digraphs. It requires `solver.optimizer`.
"""
function solve_rrsp_enum(
    net::Network,
    params::Params,
    solver::Solver;
    paths::Union{Nothing,AbstractVector{Path}} = nothing,
    lower_bounds::Union{Nothing,AbstractVector{<:Real}} = nothing,
    limit::Integer = 1_000_000,
    tol::Real = 1e-6,
)::Solution
    t0 = time_ns()
    g = net.graph
    m = Int(g.m)
    U = params.uncertainty
    if U == U_INTERVAL || U == U_NOMINAL
        return solve_rec(net, params, solver)
    end
    U in _ENUM_BUDGETS || return _not_impl(m, METHOD_MIP, t0)
    ps = paths === nothing ? enumerate_st_paths(g, net.s, net.t; limit = limit) : paths
    isempty(ps) && return solution_empty(m; status = ST_INFEASIBLE, method = METHOD_MIP, time_sec = (time_ns() - t0) / 1e9)
    k = Int(params.k)
    C = net.costs.C
    Cx = [path_cost(x, C) for x in ps]
    best = Inf
    best_i = 0
    best_adv = nothing
    if k == 0
        # No recovery: the adversary's best attack is the greedy one.
        _require_budget(params)
        for (i, x) in enumerate(ps)
            atk = _greedy_attack(net, params, x)
            adv = _solution_adv(net, x, x, sum(net.costs.c_hat[a] + atk[a] for a in x.seq; init = 0.0), t0, METHOD_COMB, 0.0)
            z = Cx[i] + adv.z
            if z < best
                best, best_i, best_adv = z, i, adv
            end
        end
    else
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        if lower_bounds !== nothing
            length(lower_bounds) == length(ps) || throw(ArgumentError("lower_bounds has the wrong length"))
        end
        lb = lower_bounds === nothing ? enum_lower_bounds(net, params, ps, solver) : Float64.(lower_bounds)
        # The greedy attack on each candidate gives a stronger bound.
        sep_params = params_rrsp(U_INTERVAL, params.neighborhood, k)
        w = similar(net.costs.c_hat)
        lb = copy(lb)
        for (i, x) in enumerate(ps)
            atk = _greedy_attack(net, params, x)
            w .= net.costs.c_hat .+ atk
            sep = solve_inc(net, sep_params, x, w, solver)
            sep.status == ST_OK || return sep
            lb[i] = max(lb[i], Cx[i] + sep.z)
        end
        full = params_rrsp(params.uncertainty, NB_INCLUSION, 2 * Int(g.n); gamma = params.gamma, delta = params.delta)
        adv_full = solve_adv_cuts(net, full, ps[1], solver)
        adv_full.status == ST_OK || return adv_full
        floor_z = minimum(Cx) + adv_full.z
        pool = Vector{Vector{Float64}}()    # attacks that were best for evaluated candidates
        for i in sortperm(lb)
            lb[i] >= best - Float64(tol) && break
            # The candidate improves the incumbent only if its adversarial value is
            # below `cutoff`; any attack that reaches it is a certificate.
            cutoff = best_i == 0 ? Inf : (Cx[i] >= Cx[best_i] ? best - Float64(tol) : best + Float64(tol)) - Cx[i]
            seeds = Vector{Float64}[_greedy_attack(net, params, ps[i])]
            append!(seeds, pool)
            r = _adv_cuts(net, params, ps[i], solver, 10_000, 1e-7, cutoff, seeds)
            r.status == ST_OK || return solution_empty(m; status = r.status, method = METHOD_MIP, time_sec = (time_ns() - t0) / 1e9)
            r.pruned && continue
            pushfirst!(pool, r.attack)
            length(pool) > _ENUM_POOL && pop!(pool)
            adv = _solution_adv(net, r.xp, r.y, r.value, t0, METHOD_MIP, 0.0)
            z = Cx[i] + adv.z
            if best_i == 0 || z < best - Float64(tol)
                best, best_i, best_adv = z, i, adv
            elseif z <= best + Float64(tol) && Cx[i] < Cx[best_i]
                best, best_i, best_adv = min(best, z), i, adv
            end
            best <= floor_z + Float64(tol) && break
        end
    end
    best_i == 0 && return solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = (time_ns() - t0) / 1e9)
    x = ps[best_i]
    dt = (time_ns() - t0) / 1e9
    return Solution(
        _copy_path(x), _copy_path(best_adv.second), best, Cx[best_i], best - Cx[best_i], ST_OK, METHOD_MIP, dt, 0.0,
    )
end
