# Build the RRSP solution from the commitment `x` found by a MIP.
#
# The MIP objective `z_mip` is the value reported. The adversarial problem is
# solved for `x` as an independent check: if it fails, its status is
# propagated; if it disagrees with `z_mip` for a proven-optimal solve, the
# result is reported as an error rather than silently replaced.
function _rrsp_from_x(
    net::Network,
    params::Params,
    x::Path,
    z_mip::Float64,
    gap::Float64,
    st::Status,
    solver::Solver,
    t0::UInt64,
)::Solution
    dt = (time_ns() - t0) / 1e9
    z1 = path_cost(x, net.costs.C)
    adv = solve_adv(net, params, x, solver)
    if adv.status != ST_OK
        return Solution(_copy_path(x), _copy_path(x), z_mip, z1, z_mip - z1, adv.status, METHOD_MIP, dt, gap)
    end
    z_eval = z1 + adv.z
    tol = max(1e-6, 10 * solver.mip_gap) * max(1.0, abs(z_eval))
    if st == ST_OK && abs(z_eval - z_mip) > tol
        @warn "RRSP MIP value and adversarial evaluation of its commitment disagree" z_mip z_eval maxlog = 3
        return Solution(_copy_path(x), _copy_path(adv.second), z_eval, z1, adv.z, ST_ERROR, METHOD_MIP, dt, gap)
    end
    # At a limit the incumbent need not be optimal: report its exact value.
    z = st == ST_OK ? z_mip : z_eval
    return _stamp_status(
        Solution(_copy_path(x), _copy_path(adv.second), z, z1, z - z1, ST_OK, METHOD_MIP, dt, gap),
        st, gap,
    )
end

function _solve_rrsp_cont_bg(
    net::Network,
    params::Params,
    solver::Solver,
    t0::UInt64,
)::Solution
    g = net.graph
    m = Int(g.m)
    K = m + 1
    C = net.costs.C
    c_hat = net.costs.c_hat
    d = net.costs.d
    gamma = params.gamma
    nb = params.neighborhood
    k = Int(params.k)
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    x = JuMP.@variable(model, [1:m], Bin)
    ys = [JuMP.@variable(model, [1:m], Bin) for _ in 1:K]
    lam = JuMP.@variable(model, [1:K], lower_bound = 0.0)
    beta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    theta = JuMP.@variable(model, lower_bound = 0.0)
    zvar = [JuMP.@variable(model, [1:m], lower_bound = 0.0) for _ in 1:K]
    dag = is_dag(g)
    add_path!(model, x, g, net.s, net.t)
    if !dag
        add_simple_path!(model, x, g, net.s, net.t)
    end
    @inbounds for i in 1:K
        add_path!(model, ys[i], g, net.s, net.t)
        if !dag && nb != NB_INCLUSION
            add_simple_path!(model, ys[i], g, net.s, net.t)
        end
        add_neighborhood!(model, x, ys[i], nb, k, m)
        for a in 1:m
            JuMP.@constraint(model, zvar[i][a] <= ys[i][a])
            JuMP.@constraint(model, zvar[i][a] <= lam[i])
            JuMP.@constraint(model, zvar[i][a] >= lam[i] - (1.0 - ys[i][a]))
        end
    end
    JuMP.@constraint(model, sum(lam[i] for i in 1:K) == 1)
    @inbounds for a in 1:m
        JuMP.@constraint(model, beta[a] + theta >= sum(zvar[i][a] for i in 1:K))
    end
    JuMP.@objective(
        model,
        Min,
        sum(C[a] * x[a] for a in 1:m) +
        sum(c_hat[a] * zvar[i][a] for i in 1:K, a in 1:m) +
        sum(d[a] * beta[a] for a in 1:m) +
        gamma * theta,
    )
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return _with_mip_diag(
            solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt, mip_gap = gap),
            model,
        )
    end
    xp = _path_from_mip(g, net.s, net.t, x, m)
    xp === nothing && return _with_mip_diag(
        solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = dt),
        model,
    )
    z_mip = Float64(JuMP.objective_value(model))
    return _with_mip_diag(_rrsp_from_x(net, params, xp, z_mip, gap, st, solver, t0), model)
end
