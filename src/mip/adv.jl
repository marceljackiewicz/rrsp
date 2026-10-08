function _delta_weights(c_hat::Vector{Float64}, delta, m::Int)::Vector{Float64}
    w = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        w[a] = c_hat[a] + Float64(JuMP.value(delta[a]))
        w[a] < 0 && (w[a] = 0.0)
    end
    return w
end

function _adv_from_lp(
    net::Network,
    x::Path,
    model::JuMP.Model,
    delta,
    nb::Neighborhood,
    k::Int,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if st != ST_OK
        return solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt)
    end
    z = Float64(JuMP.objective_value(model))
    w = _delta_weights(net.costs.c_hat, delta, m)
    y = _adv_recovery(net, x, w, nb, k, t0)
    return _solution_adv(net, x, y, z, t0, METHOD_MIP, gap)
end

function _solve_adv_incl_lp(
    net::Network,
    x::Path,
    k::Int,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    g = net.graph
    n = Int(g.n)
    m = Int(g.m)
    L = k
    c_hat = net.costs.c_hat
    d = net.costs.d
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    pi = JuMP.@variable(model, [1:n, 0:L])
    delta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    @inbounds for a in 1:m
        JuMP.set_upper_bound(delta[a], d[a])
        i = Int(g.tail[a])
        j = Int(g.head[a])
        if x.chi[a] != 0x00
            for l in 0:L
                JuMP.@constraint(model, pi[i, l] - pi[j, l] - delta[a] <= c_hat[a])
            end
        else
            for l in 0:(L - 1)
                JuMP.@constraint(model, pi[i, l] - pi[j, l + 1] - delta[a] <= c_hat[a])
            end
        end
    end
    @inbounds for v in 1:n
        for l in 0:(L - 1)
            JuMP.@constraint(model, pi[v, l] - pi[v, l + 1] <= 0)
        end
    end
    JuMP.@constraint(model, sum(delta[a] for a in 1:m) <= gamma)
    JuMP.@objective(model, Max, pi[Int(net.s), 0] - pi[Int(net.t), L])
    return _adv_from_lp(net, x, model, delta, NB_INCLUSION, L, t0)
end

function _solve_adv_excl_dag_lp(
    net::Network,
    x::Path,
    k::Int,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    g = net.graph
    n = Int(g.n)
    m = Int(g.m)
    nx = length(x.seq)
    # With k >= |x| every s-t path is a neighbour of x, so the adversary
    # maximises the shortest-path distance: the same LP with a single layer.
    R = max(nx - k, 0)
    c_hat = net.costs.c_hat
    d = net.costs.d
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    pi = JuMP.@variable(model, [1:n, 0:R])
    delta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    @inbounds for a in 1:m
        JuMP.set_upper_bound(delta[a], d[a])
        i = Int(g.tail[a])
        j = Int(g.head[a])
        if x.chi[a] != 0x00
            for l in 0:R
                lp = l + 1
                lp > R && (lp = R)
                JuMP.@constraint(model, pi[i, l] - pi[j, lp] - delta[a] <= c_hat[a])
            end
        else
            for l in 0:R
                JuMP.@constraint(model, pi[i, l] - pi[j, l] - delta[a] <= c_hat[a])
            end
        end
    end
    JuMP.@constraint(model, sum(delta[a] for a in 1:m) <= gamma)
    JuMP.@objective(model, Max, pi[Int(net.s), 0] - pi[Int(net.t), R])
    return _adv_from_lp(net, x, model, delta, NB_EXCLUSION, k, t0)
end

function _solve_adv_sym_dag_lp(
    net::Network,
    x::Path,
    k::Int,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    g = net.graph
    n = Int(g.n)
    m = Int(g.m)
    nx = length(x.seq)
    rmin = -nx
    rmax = k
    c_hat = net.costs.c_hat
    d = net.costs.d
    function idx(r::Int)::Int
        return r - rmin + 1
    end
    nR = rmax - rmin + 1
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    pi = JuMP.@variable(model, [1:n, 1:nR])
    ptau = JuMP.@variable(model)
    delta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    @inbounds for a in 1:m
        JuMP.set_upper_bound(delta[a], d[a])
        i = Int(g.tail[a])
        j = Int(g.head[a])
        if x.chi[a] != 0x00
            for r in (rmin + 1):rmax
                JuMP.@constraint(model, pi[i, idx(r)] - pi[j, idx(r - 1)] - delta[a] <= c_hat[a])
            end
        else
            for r in rmin:(rmax - 1)
                JuMP.@constraint(model, pi[i, idx(r)] - pi[j, idx(r + 1)] - delta[a] <= c_hat[a])
            end
        end
    end
    tmax = k - nx
    @inbounds for r in rmin:tmax
        JuMP.@constraint(model, pi[Int(net.t), idx(r)] - ptau <= 0)
    end
    JuMP.@constraint(model, sum(delta[a] for a in 1:m) <= gamma)
    JuMP.@objective(model, Max, pi[Int(net.s), idx(0)] - ptau)
    return _adv_from_lp(net, x, model, delta, NB_SYMDIFF, k, t0)
end

function _solve_adv_bg_mip(
    net::Network,
    x::Path,
    nb::Neighborhood,
    k::Int,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    g = net.graph
    m = Int(g.m)
    K = m + 1
    c_hat = net.costs.c_hat
    d = net.costs.d
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    ys = [JuMP.@variable(model, [1:m], Bin) for _ in 1:K]
    lam = JuMP.@variable(model, [1:K], lower_bound = 0.0)
    beta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    theta = JuMP.@variable(model, lower_bound = 0.0)
    zvar = [JuMP.@variable(model, [1:m], lower_bound = 0.0) for _ in 1:K]
    xfix = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        xfix[a] = Float64(x.chi[a])
    end
    dag = is_dag(g)
    @inbounds for i in 1:K
        add_path!(model, ys[i], g, net.s, net.t)
        if !dag && nb != NB_INCLUSION
            add_simple_path!(model, ys[i], g, net.s, net.t)
        end
        add_neighborhood!(model, xfix, ys[i], nb, k, m)
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
        sum(c_hat[a] * zvar[i][a] for i in 1:K, a in 1:m) + sum(d[a] * beta[a] for a in 1:m) + gamma * theta,
    )
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt)
    end
    z = Float64(JuMP.objective_value(model))
    best_i = 1
    best_l = -Inf
    @inbounds for i in 1:K
        lv = Float64(JuMP.value(lam[i]))
        if lv > best_l
            best_l = lv
            best_i = i
        end
    end
    yp = _path_from_mip(g, net.s, net.t, ys[best_i], m)
    yp === nothing && (yp = _copy_path(x))
    return _stamp_status(_solution_adv(net, x, yp, z, t0, METHOD_MIP, gap), st, gap)
end

function _solve_adv_cont_mip(
    net::Network,
    x::Path,
    nb::Neighborhood,
    k::Int,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    if nb == NB_INCLUSION
        return _solve_adv_incl_lp(net, x, k, gamma, solver, t0)
    end
    if is_dag(net.graph)
        if nb == NB_EXCLUSION
            return _solve_adv_excl_dag_lp(net, x, k, gamma, solver, t0)
        end
        return _solve_adv_sym_dag_lp(net, x, k, gamma, solver, t0)
    end
    return _solve_adv_bg_mip(net, x, nb, k, gamma, solver, t0)
end
