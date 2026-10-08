function _solve_rob_disc_mip(
    net::Network,
    base::Vector{Float64},
    delta::Int,
    solver::Solver,
    t0::UInt64;
    rob::Bool,
)::Solution
    m = Int(net.graph.m)
    g = net.graph
    c_hat = net.costs.c_hat
    d = net.costs.d
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    x = JuMP.@variable(model, [1:m], Bin)
    add_path!(model, x, g, net.s, net.t)
    if !is_dag(g)
        add_simple_path!(model, x, g, net.s, net.t)
    end
    theta = JuMP.@variable(model, lower_bound = 0.0)
    p = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    @inbounds for a in 1:m
        JuMP.@constraint(model, theta + p[a] >= d[a] * x[a])
    end
    JuMP.@objective(
        model,
        Min,
        sum((base[a] + c_hat[a]) * x[a] for a in 1:m) + Float64(delta) * theta + sum(p[a] for a in 1:m)
    )
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt)
    end
    path = _path_from_mip(g, net.s, net.t, x, m)
    if path === nothing
        return solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = dt)
    end
    z = path_cost(path, base) + _dump_delta(path, c_hat, d, delta)
    z1 = path_cost(path, net.costs.C)
    if rob
        return _stamp_status(Solution(path, path_empty(m), z, z1, z, ST_OK, METHOD_MIP, dt, gap), st, gap)
    end
    return _stamp_status(
        Solution(_copy_path(path), _copy_path(path), z, z1, z - z1, ST_OK, METHOD_MIP, dt, gap), st, gap,
    )
end
