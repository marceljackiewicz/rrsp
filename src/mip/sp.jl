function _solve_sp_mip(
    net::Network,
    w::Vector{Float64},
    solver::Solver,
    t0::UInt64;
    rob::Bool = false,
)::Solution
    m = Int(net.graph.m)
    g = net.graph
    if net.s == net.t
        return _finish_sp(net, ST_OK, path_empty(m), 0.0, t0; rob = rob, method = METHOD_MIP)
    end
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    x = JuMP.@variable(model, [1:m], Bin)
    add_path!(model, x, g, net.s, net.t)
    if !is_dag(g)
        add_simple_path!(model, x, g, net.s, net.t)
    end
    JuMP.@objective(model, Min, sum(w[a] * x[a] for a in 1:m))
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt)
    end
    path = _path_from_mip(g, net.s, net.t, x, m)
    if path === nothing
        return solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = dt)
    end
    z = path_cost(path, w)
    z1 = path_cost(path, net.costs.C)
    z2 = rob ? z : 0.0
    return _stamp_status(Solution(path, path_empty(m), z, z1, z2, ST_OK, METHOD_MIP, dt, gap), st, gap)
end
