function _solve_inc_mip(
    net::Network,
    x::Path,
    w::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    g = net.graph
    if net.s == net.t
        p = path_empty(m)
        return _solution_inc(net, p, p, w, t0, METHOD_MIP, 0.0)
    end
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    y = JuMP.@variable(model, [1:m], Bin)
    add_path!(model, y, g, net.s, net.t)
    if !is_dag(g)
        add_simple_path!(model, y, g, net.s, net.t)
    end
    xfix = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        xfix[a] = Float64(x.chi[a])
    end
    add_neighborhood!(model, xfix, y, nb, k, m)
    JuMP.@objective(model, Min, sum(w[a] * y[a] for a in 1:m))
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt)
    end
    yp = _path_from_mip(g, net.s, net.t, y, m)
    if yp === nothing
        return solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = dt)
    end
    return _stamp_status(_solution_inc(net, x, yp, w, t0, METHOD_MIP, gap), st, gap)
end
