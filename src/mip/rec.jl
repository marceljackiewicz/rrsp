function _solve_rec_mip(
    net::Network,
    c2::Vector{Float64},
    nb::Neighborhood,
    k::Int,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    g = net.graph
    C = net.costs.C
    if net.s == net.t
        p = path_empty(m)
        return _solution_rec(net, p, p, c2, t0, METHOD_MIP, 0.0)
    end
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    x = JuMP.@variable(model, [1:m], Bin)
    y = JuMP.@variable(model, [1:m], Bin)
    add_path!(model, x, g, net.s, net.t)
    add_path!(model, y, g, net.s, net.t)
    dag = is_dag(g)
    if !dag
        add_simple_path!(model, x, g, net.s, net.t)
        if nb != NB_INCLUSION
            add_simple_path!(model, y, g, net.s, net.t)
        end
    end
    add_neighborhood!(model, x, y, nb, k, m)
    JuMP.@objective(model, Min, sum(C[a] * x[a] + c2[a] * y[a] for a in 1:m))
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return _with_mip_diag(
            solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt, mip_gap = gap),
            model,
        )
    end
    xp = _path_from_mip(g, net.s, net.t, x, m)
    yp = _path_from_mip(g, net.s, net.t, y, m)
    if xp === nothing || yp === nothing
        return _with_mip_diag(
            solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = dt),
            model,
        )
    end
    return _with_mip_diag(_stamp_status(_solution_rec(net, xp, yp, c2, t0, METHOD_MIP, gap), st, gap), model)
end
