function _require_delta(delta::Int32)::Int
    d = Int(delta)
    d < 0 && throw(ArgumentError("delta < 0"))
    return d
end

"""
Worst-case second-stage cost of a committed path under a discrete budget
with no recovery: ``\\hat{c}(x)`` plus the ``Δ`` largest deviations on ``x``.
"""
function _dump_delta(x::Path, c_hat::Vector{Float64}, d::Vector{Float64}, delta::Int)::Float64
    z = 0.0
    ndev = 0
    seq = x.seq
    nseq = length(seq)
    if nseq == 0
        m = length(x.chi)
        @inbounds for a in 1:m
            x.chi[a] == 0x00 && continue
            z += c_hat[a]
            ndev += 1
        end
        ndev == 0 && return z
        devs = Vector{Float64}(undef, ndev)
        i = 0
        @inbounds for a in 1:m
            x.chi[a] == 0x00 && continue
            i += 1
            devs[i] = d[a]
        end
        sort!(devs; rev = true)
        k = delta < ndev ? delta : ndev
        @inbounds for j in 1:k
            z += devs[j]
        end
        return z
    end
    @inbounds for a in seq
        z += c_hat[a]
    end
    delta <= 0 && return z
    order = Vector{Float64}(undef, nseq)
    @inbounds for j in 1:nseq
        order[j] = d[Int(seq[j])]
    end
    sort!(order; rev = true)
    k = delta < nseq ? delta : nseq
    @inbounds for j in 1:k
        z += order[j]
    end
    return z
end

function _threshold_deviations(d::Vector{Float64})::Vector{Float64}
    vals = copy(d)
    push!(vals, 0.0)
    unique!(vals)
    sort!(vals; rev = true)
    return vals
end

function _bertsimas_sim_comb(
    net::Network,
    base::Vector{Float64},
    c_hat::Vector{Float64},
    d::Vector{Float64},
    delta::Int,
)::Tuple{Status,Path,Float64}
    m = Int(net.graph.m)
    g = net.graph
    if net.s == net.t
        return (ST_OK, path_empty(m), 0.0)
    end
    w = Vector{Float64}(undef, m)
    best_z = Inf
    best_st = ST_INFEASIBLE
    best_path = path_empty(m)
    saw_ok = false
    saw_infeas = false
    for dl in _threshold_deviations(d)
        @inbounds for a in 1:m
            extra = d[a] - dl
            w[a] = base[a] + c_hat[a] + (extra > 0.0 ? extra : 0.0)
        end
        _check_weights(w, m)
        st, path, _zsp = _shortest_path(g, net.s, net.t, w)
        if st == ST_INFEASIBLE
            saw_infeas = true
            continue
        elseif st != ST_OK
            return (st, path, Inf)
        end
        saw_ok = true
        z = path_cost(path, base) + _dump_delta(path, c_hat, d, delta)
        if z < best_z
            best_z = z
            best_st = ST_OK
            best_path = path
        end
    end
    saw_ok && return (best_st, best_path, best_z)
    saw_infeas && return (ST_INFEASIBLE, best_path, Inf)
    return (ST_ERROR, best_path, Inf)
end

function _zeros_m(m::Int)::Vector{Float64}
    return zeros(Float64, m)
end

function _solve_rob_disc(
    net::Network,
    delta::Int,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    c_hat = net.costs.c_hat
    d = net.costs.d
    if net.s == net.t
        meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
        return _finish_sp(net, ST_OK, path_empty(m), 0.0, t0; rob = true, method = meth)
    end
    if delta == 0
        w = _copy_weights(c_hat, m)
        if solver.method == METHOD_MIP
            solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
            return _solve_sp_mip(net, w, solver, t0; rob = true)
        end
        st, path, z = _shortest_path(net.graph, net.s, net.t, w)
        return _finish_sp(net, st, path, z, t0; rob = true)
    end
    if delta >= m
        w = Vector{Float64}(undef, m)
        @inbounds for a in 1:m
            w[a] = c_hat[a] + d[a]
        end
        _check_weights(w, m)
        if solver.method == METHOD_MIP
            solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
            return _solve_sp_mip(net, w, solver, t0; rob = true)
        end
        st, path, z = _shortest_path(net.graph, net.s, net.t, w)
        return _finish_sp(net, st, path, z, t0; rob = true)
    end
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        return _solve_rob_disc_mip(net, _zeros_m(m), delta, solver, t0; rob = true)
    end
    st, path, z = _bertsimas_sim_comb(net, _zeros_m(m), c_hat, d, delta)
    return _finish_sp(net, st, path, z, t0; rob = true, method = METHOD_COMB)
end

function _solve_rrsp_k0_disc(
    net::Network,
    delta::Int,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    C = net.costs.C
    c_hat = net.costs.c_hat
    d = net.costs.d
    meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
    if net.s == net.t
        p = path_empty(m)
        return Solution(p, p, 0.0, 0.0, 0.0, ST_OK, meth, (time_ns() - t0) / 1e9, 0.0)
    end
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        return _solve_rob_disc_mip(net, C, delta, solver, t0; rob = false)
    end
    st, path, z = _bertsimas_sim_comb(net, C, c_hat, d, delta)
    dt = (time_ns() - t0) / 1e9
    st != ST_OK && return solution_empty(m; status = st, method = METHOD_COMB, time_sec = dt)
    z1 = path_cost(path, C)
    return Solution(_copy_path(path), _copy_path(path), z, z1, z - z1, ST_OK, METHOD_COMB, dt, 0.0)
end
