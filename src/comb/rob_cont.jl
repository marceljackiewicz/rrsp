function _sum_vec(v::Vector{Float64})::Float64
    s = 0.0
    @inbounds for i in 1:length(v)
        s += v[i]
    end
    return s
end

function _require_gamma(gamma::Float64)::Float64
    gamma < 0.0 && throw(ArgumentError("gamma < 0"))
    return gamma
end

function _gamma_covers_all(d::Vector{Float64}, gamma::Float64)::Bool
    return gamma >= _sum_vec(d)
end

"""
Worst-case second-stage cost of a committed path under a continuous budget
with no recovery: ``\\hat{c}(x) + \\min(Γ, d(x))``.
"""
function _dump_gamma(x::Path, c_hat::Vector{Float64}, d::Vector{Float64}, gamma::Float64)::Float64
    rem = gamma
    z = 0.0
    seq = x.seq
    nseq = length(seq)
    if nseq == 0
        m = length(x.chi)
        @inbounds for a in 1:m
            x.chi[a] == 0x00 && continue
            take = d[a] < rem ? d[a] : rem
            rem -= take
            z += c_hat[a] + take
        end
        return z
    end
    order = collect(Int, seq)
    sort!(order; by = a -> d[a], rev = true)
    @inbounds for a in order
        take = d[a] < rem ? d[a] : rem
        rem -= take
        z += c_hat[a] + take
    end
    return z
end

function _rob_cont_cost(x::Path, c_hat::Vector{Float64}, d::Vector{Float64}, gamma::Float64)::Float64
    return _dump_gamma(x, c_hat, d, gamma)
end

function _rrsp_k0_cost(x::Path, C::Vector{Float64}, c_hat::Vector{Float64}, d::Vector{Float64}, gamma::Float64)::Float64
    return path_cost(x, C) + _dump_gamma(x, c_hat, d, gamma)
end

function _pick_better_rob(
    p1::Path,
    z1::Float64,
    st1::Status,
    p2::Path,
    z2::Float64,
    st2::Status,
)::Tuple{Status,Path,Float64}
    ok1 = st1 == ST_OK
    ok2 = st2 == ST_OK
    if ok1 && ok2
        return z1 <= z2 ? (ST_OK, p1, z1) : (ST_OK, p2, z2)
    elseif ok1
        return (ST_OK, p1, z1)
    elseif ok2
        return (ST_OK, p2, z2)
    elseif st1 == ST_INFEASIBLE && st2 == ST_INFEASIBLE
        return (ST_INFEASIBLE, p1, Inf)
    else
        return (st1 != ST_OK ? st1 : st2, p1, Inf)
    end
end

function _solve_rob_cont(
    net::Network,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    c_hat = net.costs.c_hat
    d = net.costs.d
    g = net.graph
    if net.s == net.t
        p = path_empty(m)
        return _finish_sp(net, ST_OK, p, 0.0, t0; rob = true, method = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB)
    end
    wnom = Vector{Float64}(undef, m)
    wmax = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        wnom[a] = c_hat[a]
        wmax[a] = c_hat[a] + d[a]
    end
    _check_weights(wnom, m)
    _check_weights(wmax, m)
    use_mip = solver.method == METHOD_MIP
    if use_mip
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        soln = _solve_sp_mip(net, wnom, solver, t0; rob = false)
        solx = _solve_sp_mip(net, wmax, solver, t0; rob = false)
        z1 = soln.status == ST_OK ? _rob_cont_cost(soln.first, c_hat, d, gamma) : Inf
        z2 = solx.status == ST_OK ? _rob_cont_cost(solx.first, c_hat, d, gamma) : Inf
        st, path, z = _pick_better_rob(soln.first, z1, soln.status, solx.first, z2, solx.status)
        gap = max(soln.mip_gap, solx.mip_gap)
        return _finish_sp(net, st, path, z, t0; rob = true, method = METHOD_MIP, gap = gap)
    end
    stn, pn, _zn = _shortest_path(g, net.s, net.t, wnom)
    stx, px, _zx = _shortest_path(g, net.s, net.t, wmax)
    z1 = stn == ST_OK ? _rob_cont_cost(pn, c_hat, d, gamma) : Inf
    z2 = stx == ST_OK ? _rob_cont_cost(px, c_hat, d, gamma) : Inf
    st, path, z = _pick_better_rob(pn, z1, stn, px, z2, stx)
    return _finish_sp(net, st, path, z, t0; rob = true, method = METHOD_COMB)
end

function _only_st_path(net::Network)::Union{Nothing,Path}
    m = Int(net.graph.m)
    m == 0 && return nothing
    w = zeros(Float64, m)
    st, p, _z = _shortest_path(net.graph, net.s, net.t, w)
    st == ST_OK && length(p.seq) == m || return nothing
    return p
end

function _solve_rrsp_k0_cont(
    net::Network,
    gamma::Float64,
    solver::Solver,
    t0::UInt64,
)::Solution
    m = Int(net.graph.m)
    C = net.costs.C
    c_hat = net.costs.c_hat
    d = net.costs.d
    g = net.graph
    meth = solver.method == METHOD_MIP ? METHOD_MIP : METHOD_COMB
    if net.s == net.t
        p = path_empty(m)
        return Solution(p, p, 0.0, 0.0, 0.0, ST_OK, meth, (time_ns() - t0) / 1e9, 0.0)
    end
    wnom = Vector{Float64}(undef, m)
    wmax = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        wnom[a] = C[a] + c_hat[a]
        wmax[a] = C[a] + c_hat[a] + d[a]
    end
    _check_weights(wnom, m)
    _check_weights(wmax, m)
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        soln = _solve_sp_mip(net, wnom, solver, t0; rob = false)
        solx = _solve_sp_mip(net, wmax, solver, t0; rob = false)
        z1 = soln.status == ST_OK ? _rrsp_k0_cost(soln.first, C, c_hat, d, gamma) : Inf
        z2 = solx.status == ST_OK ? _rrsp_k0_cost(solx.first, C, c_hat, d, gamma) : Inf
        st, path, z = _pick_better_rob(soln.first, z1, soln.status, solx.first, z2, solx.status)
        dt = (time_ns() - t0) / 1e9
        st != ST_OK && return solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt)
        z1c = path_cost(path, C)
        z2c = z - z1c
        gap = max(soln.mip_gap, solx.mip_gap)
        return Solution(_copy_path(path), _copy_path(path), z, z1c, z2c, ST_OK, METHOD_MIP, dt, gap)
    end
    stn, pn, _zn = _shortest_path(g, net.s, net.t, wnom)
    stx, px, _zx = _shortest_path(g, net.s, net.t, wmax)
    z1 = stn == ST_OK ? _rrsp_k0_cost(pn, C, c_hat, d, gamma) : Inf
    z2 = stx == ST_OK ? _rrsp_k0_cost(px, C, c_hat, d, gamma) : Inf
    st, path, z = _pick_better_rob(pn, z1, stn, px, z2, stx)
    dt = (time_ns() - t0) / 1e9
    st != ST_OK && return solution_empty(m; status = st, method = METHOD_COMB, time_sec = dt)
    z1c = path_cost(path, C)
    return Solution(_copy_path(path), _copy_path(path), z, z1c, z - z1c, ST_OK, METHOD_COMB, dt, 0.0)
end
