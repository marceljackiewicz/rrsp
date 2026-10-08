function add_path!(model::JuMP.Model, x, g::Graph, s::Int32, t::Int32)
    n = Int(g.n)
    @inbounds for v in Int32(1):Int32(n)
        inflow = JuMP.AffExpr(0.0)
        for a in incoming(g, v)
            JuMP.add_to_expression!(inflow, x[a])
        end
        outflow = JuMP.AffExpr(0.0)
        for a in outgoing(g, v)
            JuMP.add_to_expression!(outflow, x[a])
        end
        if v == s && v == t
            JuMP.@constraint(model, inflow == outflow)
            JuMP.@constraint(model, outflow == 0)
        elseif v == s
            JuMP.@constraint(model, outflow - inflow == 1)
        elseif v == t
            JuMP.@constraint(model, inflow - outflow == 1)
        else
            JuMP.@constraint(model, inflow == outflow)
        end
    end
    return nothing
end

"""
    add_simple_path!(model, x, g, s, t)

Miller–Tucker–Zemlin ordering constraints so that `x` is a *simple* path.
The ordering variables are continuous, as required by the thesis (not integer).
"""
function add_simple_path!(model::JuMP.Model, x, g::Graph, s::Int32, t::Int32)
    n = Int(g.n)
    m = Int(g.m)
    p = JuMP.@variable(model, [1:n], lower_bound = 0.0)
    nn = Float64(n)
    @inbounds for a in 1:m
        i = Int(g.tail[a])
        j = Int(g.head[a])
        JuMP.@constraint(model, p[i] - p[j] + nn * x[a] <= nn - 1.0)
    end
    if s != t
        for a in incoming(g, s)
            JuMP.@constraint(model, x[a] == 0)
        end
        for a in outgoing(g, t)
            JuMP.@constraint(model, x[a] == 0)
        end
    end
    return nothing
end

function add_neighborhood!(model::JuMP.Model, x, y, nb::Neighborhood, k::Integer, m::Int)
    z = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    @inbounds for a in 1:m
        JuMP.@constraint(model, z[a] <= x[a])
        JuMP.@constraint(model, z[a] <= y[a])
    end
    kf = Float64(k)
    if nb == NB_INCLUSION
        JuMP.@constraint(model, sum(y[a] - z[a] for a in 1:m) <= kf)
    elseif nb == NB_EXCLUSION
        JuMP.@constraint(model, sum(x[a] - z[a] for a in 1:m) <= kf)
    else
        JuMP.@constraint(model, sum(y[a] + x[a] - 2 * z[a] for a in 1:m) <= kf)
    end
    return nothing
end

# MOI errors that only mean "this optimizer does not support that".
const _MOI_SOFT_ERRORS = Union{JuMP.MOI.UnsupportedError,JuMP.MOI.NotAllowedError}

# Name of the relative MIP gap option, by `JuMP.solver_name`. `MOI.RelativeGap`
# itself is a result attribute and cannot be set. Solvers not listed here run
# with their own default gap.
const _GAP_ATTRIBUTE = Dict{String,String}(
    "HiGHS" => "mip_rel_gap",
    "Gurobi" => "MIPGap",
    "CPLEX" => "CPXPARAM_MIP_Tolerances_MIPGap",
    "SCIP" => "limits/gap",
    "Cbc" => "ratioGap",
    "GLPK" => "mip_gap",
    "COPT" => "RelGap",
    "Xpress" => "MIPRELSTOP",
)

function _set_optional_attribute(model::JuMP.Model, attr, value)
    try
        JuMP.set_attribute(model, attr, value)
    catch e
        e isa _MOI_SOFT_ERRORS || rethrow()
    end
    return nothing
end

function _new_mip(solver::Solver)
    solver.optimizer === nothing && return nothing
    model = JuMP.Model(solver.optimizer)
    if solver.silent
        JuMP.set_silent(model)
    end
    if isfinite(solver.time_limit)
        JuMP.set_time_limit_sec(model, solver.time_limit)
    end
    gap_attr = get(_GAP_ATTRIBUTE, JuMP.solver_name(model), nothing)
    if gap_attr !== nothing
        _set_optional_attribute(model, JuMP.MOI.RawOptimizerAttribute(gap_attr), solver.mip_gap)
    end
    if solver.threads > 0
        _set_optional_attribute(model, JuMP.MOI.NumberOfThreads(), Int(solver.threads))
    end
    return model
end

function _mip_status(model::JuMP.Model)::Status
    st = JuMP.termination_status(model)
    st == JuMP.MOI.OPTIMAL && return ST_OK
    st == JuMP.MOI.ALMOST_OPTIMAL && return ST_OK
    st == JuMP.MOI.INFEASIBLE && return ST_INFEASIBLE
    st == JuMP.MOI.INFEASIBLE_OR_UNBOUNDED && return ST_INFEASIBLE
    st == JuMP.MOI.DUAL_INFEASIBLE && return ST_UNBOUNDED
    st == JuMP.MOI.TIME_LIMIT && return ST_TIME_LIMIT
    st == JuMP.MOI.NODE_LIMIT && return ST_TIME_LIMIT
    @warn "MIP solver stopped with status $st" raw_status = JuMP.raw_status(model) maxlog = 3
    return ST_ERROR
end

function _mip_gap(model::JuMP.Model)::Float64
    try
        g = JuMP.relative_gap(model)
        return isfinite(g) ? Float64(g) : 0.0
    catch e
        # An LP has no MIP gap.
        (e isa _MOI_SOFT_ERRORS || e isa ErrorException) || rethrow()
        return 0.0
    end
end

# `true` if the solve produced a usable point: optimal, or stopped at a limit
# with a feasible incumbent.
function _mip_usable(model::JuMP.Model, st::Status)::Bool
    st == ST_OK && return true
    st == ST_TIME_LIMIT || return false
    return JuMP.primal_status(model) == JuMP.MOI.FEASIBLE_POINT
end

function _chi_from_vars(x, m::Int)::Vector{UInt8}
    chi = Vector{UInt8}(undef, m)
    @inbounds for a in 1:m
        chi[a] = JuMP.value(x[a]) > 0.5 ? 0x01 : 0x00
    end
    return chi
end

# Extract the s–t path from a MIP solution. If the selected arcs also contain
# a circulation (a zero-cost cycle that the model does not forbid), erase the
# loops while walking from `s`; for nonnegative costs this does not increase
# the cost or leave the neighborhood.
function _path_from_mip(g::Graph, s::Int32, t::Int32, x, m::Int)::Union{Path,Nothing}
    chi = _chi_from_vars(x, m)
    try
        return path_from_chi(g, s, t, chi)
    catch e
        e isa ArgumentError || rethrow()
    end
    return _loop_erased_path(g, s, t, chi)
end

function _loop_erased_path(g::Graph, s::Int32, t::Int32, chi::Vector{UInt8})::Union{Path,Nothing}
    n = Int(g.n)
    pos = zeros(Int, n)
    used = falses(Int(g.m))
    stack = Int32[s]
    seq = Int32[]
    pos[Int(s)] = 1
    v = s
    while v != t
        nxt = Int32(0)
        for a in outgoing(g, v)
            if chi[a] != 0x00 && !used[a]
                nxt = a
                break
            end
        end
        nxt == 0 && return nothing
        used[nxt] = true
        w = g.head[nxt]
        if pos[Int(w)] != 0
            k = pos[Int(w)]
            for u in stack[(k + 1):end]
                pos[Int(u)] = 0
            end
            resize!(stack, k)
            resize!(seq, k - 1)
        else
            push!(stack, w)
            push!(seq, nxt)
            pos[Int(w)] = length(stack)
        end
        v = w
    end
    return path_from_seq(g, s, t, seq)
end

function _not_impl(m::Int, method::Method, t0::UInt64)::Solution
    dt = (time_ns() - t0) / 1e9
    return solution_empty(m; status = ST_NOT_IMPL, method = method, time_sec = dt)
end

function _optimize_mip(model::JuMP.Model)::Tuple{Status,Float64}
    try
        JuMP.optimize!(model)
    catch e
        e isa _MOI_SOFT_ERRORS || rethrow()
        @warn "MIP solver cannot handle the model" exception = e maxlog = 3
        return (ST_NOT_IMPL, 0.0)
    end
    return (_mip_status(model), _mip_gap(model))
end

# Result of a solve that stopped at a limit with an incumbent keeps the status.
function _stamp_status(sol::Solution, st::Status, gap::Float64)::Solution
    st == ST_OK && return sol
    return Solution(
        sol.first, sol.second, sol.z, sol.z_first, sol.z_second, st, sol.method_used,
        sol.time_sec, gap, sol.mip_nodes, sol.n_binaries, sol.n_constraints,
    )
end

function _mip_diag(model::JuMP.Model)::Tuple{Int32,Int32,Int32}
    n_bin = 0
    try
        for v in JuMP.all_variables(model)
            try
                JuMP.is_binary(v) && (n_bin += 1)
            catch
            end
        end
    catch
    end
    n_cons = 0
    try
        n_cons = JuMP.num_constraints(model; count_variable_in_set_constraints = false)
    catch
        try
            n_cons = JuMP.num_constraints(model)
        catch
        end
    end
    n_nodes = 0
    try
        n_nodes = Int(JuMP.MOI.get(model, JuMP.MOI.NodeCount()))
    catch
        try
            n_nodes = Int(JuMP.MOI.get(JuMP.unsafe_backend(model), JuMP.MOI.NodeCount()))
        catch
        end
    end
    return (Int32(n_nodes), Int32(n_bin), Int32(n_cons))
end

function _with_mip_diag(sol::Solution, model::JuMP.Model)::Solution
    nodes, nbin, ncons = _mip_diag(model)
    return Solution(
        sol.first,
        sol.second,
        sol.z,
        sol.z_first,
        sol.z_second,
        sol.status,
        sol.method_used,
        sol.time_sec,
        sol.mip_gap,
        nodes,
        nbin,
        ncons,
    )
end

function _copy_path(p::Path)::Path
    return Path(copy(p.chi), copy(p.seq))
end

function _solution_inc(
    net::Network,
    x::Path,
    y::Path,
    w::Vector{Float64},
    t0::UInt64,
    method::Method,
    gap::Float64,
)::Solution
    dt = (time_ns() - t0) / 1e9
    z = path_cost(y, w)
    z1 = path_cost(x, net.costs.C)
    return Solution(_copy_path(x), _copy_path(y), z, z1, z, ST_OK, method, dt, gap)
end

function _solution_rec(
    net::Network,
    x::Path,
    y::Path,
    c2::Vector{Float64},
    t0::UInt64,
    method::Method,
    gap::Float64,
)::Solution
    dt = (time_ns() - t0) / 1e9
    z1 = path_cost(x, net.costs.C)
    z2 = path_cost(y, c2)
    return Solution(_copy_path(x), _copy_path(y), z1 + z2, z1, z2, ST_OK, method, dt, gap)
end
