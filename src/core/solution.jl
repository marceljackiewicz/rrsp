"""
    solution_empty(m; status=ST_INFEASIBLE, method=METHOD_AUTO, time_sec=0) -> Solution

An empty result for a graph with `m` arcs: both paths empty, costs `Inf`.
"""
function solution_empty(
    m::Integer;
    status::Status = ST_INFEASIBLE,
    method::Method = METHOD_AUTO,
    time_sec::Real = 0.0,
    mip_gap::Real = 0.0,
)::Solution
    return Solution(path_empty(m), path_empty(m), Inf, Inf, Inf, status, method, Float64(time_sec), Float64(mip_gap))
end

"""
    Solver(; optimizer=nothing, method=METHOD_AUTO, time_limit=Inf,
            mip_gap=1e-6, threads=0, silent=true) -> Solver

Construct a solver handle. There is no process-wide optimizer.

`threads = 0` (the default) leaves the thread count to the MIP solver. Set a
positive value only if every model in the process uses the same value: HiGHS
keeps one global thread pool, and a model that asks for a different number of
threads than an earlier one fails with an error status.
"""
function Solver(;
    optimizer = nothing,
    method::Method = METHOD_AUTO,
    time_limit::Real = Inf,
    mip_gap::Real = 1e-6,
    threads::Integer = 0,
    silent::Bool = true,
)
    return Solver(optimizer, method, Float64(time_limit), Float64(mip_gap), Int32(threads), silent)
end
