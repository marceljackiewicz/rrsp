# Budgets at which the adversary can no longer be recovered from (continuous-budget and discrete-budget experiments).
#
# For zero first-stage costs, let Z_int be the interval-robust value
# min_x (ĉ + d)(x). A budget "forces" the interval value when one attack
# raises every s–t path to at least Z_int. Beyond that budget every
# commitment has the same adversarial value, so the budget axes of the figures
# are fractions of it.

using JuMP

"""Nominal and interval costs of every path, and the derived constants."""
struct PathCosts
    nominal::Vector{Float64}     # ĉ(p)
    interval::Vector{Float64}    # (ĉ + d)(p)
    z_nom::Float64               # min ĉ(p)
    z_int::Float64               # min (ĉ + d)(p)
end

function path_costs(paths::AbstractVector{Path}, c_hat::AbstractVector{Float64}, d::AbstractVector{Float64})
    nominal = [sum(c_hat[a] for a in p.seq; init = 0.0) for p in paths]
    interval = [nominal[i] + sum(d[a] for a in paths[i].seq; init = 0.0) for i in eachindex(paths)]
    return PathCosts(nominal, interval, minimum(nominal), minimum(interval))
end

"""
    saturation_budget(pc) -> Float64

Largest single-path shortfall `Z_int − ĉ(p)`: a budget below it cannot raise
every path to `Z_int` when one attack gives it all to the cheapest path.
"""
function saturation_budget(pc::PathCosts)::Float64
    g = 0.0
    for c in pc.nominal
        c < pc.z_int - 1e-9 && (g = max(g, pc.z_int - c))
    end
    return g
end

"""
    continuous_forcing_budget(paths, pc, c_hat, d, optimizer) -> Float64

The least continuous budget Γ for which some attack `0 ≤ δ ≤ d`, `sum(δ) ≤ Γ`
gives every path cost at least `Z_int`: the optimum of a covering LP, solved
by adding the most violated path constraint until none is violated.
"""
function continuous_forcing_budget(
    paths::AbstractVector{Path},
    pc::PathCosts,
    d::AbstractVector{Float64},
    optimizer;
    tol::Real = 1e-6,
)::Float64
    need = [max(pc.z_int - c, 0.0) for c in pc.nominal]
    active = [i for i in eachindex(paths) if need[i] > tol]
    isempty(active) && return 0.0
    m = length(d)
    model = JuMP.Model(optimizer)
    JuMP.set_silent(model)
    delta = JuMP.@variable(model, [a = 1:m], lower_bound = 0.0, upper_bound = d[a])
    JuMP.@objective(model, Min, sum(delta))
    seed = active[argmax([need[i] for i in active])]
    added = Set{Int}([seed])
    JuMP.@constraint(model, sum(delta[a] for a in paths[seed].seq) >= need[seed])
    for _ in 1:(length(active) + 1)
        JuMP.optimize!(model)
        JuMP.termination_status(model) == JuMP.MOI.OPTIMAL || error("forcing-budget LP: $(JuMP.termination_status(model))")
        vals = JuMP.value.(delta)
        worst, worst_slack = 0, 0.0
        for i in active
            slack = need[i] - sum(vals[a] for a in paths[i].seq)
            if slack > worst_slack
                worst, worst_slack = i, slack
            end
        end
        worst_slack <= 1e-5 && return JuMP.objective_value(model)
        worst in added && error("forcing-budget LP: constraint repeated (slack $worst_slack)")
        push!(added, worst)
        JuMP.@constraint(model, sum(delta[a] for a in paths[worst].seq) >= need[worst])
    end
    error("forcing-budget LP did not converge")
end

"""
    discrete_forcing_budget(net, paths, pc, solver, floor_of) -> Int

The least number Δ of arcs to raise so that every path costs at least `Z_int`:
a binary search on the adversarial value over all paths (`floor_of(Δ)`),
which is nondecreasing in Δ.
"""
function discrete_forcing_budget(net::Network, pc::PathCosts, floor_of::Function; tol::Real = 1e-4)::Int
    floor_of(0) >= pc.z_int - tol && return 0
    lo, hi = 0, Int(net.graph.m)
    while lo < hi
        mid = (lo + hi) >>> 1
        if floor_of(mid) >= pc.z_int - tol
            hi = mid
        else
            lo = mid + 1
        end
    end
    return lo
end
