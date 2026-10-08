# From the raw result CSVs to the curves and table entries of the thesis.
#
# A *curve* is a piecewise linear function of the relative recovery budget
# x = k/ℓ (k/(3ℓ) for the arc-series-parallel family), given by its break
# points. The value of recovery of an instance is
#     VoR(k) = 100 (Z(0) − Z(k)) / Z(0).
# Fixed families average the cost draws of their one digraph. The random-DAG
# family averages the cost draws of each digraph first and then the digraphs,
# after bringing every digraph to the common grid j/12 by linear interpolation
# (the digraphs have different ℓ).

const GRID_STEPS = 12
const GRID = [j / GRID_STEPS for j in 0:GRID_STEPS]
const FAMILY_ORDER = ("layered", "random_dag", "grid", "asp")
const FAMILY_TITLE = Dict("layered" => "layered", "random_dag" => "random DAG", "grid" => "grid", "asp" => "ASP")
const FAMILY_ROW = Dict("layered" => "Layered", "random_dag" => "Random DAG", "grid" => "Grid", "asp" => "ASP")

struct Curve
    x::Vector{Float64}
    y::Vector{Float64}
end

"""Value of the piecewise linear curve at `t` (constant outside its range)."""
function interpolate(c::Curve, t::Real)::Float64
    x, y = c.x, c.y
    t <= x[1] && return y[1]
    t >= x[end] && return y[end]
    i = searchsortedlast(x, t)
    x[i] == t && return y[i]
    w = (t - x[i]) / (x[i + 1] - x[i])
    return (1 - w) * y[i] + w * y[i + 1]
end

on_grid(c::Curve) = Curve(copy(GRID), [interpolate(c, t) for t in GRID])

"""Value of recovery in percent from the values `z` for `k = 0, 1, ...`."""
function recovery_percent(z::AbstractVector{<:Real})::Vector{Float64}
    z0 = z[1]
    z0 > 0 || return zeros(length(z))
    return [100 * (z0 - v) / z0 for v in z]
end

"""Budget unit of an instance: 3ℓ for the ASP family, ℓ otherwise."""
budget_unit(family::AbstractString, ell::Integer) = family == "asp" ? 3 * ell : ell

"""One instance (a digraph and a cost draw) with a VoR curve."""
struct InstanceCurve
    family::String
    graph::Int
    cost::Int
    curve::Curve
end

mean_vec(vs) = sum(vs) / length(vs)

"""Pointwise mean of curves that share their break points."""
function mean_curve(cs::AbstractVector{Curve})::Curve
    x = cs[1].x
    all(c -> c.x == x, cs) || error("curves with different break points")
    return Curve(copy(x), [mean_vec([c.y[i] for c in cs]) for i in eachindex(x)])
end

"""
    family_curve(ics; grid) -> Curve

Average VoR curve of one family from its instance curves. With `grid = true`
(or for the random-DAG family, always) the result is on the grid `j/12`.
"""
function family_curve(ics::AbstractVector{InstanceCurve}; grid::Bool = false)::Curve
    family = ics[1].family
    graphs = sort(unique(ic.graph for ic in ics))
    per_graph = Curve[]
    for gi in graphs
        cs = [ic.curve for ic in ics if ic.graph == gi]
        push!(per_graph, mean_curve(cs))
    end
    if family == "random_dag"
        return mean_curve(on_grid.(per_graph))
    end
    length(per_graph) == 1 || error("$family: expected one digraph")
    return grid ? on_grid(per_graph[1]) : per_graph[1]
end

# --- grouping raw rows ----------------------------------------------------------

"""Group `rows` by the instance key `(family, graph, cost)`, in sorted order."""
function by_instance(rows)
    groups = Dict{Tuple{String,Int,Int},Vector{eltype(rows)}}()
    for r in rows
        push!(get!(groups, (r.family, r.graph, r.cost), eltype(rows)[]), r)
    end
    return [k => sort(groups[k]; by = r -> r.k) for k in sort(collect(keys(groups)); by = k -> (findfirst(==(k[1]), FAMILY_ORDER), k[2], k[3]))]
end

check_ks(rs) = (collect(r.k for r in rs) == collect(0:length(rs) - 1) || error("missing recovery budgets for $(rs[1].family)"))

# --- interval experiment -------------------------------------------------------------

"""VoR curves of the interval experiment (independent costs)."""
function interval_instance_curves(rows)::Vector{InstanceCurve}
    out = InstanceCurve[]
    for ((family, graph, cost), rs) in by_instance(rows)
        check_ks(rs)
        unit = budget_unit(family, rs[1].ell)
        push!(out, InstanceCurve(family, graph, cost, Curve([r.k / unit for r in rs], recovery_percent([r.z for r in rs]))))
    end
    return out
end

"""VoR curve of every family for the designated-path costs (one instance each)."""
function bottleneck_curves(rows)::Dict{String,Curve}
    out = Dict{String,Curve}()
    for family in FAMILY_ORDER
        rs = sort([r for r in rows if r.family == family]; by = r -> r.k)
        isempty(rs) && continue
        check_ks(rs)
        unit = budget_unit(family, rs[1].ell)
        out[family] = Curve([r.k / unit for r in rs], recovery_percent([r.z for r in rs]))
    end
    return out
end

# --- continuous-budget and discrete-budget experiments ---------------------------------

"""Fractions of the forcing budget that are shown."""
const BUDGET_FRACTIONS = (0.2, 0.4, 0.6, 0.8)

function continuous_instance_curves(rows, frac::Real)::Vector{InstanceCurve}
    out = InstanceCurve[]
    for ((family, graph, cost), rs) in by_instance(rows)
        sel = [r for r in rs if abs(r.gamma_frac - frac) < 1e-9]
        isempty(sel) && error("no rows for fraction $frac ($family $graph $cost)")
        check_ks(sel)
        unit = budget_unit(family, sel[1].ell)
        push!(out, InstanceCurve(family, graph, cost, Curve([r.k / unit for r in sel], recovery_percent([r.z for r in sel]))))
    end
    return out
end

"""
Discrete budget `Δ = frac · Δ_force`; a fractional value is the linear
interpolation of the values at the two neighboring integer budgets.
"""
function discrete_instance_curves(rows, frac::Real)::Vector{InstanceCurve}
    out = InstanceCurve[]
    for ((family, graph, cost), rs) in by_instance(rows)
        dforce = rs[1].delta_force
        target = round(frac * dforce; digits = 9)
        lo, hi = floor(Int, target), ceil(Int, target)
        w = target - lo
        zs(delta) = begin
            sel = [r for r in rs if r.delta == delta]
            isempty(sel) && error("no rows for Δ=$delta ($family $graph $cost)")
            check_ks(sel)
            [r.z for r in sel]
        end
        zlo = zs(lo)
        z = lo == hi ? zlo : (1 - w) .* zlo .+ w .* zs(hi)
        unit = budget_unit(family, rs[1].ell)
        push!(out, InstanceCurve(family, graph, cost, Curve([k / unit for k in 0:length(z) - 1], recovery_percent(z))))
    end
    return out
end

"""Average curves of all four families, on the grid `j/12`, for each fraction."""
function budget_family_curves(make_instances, rows; grid::Bool = true)
    out = Dict{Float64,Dict{String,Curve}}()
    for frac in BUDGET_FRACTIONS
        ics = make_instances(rows, frac)
        out[frac] = Dict(f => family_curve(filter(ic -> ic.family == f, ics); grid = grid) for f in FAMILY_ORDER if any(ic -> ic.family == f, ics))
    end
    return out
end

# --- instance tables -------------------------------------------------------------------

struct SizeRange
    lo::Int
    hi::Int
end
SizeRange(v::Integer) = SizeRange(v, v)
Base.union(a::SizeRange, b::SizeRange) = SizeRange(min(a.lo, b.lo), max(a.hi, b.hi))

"""Per-family ranges of `n`, `m`, number of paths (if present), `ℓ` and the maximum number of arcs."""
function instance_sizes(rows)
    out = Dict{String,Dict{Symbol,SizeRange}}()
    for r in rows
        r.k == 0 || continue
        d = get!(out, r.family, Dict{Symbol,SizeRange}())
        fields = hasproperty(r, :npaths) ? (:n, :m, :npaths, :ell, :longest) : (:n, :m, :ell, :longest)
        for f in fields
            v = SizeRange(getproperty(r, f))
            d[f] = haskey(d, f) ? union(d[f], v) : v
        end
    end
    return out
end
