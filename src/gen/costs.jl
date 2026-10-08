"""
    with_first_stage(costs, mode) -> Costs

Copy of `costs` with first-stage vector set by `mode`:

- `:zero`: ``C = 0`` (isolates robustness).
- `:c_hat`: ``C = \\hat{c}`` (pay now and after recovery).
"""
function with_first_stage(costs::Costs, mode::Symbol)::Costs
    m = length(costs.C)
    if mode === :zero
        return costs_new(zeros(Float64, m), costs.c_hat, costs.d)
    elseif mode === :c_hat
        return costs_new(costs.c_hat, costs.c_hat, costs.d)
    end
    throw(ArgumentError("unknown first-stage mode $mode"))
end

function _rand_c_hat(m::Int, rng::AbstractRNG, c_hat_lo::Float64, c_hat_hi::Float64)::Vector{Float64}
    c_hat_hi < c_hat_lo && throw(ArgumentError("c_hat_hi < c_hat_lo"))
    out = Vector{Float64}(undef, m)
    span = c_hat_hi - c_hat_lo
    @inbounds for a in 1:m
        out[a] = c_hat_lo + span * rand(rng)
    end
    return out
end

"""
    overlay_uniform(m; rng, c_hat_lo, c_hat_hi, d_lo, d_hi) -> Costs
    overlay_uniform(g; ...) -> Costs

Independent uniform nominal costs in `[c_hat_lo, c_hat_hi]` and deviations in
`[d_lo, d_hi]`. First-stage cost is zero; apply [`with_first_stage`](@ref)
afterwards.
"""
function overlay_uniform(
    m::Integer;
    rng::AbstractRNG = Random.default_rng(),
    c_hat_lo::Real = 1.0,
    c_hat_hi::Real = 10.0,
    d_lo::Real = 0.0,
    d_hi::Real = 10.0,
)::Costs
    mm = Int(m)
    c_hat = _rand_c_hat(mm, rng, Float64(c_hat_lo), Float64(c_hat_hi))
    d = _rand_c_hat(mm, rng, Float64(d_lo), Float64(d_hi))
    return costs_new(zeros(Float64, mm), c_hat, d)
end

overlay_uniform(g::Graph; kwargs...) = overlay_uniform(Int(g.m); kwargs...)

"""
    overlay_alpha(m; alpha, rng) -> Costs

Cost structure with ``\\hat{c}_a = α(\\hat{c}_a + d_a)`` for every arc,
`alpha ∈ (0, 1]`. First-stage cost is zero.
"""
function overlay_alpha(
    m::Integer;
    alpha::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
    c_hat_lo::Real = 2.0,
    c_hat_hi::Real = 10.0,
)::Costs
    a = Float64(alpha)
    (0.0 < a <= 1.0) || throw(ArgumentError("alpha must lie in (0, 1]"))
    mm = Int(m)
    c_hat = _rand_c_hat(mm, rng, Float64(c_hat_lo), Float64(c_hat_hi))
    scale = (1.0 - a) / a
    d = Vector{Float64}(undef, mm)
    @inbounds for i in 1:mm
        d[i] = c_hat[i] * scale
    end
    return costs_new(zeros(Float64, mm), c_hat, d)
end

overlay_alpha(g::Graph; kwargs...) = overlay_alpha(Int(g.m); kwargs...)

"""
    overlay_zero_nominal(m; n_zero, rng) -> Costs

`n_zero` arcs have ``\\hat{c} = 0`` and large deviation (the ``α``-condition
fails). Remaining arcs have positive nominal cost. First-stage cost is zero.
"""
function overlay_zero_nominal(
    m::Integer;
    n_zero::Integer = 1,
    d_zero::Real = 10.0,
    rng::AbstractRNG = Random.default_rng(),
    c_hat_lo::Real = 2.0,
    c_hat_hi::Real = 10.0,
    d_rest::Real = 1.0,
)::Costs
    mm = Int(m)
    nz = Int(n_zero)
    nz < 0 && throw(ArgumentError("n_zero < 0"))
    nz > mm && throw(ArgumentError("n_zero > m"))
    c_hat = _rand_c_hat(mm, rng, Float64(c_hat_lo), Float64(c_hat_hi))
    d = fill(Float64(d_rest), mm)
    if nz > 0
        hot = randperm(rng, mm)[1:nz]
        @inbounds for a in hot
            c_hat[a] = 0.0
            d[a] = Float64(d_zero)
        end
    end
    return costs_new(zeros(Float64, mm), c_hat, d)
end

overlay_zero_nominal(g::Graph; kwargs...) = overlay_zero_nominal(Int(g.m); kwargs...)

"""
    overlay_designated_bottleneck(g, hot_arcs; c_hat_hot=1, d_hot=100, c_hat_cold=8, d_cold=1) -> Costs

Designated-path costs of the "bottleneck" experiment. The arcs in `hot_arcs`
(a chosen ``s``–``t`` path) are nominally cheap but have a large deviation;
every other arc is nominally expensive with a small deviation, so the path
that is best nominally and the path that is best in the worst case are
arc-disjoint. The first-stage cost equals the nominal cost.
"""
function overlay_designated_bottleneck(
    g::Graph,
    hot_arcs::AbstractVector{<:Integer};
    c_hat_hot::Real = 1.0,
    d_hot::Real = 100.0,
    c_hat_cold::Real = 8.0,
    d_cold::Real = 1.0,
)::Costs
    m = Int(g.m)
    c_hat = fill(Float64(c_hat_cold), m)
    d = fill(Float64(d_cold), m)
    for a in hot_arcs
        (1 <= a <= m) || throw(ArgumentError("arc $a is outside 1:$m"))
        c_hat[a] = Float64(c_hat_hot)
        d[a] = Float64(d_hot)
    end
    return costs_new(c_hat, c_hat, d)
end

"""
    overlay_asp_detours(H; rng, c_hat_lo=1, c_hat_hi=10, d_lo=0, d_hi=10) -> Costs

Independent costs for the series of `H` gadgets of [`gen_asp_detours`](@ref),
drawn *per alternative*. Drawing per arc would give the three-arc detour three
independent positive draws against one for the direct arc, so the direct arc
would win in every gadget. Instead each gadget draws one pair
``(\\hat{c}, d)`` for the direct arc and one for the detour as a whole, and
the detour pair is divided equally among its three arcs. The first-stage cost
equals the nominal cost.
"""
function overlay_asp_detours(
    H::Integer;
    rng::AbstractRNG = Random.default_rng(),
    c_hat_lo::Real = 1.0,
    c_hat_hi::Real = 10.0,
    d_lo::Real = 0.0,
    d_hi::Real = 10.0,
)::Costs
    hh = Int(H)
    hh >= 1 || throw(ArgumentError("H < 1"))
    cspan = Float64(c_hat_hi) - Float64(c_hat_lo)
    dspan = Float64(d_hi) - Float64(d_lo)
    c_hat = zeros(4 * hh)
    d = zeros(4 * hh)
    for gix in 0:(hh - 1)
        base = 4 * gix
        c_hat[base + 1] = Float64(c_hat_lo) + cspan * rand(rng)
        d[base + 1] = Float64(d_lo) + dspan * rand(rng)
        bypass_c = Float64(c_hat_lo) + cspan * rand(rng)
        bypass_d = Float64(d_lo) + dspan * rand(rng)
        for j in 2:4
            c_hat[base + j] = bypass_c / 3
            d[base + j] = bypass_d / 3
        end
    end
    return costs_new(c_hat, c_hat, d)
end
