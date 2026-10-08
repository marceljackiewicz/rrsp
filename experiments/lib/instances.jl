# The thesis instances: digraph families, sizes, and seeds for the interval, continuous-budget and discrete-budget experiments.
#
# One fixed digraph per family (layered, grid, ASP) with 30 independent cost
# draws, and ten random DAGs with three cost draws each. Seeds are fixed so
# that every run reproduces the same instances.

const GRAPH_SEED = 20260910
const COST_SEED_BASE = Dict("layered" => 12, "grid" => 13, "asp" => 14)
const FAMILIES = ("layered", "random_dag", "grid", "asp")
const N_COST_FIXED = 30
const N_GRAPHS_RANDOM = 10
const N_COST_RANDOM = 3
const SKIP_P = 0.5

"""Sizes per experiment: `H` arcs per s–t path (ℓ for layered/grid/ASP), layer width `W`, grid side."""
const EXPERIMENT_SPECS = Dict(
    :interval => (H = 12, W = 4, grid = 7),
    :cont => (H = 12, W = 2, grid = 7),
    :disc => (H = 8, W = 2, grid = 5),
    # Reduced digraphs for the `quick` scale (a few minutes in total).
    :cont_quick => (H = 6, W = 2, grid = 4),
    :disc_quick => (H = 4, W = 2, grid = 3),
)

"""The size specification used for `experiment` at `scale`."""
experiment_for_scale(experiment::Symbol, scale::AbstractString) =
    scale == "quick" && haskey(EXPERIMENT_SPECS, Symbol(experiment, :_quick)) ? Symbol(experiment, :_quick) : experiment

struct Instance
    experiment::Symbol
    family::String
    graph::Int          # graph index (0 for the fixed families)
    cost::Int           # cost-draw index
    graph_seed::Int
    cost_seed::Int
    g::Graph
    s::Int
    t::Int
end

function fixed_graph(experiment::Symbol, family::AbstractString)::Graph
    sp = EXPERIMENT_SPECS[experiment]
    family == "layered" && return gen_wide_layered(sp.H, sp.W)
    family == "grid" && return gen_grid(sp.grid, sp.grid)
    family == "asp" && return gen_asp_detours(sp.H)
    error("not a fixed family: $family")
end

random_dag_graph(experiment::Symbol, i::Integer)::Graph =
    gen_layered_skips(EXPERIMENT_SPECS[experiment].H, EXPERIMENT_SPECS[experiment].W; p = SKIP_P, rng = MersenneTwister(GRAPH_SEED + i))

"""Instances of `experiment` and `family` with cost indices in the given ranges.
`graphs` indexes the random DAGs and is ignored for the fixed families."""
function family_instances(experiment::Symbol, family::AbstractString; costs::AbstractVector{<:Integer}, graphs::AbstractVector{<:Integer} = 0:0)
    out = Instance[]
    if family == "random_dag"
        for i in graphs
            g = random_dag_graph(experiment, i)
            for r in costs
                push!(out, Instance(experiment, family, i, r, GRAPH_SEED + i, 15 + 100 * i + r, g, 1, Int(g.n)))
            end
        end
    else
        g = fixed_graph(experiment, family)
        for r in costs
            push!(out, Instance(experiment, family, 0, r, 0, COST_SEED_BASE[family] + 1000 * r, g, 1, Int(g.n)))
        end
    end
    return out
end

"""Index ranges `(costs, graphs)` for a family at the given scale."""
function scale_ranges(family::AbstractString, scale::AbstractString)
    if family == "random_dag"
        return scale == "full" ? (0:(N_COST_RANDOM - 1), 0:(N_GRAPHS_RANDOM - 1)) : (0:0, 0:1)
    end
    return scale == "full" ? (0:(N_COST_FIXED - 1), 0:0) : (0:1, 0:0)
end

function scale_instances(experiment::Symbol, family::AbstractString, scale::AbstractString)
    costs, graphs = scale_ranges(family, scale)
    return family_instances(experiment, family; costs = costs, graphs = graphs)
end

"""Independent second-stage costs of an instance, drawn per arc, or per
alternative on the ASP digraph. `first_stage` is `:c_hat` (C = ĉ, interval experiment) or
`:zero` (C = 0, continuous-budget and discrete-budget experiments)."""
function instance_costs(inst::Instance; first_stage::Symbol)::Costs
    rng = MersenneTwister(inst.cost_seed)
    base = if inst.family == "asp"
        overlay_asp_detours(EXPERIMENT_SPECS[inst.experiment].H; rng = rng)
    else
        with_first_stage(overlay_uniform(inst.g; rng = rng), :c_hat)
    end
    return first_stage === :c_hat ? base : with_first_stage(base, first_stage)
end

# --- arc counts ---------------------------------------------------------------

function longest_arc_count(g::Graph, s::Integer, t::Integer)::Int
    is_dag(g) || error("longest_arc_count needs a DAG")
    n = Int(g.n)
    best = fill(-1, n)
    best[s] = 0
    indeg = [Int(in_degree(g, v)) for v in 1:n]
    queue = [v for v in 1:n if indeg[v] == 0]
    head = 1
    while head <= length(queue)
        v = queue[head]
        head += 1
        for a in outgoing(g, v)
            w = Int(g.head[a])
            best[v] >= 0 && best[v] + 1 > best[w] && (best[w] = best[v] + 1)
            indeg[w] -= 1
            indeg[w] == 0 && push!(queue, w)
        end
    end
    return best[t]
end

function shortest_arc_count(g::Graph, s::Integer, t::Integer)::Int
    m = Int(g.m)
    w = ones(Float64, m)
    sol = solve_sp(network_new(g, s, t, costs_new(zeros(m), w, zeros(m))), w, Solver(; method = METHOD_COMB))
    sol.status == ST_OK || error("no s-t path")
    return Int(round(sol.z))
end

# --- designated bottleneck path (interval experiment) -------------------------------------

function arc_between(g::Graph, u::Int, v::Int)::Int
    for a in outgoing(g, Int32(u))
        Int(g.head[a]) == v && return Int(a)
    end
    error("no arc $u → $v")
end

"""Arcs of the designated cheap-hot path: the first vertex of each layer
(layered, random DAG), the first row then the last column (grid), the direct
arc of every gadget (ASP)."""
function designated_arcs(experiment::Symbol, family::AbstractString, g::Graph)::Vector{Int}
    sp = EXPERIMENT_SPECS[experiment]
    verts = if family == "asp"
        collect(1:3:Int(g.n))
    elseif family == "grid"
        nc = sp.grid
        v = collect(1:nc)
        for r in 2:nc
            push!(v, (r - 1) * nc + nc)
        end
        v
    elseif family == "layered" || family == "random_dag"
        v = [1]
        for layer in 1:(sp.H - 1)
            push!(v, first(layer_vertex_ids(sp.H, sp.W, layer)))
        end
        push!(v, Int(g.n))
        v
    else
        error("no designated path for $family")
    end
    return [arc_between(g, verts[i], verts[i + 1]) for i in 1:(length(verts) - 1)]
end
