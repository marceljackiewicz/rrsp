# Interval experiment: value of recovery under interval uncertainty (Rec, exact).
#
#   julia --project=experiments experiments/interval.jl [--scale=quick|full]
#
# Writes (below `run_root()`, default experiments/out/data):
#   interval/replicates.csv  independent costs, C = ĉ: 30 draws on each fixed digraph,
#                          3 draws on each of 10 random DAGs
#   interval/bottleneck.csv  designated-path costs on one draw of every family

@isdefined(EXPERIMENTS_DIR) || include(joinpath(@__DIR__, "lib", "common.jl"))
@isdefined(EXPERIMENT_SPECS) || include(joinpath(@__DIR__, "lib", "instances.jl"))

const REPLICATE_HEADER = ["family", "graph", "cost", "graph_seed", "cost_seed", "n", "m", "ell", "longest", "k", "z"]
const BOTTLENECK_HEADER = ["family", "graph", "graph_seed", "n", "m", "ell", "longest", "k", "z"]

function interval_value(net::Network, k::Int, solver::Solver)::Float64
    sol = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, k), solver)
    sol.status == ST_OK || error("Rec failed: $(sol.status) (k = $k)")
    return sol.z
end

# ASP detours reach k = 3H (all gadgets replaced); other families saturate at ℓ.
budgets(family::AbstractString, H::Int, ell::Int) = family == "asp" ? (0:(3 * H)) : (0:ell)

function run_replicates(scale::AbstractString; root::AbstractString = run_root())
    solver = Solver(; method = METHOD_COMB)
    H = EXPERIMENT_SPECS[:interval].H
    rows = Tuple[]
    for family in FAMILIES
        for inst in scale_instances(:interval, family, scale)
            g = inst.g
            net = network_new(g, inst.s, inst.t, instance_costs(inst; first_stage = :c_hat))
            ell = shortest_arc_count(g, inst.s, inst.t)
            longest = longest_arc_count(g, inst.s, inst.t)
            for k in budgets(family, H, ell)
                push!(rows, (family, inst.graph, inst.cost, inst.graph_seed, inst.cost_seed,
                             Int(g.n), Int(g.m), ell, longest, k, interval_value(net, k, solver)))
            end
        end
        log_line("replicates: $family done")
    end
    return write_csv(joinpath(root, "interval", "replicates.csv"), REPLICATE_HEADER, rows)
end

function run_bottleneck(; root::AbstractString = run_root())
    solver = Solver(; method = METHOD_COMB)
    H = EXPERIMENT_SPECS[:interval].H
    rows = Tuple[]
    for family in FAMILIES
        inst = first(family_instances(:interval, family; costs = 0:0, graphs = 0:0))
        g = inst.g
        hot = designated_arcs(:interval, family, g)
        length(hot) == H || error("$family: designated path has $(length(hot)) arcs, expected $H")
        net = network_new(g, inst.s, inst.t, overlay_designated_bottleneck(g, hot))
        ell = shortest_arc_count(g, inst.s, inst.t)
        longest = longest_arc_count(g, inst.s, inst.t)
        for k in budgets(family, H, ell)
            push!(rows, (family, inst.graph, inst.graph_seed, Int(g.n), Int(g.m), ell, longest, k,
                         interval_value(net, k, solver)))
        end
        log_line("bottleneck: $family done")
    end
    return write_csv(joinpath(root, "interval", "bottleneck.csv"), BOTTLENECK_HEADER, rows)
end

function run_interval(scale::AbstractString; root::AbstractString = run_root())
    return (run_replicates(scale; root = root), run_bottleneck(; root = root))
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_options()
    run_interval(opt_scale(opts))
end
