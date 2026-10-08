# Driver loop for the per-instance experiments (continuous-budget and discrete-budget experiments).
#
# Every (digraph, cost draw) is solved on its own and written to
# `<root>/<experiment>/<family>/g<graph>_c<cost>.csv`. An existing file is skipped, so
# an interrupted run resumes, and disjoint slices (`--family`, `--from`, `--to`)
# can run in separate processes.

const INSTANCE_META_HEADER = ["family", "graph", "cost", "graph_seed", "cost_seed"]

instance_path(root, experiment::AbstractString, inst::Instance) =
    joinpath(root, experiment, inst.family, "g$(inst.graph)_c$(inst.cost).csv")

"""
    run_instances(experiment, name, header, solve_instance; opts, root)

Solve the instances selected by the command-line options `opts`:
`--scale=quick|full`, `--family=a,b` (default: all), `--from=i --to=j`
(cost draws of the fixed families, graphs of the random DAGs), `--force`
(recompute existing files). `solve_instance(inst)` returns tuples that follow
`header` after the five instance columns.
"""
function run_instances(experiment::Symbol, name::AbstractString, header, solve_instance::Function; opts, root::AbstractString = run_root())
    scale = opt_scale(opts)
    families = String.(split(opt(opts, "family", join(FAMILIES, ",")), ','))
    all(f -> f in FAMILIES, families) || error("--family must be a subset of $(join(FAMILIES, ','))")
    force = haskey(opts, "force")
    written = String[]
    for family in families
        costs, graphs = scale_ranges(family, scale)
        rng = family == "random_dag" ? graphs : costs
        lo = opt_int(opts, "from", first(rng))
        hi = opt_int(opts, "to", last(rng))
        keep = intersect(rng, lo:hi)
        isempty(keep) && continue
        sel_costs, sel_graphs = family == "random_dag" ? (costs, keep) : (keep, graphs)
        for inst in family_instances(experiment_for_scale(experiment, scale), family; costs = sel_costs, graphs = sel_graphs)
            dest = instance_path(root, name, inst)
            if isfile(dest) && !force
                log_line("skip $dest")
                continue
            end
            t0 = time()
            log_line("$name $family graph=$(inst.graph) cost=$(inst.cost) seed=$(inst.cost_seed)")
            meta = (inst.family, inst.graph, inst.cost, inst.graph_seed, inst.cost_seed)
            rows = [(meta..., r...) for r in solve_instance(inst)]
            write_csv(dest, vcat(INSTANCE_META_HEADER, header), rows)
            log_line("  wrote $(length(rows)) rows in $(round(time() - t0; digits = 1)) s")
            push!(written, dest)
        end
    end
    return written
end
