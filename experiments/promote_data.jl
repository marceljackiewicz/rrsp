# Copy a finished run into the tracked data directory.
#
#   julia --project=experiments experiments/promote_data.jl [--from=DIR] [--to=DIR]
#
# Merges the per-instance files of the continuous-budget and discrete-budget experiments into one CSV each
# (cont/all.csv, disc/all.csv, rows in a fixed order) and copies the two
# interval files. Defaults: --from experiments/out/data, --to experiments/data.
# Refuses to mix scales: every family needs the full set of instances.

@isdefined(EXPERIMENTS_DIR) || include(joinpath(@__DIR__, "lib", "common.jl"))
@isdefined(EXPERIMENT_SPECS) || include(joinpath(@__DIR__, "lib", "instances.jl"))

function promote(from::AbstractString, to::AbstractString)
    for name in ("replicates", "bottleneck")
        src = joinpath(from, "interval", name * ".csv")
        isfile(src) || error("missing $src")
        mkpath(joinpath(to, "interval"))
        cp(src, joinpath(to, "interval", name * ".csv"); force = true)
    end
    order = Dict(f => i for (i, f) in enumerate(FAMILIES))
    for experiment in ("cont", "disc")
        files = String[]
        for (root, _, names) in walkdir(joinpath(from, experiment))
            append!(files, joinpath.(root, filter(endswith(".csv"), names)))
        end
        isempty(files) && error("no result files in $(joinpath(from, experiment))")
        header = ""
        rows = NamedTuple[]
        for f in files
            header = readline(f)
            append!(rows, read_csv(f))
        end
        cols = Tuple(Symbol.(split(header, ',')))
        for fam in FAMILIES
            costs, graphs = scale_ranges(fam, "full")
            expected = fam == "random_dag" ? length(costs) * length(graphs) : length(costs)
            got = length(unique((r.graph, r.cost) for r in rows if r.family == fam))
            got == expected || error("$experiment/$fam has $got instances, the thesis scale has $expected; run --scale=full first")
        end
        sort!(rows; by = r -> (order[r.family], r.graph, r.cost, get(r, :gamma_frac, 0.0), get(r, :delta, 0), r.k))
        dest = joinpath(to, experiment, "all.csv")
        write_csv(dest, collect(String.(cols)), (Tuple(getproperty(r, c) for c in cols) for r in rows))
        log_line("wrote $dest ($(length(rows)) rows)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_options()
    promote(opt(opts, "from", run_root()), opt(opts, "to", data_root()))
end
