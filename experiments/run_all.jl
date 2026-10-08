# Every experiment, then every figure and table, in one command.
#
#   julia --project=experiments experiments/run_all.jl [--scale=quick|full] [--pdf]
#
# Raw results go to experiments/out/data and the generated TeX/SVG files to
# experiments/out/artifacts, so the tracked data/ stays as it is (see
# promote_data.jl). `quick` (default) takes about a minute on reduced
# digraphs; `full` is the thesis scale and takes much longer (the discrete-budget experiment is by
# far the slowest). Interrupted runs resume: finished instances are skipped.
# Slice a run over several processes with --family, --from, --to of cont.jl
# and disc.jl, then call make_artifacts.jl.

include(joinpath(@__DIR__, "interval.jl"))
include(joinpath(@__DIR__, "cont.jl"))
include(joinpath(@__DIR__, "disc.jl"))
include(joinpath(@__DIR__, "make_artifacts.jl"))

function run_all(opts)
    scale = opt_scale(opts)
    log_line("== interval ($scale)")
    run_interval(scale)
    log_line("== continuous budget ($scale)")
    run_cont(opts)
    log_line("== discrete budget ($scale)")
    run_disc(opts)
    log_line("== figures and tables")
    out = opt(opts, "out", joinpath(EXPERIMENTS_DIR, "out", "artifacts"))
    files = make_artifacts(; data = run_root(), out = out)
    haskey(opts, "pdf") && make_preview(out, files)
    return files
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all(parse_options())
end
