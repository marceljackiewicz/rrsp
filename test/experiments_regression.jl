# Opt-in regression test of the experiment pipeline:
#
#   RRSP_TEST_EXPERIMENTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
#
# Runs the scripts in `experiments/` as subprocesses (they have their own
# environment) and checks that
#   1. the figures and tables build from the tracked data,
#   2. a quick run of the interval experiment reproduces rows of the tracked data,
#   3. quick runs of the continuous-budget and discrete-budget experiments finish and write every instance.

const REPO = normpath(joinpath(@__DIR__, ".."))
const EXP = joinpath(REPO, "experiments")

# Pkg.test restricts the load path of its own process; the experiment scripts
# need their normal environment.
function exp_run(cmd::Cmd; env = Dict{String,String}(), quiet = true)
    cmd = addenv(Cmd(cmd; dir = REPO), "JULIA_LOAD_PATH" => nothing, "JULIA_PROJECT" => nothing, env...)
    run(quiet ? pipeline(cmd; stdout = devnull) : cmd)
end

julia_script(script, args...; env = Dict{String,String}()) =
    exp_run(`$(Base.julia_cmd()) --project=$EXP $(joinpath(EXP, script)) $args`; env = env)

@testset "experiments pipeline" begin
    # instantiate once; needs network on a fresh machine
    exp_run(`$(Base.julia_cmd()) --project=$EXP -e "using Pkg; Pkg.instantiate()"`; quiet = false)

    mktempdir() do tmp
        @testset "figures and tables build from the tracked data" begin
            out = joinpath(tmp, "artifacts")
            julia_script("make_artifacts.jl", "--out=$out", "--no-svg")
            figures = filter(endswith(".tex"), readdir(joinpath(out, "figures")))
            tables = filter(endswith(".tex"), readdir(joinpath(out, "tables")))
            @test length(figures) == 4
            @test length(tables) == 7
            @test all(f -> filesize(joinpath(out, "figures", f)) > 0, figures)
            @test all(f -> filesize(joinpath(out, "tables", f)) > 0, tables)
        end

        @testset "interval quick-run rows are rows of the tracked data" begin
            run_dir = joinpath(tmp, "run_interval")
            julia_script("interval.jl", "--scale=quick"; env = Dict("RRSP_RUN_DIR" => run_dir))
            for name in ("replicates", "bottleneck")
                tracked = Set(readlines(joinpath(EXP, "data", "interval", name * ".csv")))
                fresh = readlines(joinpath(run_dir, "interval", name * ".csv"))
                @test length(fresh) > 1
                @test all(in(tracked), fresh)
            end
        end

        @testset "continuous-budget and discrete-budget quick runs" begin
            run_dir = joinpath(tmp, "run_cont_disc")
            for experiment in ("cont.jl", "disc.jl")
                julia_script(experiment, "--scale=quick"; env = Dict("RRSP_RUN_DIR" => run_dir))
            end
            for experiment in ("cont", "disc")
                files = String[]
                for (root, _, names) in walkdir(joinpath(run_dir, experiment))
                    append!(files, joinpath.(root, filter(endswith(".csv"), names)))
                end
                @test length(files) == 8       # 2 draws x 3 fixed families + 2 random DAGs
                @test all(f -> countlines(f) > 1, files)
            end
        end
    end
end
