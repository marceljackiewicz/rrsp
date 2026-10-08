using Test
using Random
using HiGHS
import JuMP
using Rrsp

include("setup.jl")

@testset "Rrsp" begin
    @testset "core/graph" begin
        include("core/graph.jl")
    end
    @testset "core/path" begin
        include("core/path.jl")
    end
    @testset "core/costs_network" begin
        include("core/costs_network.jl")
    end
    @testset "core/params" begin
        include("core/params.jl")
    end
    @testset "core/io" begin
        include("core/io.jl")
    end
    @testset "core/asp" begin
        include("core/asp.jl")
    end
    @testset "eval/eval_costs" begin
        include("eval/eval_costs.jl")
    end
    @testset "eval/sp" begin
        include("eval/sp.jl")
    end
    @testset "solve/rob" begin
        include("solve/rob.jl")
    end
    @testset "eval/inc" begin
        include("eval/inc.jl")
    end
    @testset "eval/adv" begin
        include("eval/adv.jl")
    end
    @testset "solve/rec" begin
        include("solve/rec.jl")
    end
    @testset "solve/rrsp" begin
        include("solve/rrsp.jl")
    end
    @testset "solve/audit_regressions" begin
        include("solve/audit_regressions.jl")
    end
    @testset "solve/random_oracle" begin
        include("solve/random_oracle.jl")
    end
    @testset "solve/rrsp_enum" begin
        include("solve/rrsp_enum.jl")
    end
    @testset "approx/alpha" begin
        include("approx/alpha.jl")
    end
    @testset "identities" begin
        include("identities.jl")
    end
    @testset "regression" begin
        include("regression.jl")
    end
    @testset "gen" begin
        include("gen.jl")
    end
    @testset "experiments" begin
        include("experiments.jl")
    end
    if get(ENV, "RRSP_TEST_EXPERIMENTS", "0") == "1"
        include("experiments_regression.jl")
    else
        @info "Skipping the experiment pipeline regression test (set RRSP_TEST_EXPERIMENTS=1 to run it)"
    end
end
