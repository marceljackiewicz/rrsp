@testset "enum integer values" begin
    @test Int(NB_INCLUSION) == 0
    @test Int(NB_EXCLUSION) == 1
    @test Int(NB_SYMDIFF) == 2
    @test Int(U_NOMINAL) == 0
    @test Int(U_INTERVAL) == 1
    @test Int(U_CONT_BUDGET) == 2
    @test Int(U_DISC_BUDGET) == 3
    @test Int(METHOD_AUTO) == 0
    @test Int(METHOD_COMB) == 1
    @test Int(METHOD_MIP) == 2
    @test Int(ST_OK) == 0
    @test Int(ST_INFEASIBLE) == 1
    @test Int(ST_TIME_LIMIT) == 2
    @test Int(ST_UNBOUNDED) == 3
    @test Int(ST_NOT_IMPL) == 4
    @test Int(ST_ERROR) == 5
end

@testset "params_nominal" begin
    p = params_nominal()
    @test p.uncertainty == U_NOMINAL
    @test p.neighborhood == NB_INCLUSION
    @test p.k == 0
    @test p.gamma == 0.0
    @test p.delta == 0
end

@testset "params_rob forces k = 0" begin
    p = params_rob(U_INTERVAL)
    @test p.k == 0
    @test p.uncertainty == U_INTERVAL
    p2 = params_rob(U_CONT_BUDGET; gamma = 3.5)
    @test p2.k == 0
    @test p2.gamma == 3.5
    @test p2.delta == 0
    p3 = params_rob(U_DISC_BUDGET; delta = 4)
    @test p3.k == 0
    @test p3.delta == 4
    @test p3.gamma == 0.0
end

@testset "params_rrsp records neighborhood and k" begin
    p = params_rrsp(U_INTERVAL, NB_EXCLUSION, 3)
    @test p.uncertainty == U_INTERVAL
    @test p.neighborhood == NB_EXCLUSION
    @test p.k == 3
    p2 = params_rrsp(U_CONT_BUDGET, NB_SYMDIFF, 2; gamma = 10.0)
    @test p2.gamma == 10.0
    @test p2.neighborhood == NB_SYMDIFF
    p3 = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 1; delta = 2)
    @test p3.delta == 2
end

@testset "params_rrsp rejects negative k" begin
    @test_throws ArgumentError params_rrsp(U_INTERVAL, NB_INCLUSION, -1)
end

@testset "interval params may still carry unused budgets" begin
    p = params_rrsp(U_INTERVAL, NB_INCLUSION, 1; gamma = 9.0, delta = 3)
    @test p.uncertainty == U_INTERVAL
    @test p.gamma == 9.0
    @test p.delta == 3
    @test p.k == 1
end

@testset "solution_empty" begin
    sol = solution_empty(4)
    @test sol.status == ST_INFEASIBLE
    @test sol.z == Inf
    @test sol.z_first == Inf
    @test sol.z_second == Inf
    @test isempty(sol.first.seq)
    @test isempty(sol.second.seq)
    @test length(sol.first.chi) == 4
    @test length(sol.second.chi) == 4
    @test all(==(0x00), sol.first.chi)
    @test sol.mip_gap == 0.0
    @test sol.mip_nodes == 0
    @test sol.n_binaries == 0
    @test sol.n_constraints == 0
end

@testset "Solver defaults" begin
    slv = Solver()
    @test slv.optimizer === nothing
    @test slv.method == METHOD_AUTO
    @test slv.time_limit == Inf
    @test slv.mip_gap == 1e-6
    @test slv.threads == 0
    @test slv.silent == true
end

@testset "Solver keyword overrides" begin
    slv = Solver(; method = METHOD_COMB, time_limit = 2.5, silent = false, threads = 3)
    @test slv.method == METHOD_COMB
    @test slv.time_limit == 2.5
    @test slv.silent == false
    @test slv.threads == 3
    @test slv.optimizer === nothing
end
