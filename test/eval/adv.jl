function adv_diamond()
    g = fixture_two_paths()
    costs = costs_new(
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 10.0, 10.0],
        [4.0, 4.0, 1.0, 1.0],
    )
    return network_new(g, 1, 4, costs)
end

@testset "interval ADV equals INC under c_hat+d" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    w = cmax_of(net)
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        for k in (0, 1, 2)
            params = params_rrsp(U_INTERVAL, nb, k)
            adv = solve_adv(net, params, x, COMB_SOLVER)
            inc = solve_inc(net, params, x, w, COMB_SOLVER)
            @test adv.status == inc.status == ST_OK
            @test adv.z == inc.z
            @test eval_worstcase(net, params, x, COMB_SOLVER) == adv.z
        end
    end
end

@testset "continuous ADV k = 0 dumps Gamma on x" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    for gamma in (0.0, 3.0, 8.0, 100.0)
        params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 0; gamma = gamma)
        adv = solve_adv(net, params, x, COMB_SOLVER)
        @test adv.status == ST_OK
        @test adv.second.seq == x.seq
        @test adv.z == dump_gamma_on_path(x, net.costs.c_hat, net.costs.d, gamma)
        @test adv.method_used == METHOD_COMB
    end
end

@testset "continuous ADV Gamma = 0 equals INC under c_hat" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 2; gamma = 0.0)
    adv = solve_adv(net, params, x, COMB_SOLVER)
    inc = solve_inc(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 2), x, net.costs.c_hat, COMB_SOLVER)
    @test adv.status == ST_OK
    @test adv.z == inc.z
end

@testset "continuous ADV large Gamma equals INC under c_hat+d" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    gamma = sum(net.costs.d)
    params = params_rrsp(U_CONT_BUDGET, NB_SYMDIFF, 2; gamma = gamma)
    adv = solve_adv(net, params, x, COMB_SOLVER)
    inc = solve_inc(net, params_rrsp(U_INTERVAL, NB_SYMDIFF, 2), x, cmax_of(net), COMB_SOLVER)
    @test adv.status == ST_OK
    @test adv.z == inc.z
end

@testset "continuous inclusion ADV LP vs enumeration" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    for k in 1:2
        for gamma in (2.0, 6.0)
            params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, k; gamma = gamma)
            z_or = oracle_adv_cont_z(
                net.graph, net.s, net.t, x, net.costs.c_hat, net.costs.d, NB_INCLUSION, k, gamma
            )
            sol = solve_adv(net, params, x, TEST_SOLVER)
            @test sol.status == ST_OK
            @test approxz(sol.z, z_or)
            @test sol.method_used == METHOD_MIP
        end
    end
end

@testset "continuous ADV on DAG exclusion and symdiff vs enumeration" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    for (nb, k) in ((NB_EXCLUSION, 1), (NB_SYMDIFF, 2))
        params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = 5.0)
        z_or = oracle_adv_cont_z(
            net.graph, net.s, net.t, x, net.costs.c_hat, net.costs.d, nb, k, 5.0
        )
        sol = solve_adv(net, params, x, TEST_SOLVER)
        @test sol.status == ST_OK
        @test approxz(sol.z, z_or)
    end
end

@testset "continuous inclusion ADV on a cyclic graph" begin
    g = fixture_cycle_reachable()
    costs = costs_new(ones(4), [1.0, 1.0, 0.0, 1.0], [2.0, 2.0, 2.0, 2.0])
    net = network_new(g, 1, 4, costs)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 4])
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 2.0)
    z_or = oracle_adv_cont_z(g, 1, 4, x, net.costs.c_hat, net.costs.d, NB_INCLUSION, 1, 2.0)
    sol = solve_adv(net, params, x, MIP_SOLVER)
    @test sol.status == ST_OK
    @test approxz(sol.z, z_or)
end

@testset "discrete ADV k = 0 dumps Delta on x" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    for delta in (0, 1, 2, 8)
        params = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 0; delta = delta)
        adv = solve_adv(net, params, x, COMB_SOLVER)
        @test adv.status == ST_OK
        @test adv.second.seq == x.seq
        @test adv.z == dump_delta_on_path(x, net.costs.c_hat, net.costs.d, delta)
        @test eval_worstcase(net, params, x, COMB_SOLVER) == adv.z
    end
end

@testset "discrete ADV Delta = 0 equals INC under c_hat" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    params = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 2; delta = 0)
    adv = solve_adv(net, params, x, COMB_SOLVER)
    inc = solve_inc(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 2), x, net.costs.c_hat, COMB_SOLVER)
    @test adv.status == ST_OK
    @test adv.z == inc.z
end

@testset "discrete ADV large Delta equals INC under c_hat+d" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    params = params_rrsp(U_DISC_BUDGET, NB_SYMDIFF, 2; delta = Int(net.graph.m))
    adv = solve_adv(net, params, x, COMB_SOLVER)
    inc = solve_inc(net, params_rrsp(U_INTERVAL, NB_SYMDIFF, 2), x, cmax_of(net), COMB_SOLVER)
    @test adv.status == ST_OK
    @test adv.z == inc.z
end

@testset "METHOD_COMB + continuous ADV with k > 0 is not implemented" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 2.0)
    sol = solve_adv(net, params, x, COMB_SOLVER)
    @test sol.status == ST_NOT_IMPL
end

@testset "mid-budget discrete ADV with k > 0 is not implemented" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    sol = solve_adv(net, params_rrsp(U_DISC_BUDGET, NB_EXCLUSION, 1; delta = 1), x, COMB_SOLVER)
    @test sol.status == ST_NOT_IMPL
end

@testset "ADV rejects a non s-t path" begin
    net = adv_diamond()
    bad = path_empty(Int(net.graph.m))
    @test_throws ArgumentError solve_adv(
        net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), bad, COMB_SOLVER
    )
end
