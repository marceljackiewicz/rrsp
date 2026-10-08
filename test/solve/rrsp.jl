function rrsp_beads()
    g = fixture_two_beads()
    costs = costs_new(
        [1.0, 2.0, 1.0, 2.0],
        [50.0, 0.0, 20.0, 40.0],
        [10.0, 100.0, 30.0, 0.0],
    )
    return network_new(g, 1, 3, costs)
end

@testset "RRSP interval dispatches to Rec" begin
    net = rrsp_beads()
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    a = solve_rrsp(net, params, COMB_SOLVER)
    b = solve_rec(net, params, COMB_SOLVER)
    @test a.status == b.status == ST_OK
    @test a.z == b.z
end

@testset "RRSP k = 0 and C = 0 equals ROB continuous" begin
    net = rrsp_beads()
    costs0 = costs_new(zeros(4), net.costs.c_hat, net.costs.d)
    net0 = network_new(net.graph, net.s, net.t, costs0)
    for gamma in (0.0, 25.0, 200.0)
        rrsp = solve_rrsp(net0, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 0; gamma = gamma), COMB_SOLVER)
        rob = solve_rob(net0, params_rob(U_CONT_BUDGET; gamma = gamma), COMB_SOLVER)
        @test rrsp.status == rob.status == ST_OK
        @test rrsp.z == rob.z
        @test rrsp.first.seq == rrsp.second.seq
    end
end

@testset "RRSP Gamma = 0 equals Rec under c_hat" begin
    net = rrsp_beads()
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 0.0)
    a = solve_rrsp(net, params, COMB_SOLVER)
    b = solve_rec(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), COMB_SOLVER)
    @test a.status == b.status == ST_OK
    @test a.z == b.z
end

@testset "RRSP large Gamma equals Rec under c_hat+d" begin
    net = rrsp_beads()
    gamma = sum(net.costs.d)
    params = params_rrsp(U_CONT_BUDGET, NB_EXCLUSION, 1; gamma = gamma)
    a = solve_rrsp(net, params, COMB_SOLVER)
    b = solve_rec(net, params_rrsp(U_INTERVAL, NB_EXCLUSION, 1), COMB_SOLVER)
    @test a.status == b.status == ST_OK
    @test a.z == b.z
end

@testset "RRSP unique path both stages" begin
    g = fixture_unique_path(4)
    net = fixture_network(g, 1, 4; C = 2.0, c_hat = 3.0, d = 1.0)
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 2.0)
    sol = solve_rrsp(net, params, COMB_SOLVER)
    @test sol.status == ST_OK
    @test sol.first.seq == sol.second.seq == Int32[1, 2, 3]
    @test sol.z == 2.0 * 3 + dump_gamma_on_path(sol.first, net.costs.c_hat, net.costs.d, 2.0)
end

@testset "RRSP k increasing is nonincreasing in z" begin
    net = rrsp_beads()
    prev = Inf
    for k in 0:3
        z = solve_rrsp(
            net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, k; gamma = 20.0), TEST_SOLVER
        ).z
        @test z <= prev + 1e-8
        prev = z
    end
end

@testset "RRSP continuous MIP vs enumeration on two_beads" begin
    net = rrsp_beads()
    @test asp_decompose(net.graph, net.s, net.t) !== nothing
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        for k in 1:2
            params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = 25.0)
            z_or = oracle_rrsp_cont_z(
                net.graph, 1, 3, net.costs.C, net.costs.c_hat, net.costs.d, nb, k, 25.0
            )
            sol = solve_rrsp(net, params, TEST_SOLVER)
            @test sol.status == ST_OK
            @test approxz(sol.z, z_or)
        end
    end
end

@testset "RRSP discrete k = 0 and C = 0 equals ROB discrete" begin
    net = rrsp_beads()
    costs0 = costs_new(zeros(4), net.costs.c_hat, net.costs.d)
    net0 = network_new(net.graph, net.s, net.t, costs0)
    for delta in (0, 1, 4)
        rrsp = solve_rrsp(net0, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 0; delta = delta), COMB_SOLVER)
        rob = solve_rob(net0, params_rob(U_DISC_BUDGET; delta = delta), COMB_SOLVER)
        @test rrsp.status == rob.status == ST_OK
        @test rrsp.z == rob.z
        @test rrsp.first.seq == rrsp.second.seq
    end
end

@testset "RRSP discrete Delta = 0 equals Rec under c_hat" begin
    net = rrsp_beads()
    params = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 1; delta = 0)
    a = solve_rrsp(net, params, COMB_SOLVER)
    b = solve_rec(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), COMB_SOLVER)
    @test a.status == b.status == ST_OK
    @test a.z == b.z
end

@testset "RRSP discrete large Delta equals Rec under c_hat+d" begin
    net = rrsp_beads()
    params = params_rrsp(U_DISC_BUDGET, NB_EXCLUSION, 1; delta = Int(net.graph.m))
    a = solve_rrsp(net, params, COMB_SOLVER)
    b = solve_rec(net, params_rrsp(U_INTERVAL, NB_EXCLUSION, 1), COMB_SOLVER)
    @test a.status == b.status == ST_OK
    @test a.z == b.z
end

@testset "RRSP discrete MIP vs enumeration on two_beads k = 0" begin
    net = rrsp_beads()
    for delta in 0:2
        z_or = oracle_rrsp_disc_z(
            net.graph, 1, 3, net.costs.C, net.costs.c_hat, net.costs.d, NB_INCLUSION, 0, delta
        )
        sol = solve_rrsp(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 0; delta = delta), TEST_SOLVER)
        @test sol.status == ST_OK
        @test approxz(sol.z, z_or)
        mip = solve_rrsp(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 0; delta = delta), MIP_SOLVER)
        @test mip.status == ST_OK
        @test approxz(mip.z, z_or)
    end
end

@testset "RRSP continuous COMB on a mid-budget instance is not implemented" begin
    net = rrsp_beads()
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 25.0)
    sol = solve_rrsp(net, params, COMB_SOLVER)
    @test sol.status == ST_NOT_IMPL
end

@testset "RRSP discrete mid-budget with k > 0 is not implemented" begin
    net = rrsp_beads()
    sol = solve_rrsp(net, params_rrsp(U_DISC_BUDGET, NB_EXCLUSION, 1; delta = 1), MIP_SOLVER)
    @test sol.status == ST_NOT_IMPL
end

@testset "RRSP MIP with no optimizer is not implemented" begin
    net = rrsp_beads()
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 25.0)
    sol = solve_rrsp(net, params, MIP_SOLVER_NO_OPT)
    @test sol.status == ST_NOT_IMPL
end

@testset "RRSP tiny time limit does not crash" begin
    net = rrsp_beads()
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 25.0)
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_MIP, silent = true, time_limit = 1e-6)
    sol = solve_rrsp(net, params, solver)
    @test sol.status == ST_OK || sol.status == ST_TIME_LIMIT
end

@testset "layered RRSP continuous neighborhood identity" begin
    g = fixture_layered(2, 2)
    m = Int(g.m)
    C = [Float64(a) for a in 1:m]
    c_hat = [Float64(m - a + 1) for a in 1:m]
    d = fill(3.0, m)
    net = network_new(g, 1, 3, costs_new(C, c_hat, d))
    k = 1
    gamma = 4.0
    zi = solve_rrsp(net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, k; gamma = gamma), TEST_SOLVER)
    ze = solve_rrsp(net, params_rrsp(U_CONT_BUDGET, NB_EXCLUSION, k; gamma = gamma), TEST_SOLVER)
    zs = solve_rrsp(net, params_rrsp(U_CONT_BUDGET, NB_SYMDIFF, 2 * k; gamma = gamma), TEST_SOLVER)
    @test zi.status == ze.status == zs.status == ST_OK
    @test approxz(zi.z, ze.z)
    @test approxz(ze.z, zs.z)
end
