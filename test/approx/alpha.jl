function approx_beads()
    g = fixture_two_beads()
    costs = costs_new(
        [1.0, 2.0, 1.0, 2.0],
        [50.0, 10.0, 20.0, 40.0],
        [10.0, 100.0, 30.0, 0.0],
    )
    return network_new(g, 1, 3, costs)
end

function approx_zero_nominal()
    g = fixture_two_paths()
    costs = costs_new(
        [1.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 0.0],
    )
    return network_new(g, 1, 4, costs)
end

@testset "cost_structure_alpha on single_arc_paths is 1/6" begin
    net, _ = load_instance("single_arc_paths")
    a = cost_structure_alpha(net)
    @test a !== nothing
    @test a == 1.0 / 6.0
end

@testset "cost_structure_alpha is nothing when a zero-nominal arc has positive d" begin
    net = approx_zero_nominal()
    @test cost_structure_alpha(net) === nothing
end

@testset "cost_structure_alpha is 1 when d = 0" begin
    g = fixture_unique_path(3)
    net = fixture_network(g, 1, 3; c_hat = 4.0, d = 0.0)
    @test cost_structure_alpha(net) == 1.0
end

@testset "1/alpha heuristic on single_arc_paths has realized ratio 8/7" begin
    net, _ = load_instance("single_arc_paths")
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 0; gamma = 100.0)
    heur = approx_rrsp(net, params, COMB_SOLVER)
    opt = solve_rrsp(net, params, COMB_SOLVER)
    @test heur.status == ST_OK
    @test opt.status == ST_OK
    @test heur.z == 8.0
    @test opt.z == 7.0
    a = cost_structure_alpha(net)
    @test heur.z / opt.z <= 1.0 / a + 1e-12
    @test heur.z / opt.z == 8.0 / 7.0
    fac = approx_bound_factors(net, params)
    @test fac.alpha == 1.0 / 6.0
    @test fac.bound_alpha == 6.0
    @test heur.z / opt.z <= fac.bound_alpha + 1e-12
end

@testset "1/alpha heuristic vs MIP on alpha_ok two_beads, continuous" begin
    net = approx_beads()
    a = cost_structure_alpha(net)
    @test a !== nothing
    @test a == 10.0 / 110.0
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 25.0)
    heur = approx_rrsp(net, params, TEST_SOLVER)
    opt = solve_rrsp(net, params, TEST_SOLVER)
    @test heur.status == ST_OK
    @test opt.status == ST_OK
    @test heur.z + 1e-8 >= opt.z
    @test heur.z / opt.z <= 1.0 / a + 1e-8
    m = Float64(net.graph.m)
    @test heur.z / opt.z <= m + 1e-8
    Dcap = 10.0 + 25.0 + 25.0
    kappa = 25.0 / Dcap
    @test heur.z / opt.z <= 1.0 / kappa + 1e-8
    fac = approx_bound_factors(net, params, TEST_SOLVER)
    @test isapprox(fac.kappa, kappa; atol = 1e-12)
    @test heur.z / opt.z <= fac.bound + 1e-8
    @test isfinite(fac.nu) || isnan(fac.nu)
end

@testset "1/alpha heuristic vs enumeration on alpha_ok two_beads, discrete" begin
    net = approx_beads()
    a = cost_structure_alpha(net)
    params = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 1; delta = 1)
    heur = approx_rrsp(net, params, COMB_SOLVER)
    z_or = oracle_rrsp_disc_z(
        net.graph, 1, 3, net.costs.C, net.costs.c_hat, net.costs.d, NB_INCLUSION, 1, 1
    )
    @test heur.status == ST_OK
    @test isfinite(heur.z)
    @test heur.z + 1e-8 >= z_or
    @test heur.z / z_or <= 1.0 / a + 1e-8
end

@testset "zero-nominal positive deviation: alpha missing, heuristic still finite" begin
    net = approx_zero_nominal()
    @test cost_structure_alpha(net) === nothing
    params = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 1; delta = 1)
    heur = approx_rrsp(net, params, COMB_SOLVER)
    @test heur.status == ST_OK
    @test isfinite(heur.z)
    @test heur.z > 0.0
end

@testset "approx_rrsp interval is exact Rec" begin
    net = approx_beads()
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    a = approx_rrsp(net, params, COMB_SOLVER)
    b = solve_rec(net, params, COMB_SOLVER)
    @test a.status == b.status == ST_OK
    @test a.z == b.z
end

@testset "approx_rrsp COMB on a cyclic graph with k > 0 is not implemented" begin
    g = fixture_cycle_reachable()
    costs = costs_new(ones(4), ones(4), ones(4))
    net = network_new(g, 1, 4, costs)
    sol = approx_rrsp(net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 1.0), COMB_SOLVER)
    @test sol.status == ST_NOT_IMPL
end

@testset "approx_rrsp rejects negative gamma" begin
    net = approx_beads()
    @test_throws ArgumentError approx_rrsp(
        net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = -1.0), COMB_SOLVER
    )
end
