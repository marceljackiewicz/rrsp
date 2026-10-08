# Cheap identities from the test plan (section 5.6), interval Rec/ROB/INC.

@testset "k = 0 and C = 0: Rec equals ROB" begin
    g = fixture_two_beads()
    costs = costs_new(zeros(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        rec = solve_rec(net, params_rrsp(U_INTERVAL, nb, 0), COMB_SOLVER)
        rob = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
        @test rec.z == rob.z
    end
end

@testset "Rec z is nonincreasing in k" begin
    net, _ = load_instance("two_beads")
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        prev = Inf
        for k in 0:4
            z = solve_rec(net, params_rrsp(U_INTERVAL, nb, k), COMB_SOLVER).z
            @test z <= prev + 1e-12
            prev = z
        end
    end
end

@testset "eval_max equals c_hat plus d inner products" begin
    g = fixture_unique_path(4)
    net = fixture_network(g, 1, 4; c_hat = 2.0, d = 3.0)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 3])
    @test eval_max(net, x) == path_cost(x, net.costs.c_hat) + path_cost(x, net.costs.d)
end

@testset "k = 0 and C = 0: RRSP continuous equals ROB continuous" begin
    g = fixture_two_beads()
    costs = costs_new(zeros(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    for gamma in (0.0, 40.0, 200.0)
        rrsp = solve_rrsp(net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 0; gamma = gamma), COMB_SOLVER)
        rob = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = gamma), COMB_SOLVER)
        @test rrsp.z == rob.z
    end
end

@testset "eval_worstcase equals solve_adv" begin
    g = fixture_two_paths()
    net = fixture_network(g, 1, 4; c_hat = 2.0, d = 3.0)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2])
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    @test eval_worstcase(net, params, x, COMB_SOLVER) == solve_adv(net, params, x, COMB_SOLVER).z
end

@testset "k = 0 and C = 0: RRSP discrete equals ROB discrete" begin
    g = fixture_two_beads()
    costs = costs_new(zeros(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    for delta in (0, 1, 4)
        rrsp = solve_rrsp(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 0; delta = delta), COMB_SOLVER)
        rob = solve_rob(net, params_rob(U_DISC_BUDGET; delta = delta), COMB_SOLVER)
        @test rrsp.z == rob.z
    end
end
