@testset "C1: exclusion ADV with k ≥ |x| is the max-min value, not the robust value" begin
    # Two parallel s–t arcs, ĉ = 0, d = 1, Γ = 1, x = arc 1, k = 1.
    # Every path is a neighbour of x, so the adversary raises one arc by 1/2
    # and the recovery still takes the other: value 1/2 (not 1).
    g = graph_new(2, Int32[1, 1], Int32[2, 2])
    net = network_new(g, 1, 2, costs_new([0.0, 0.0], [0.0, 0.0], [1.0, 1.0]))
    x = path_from_seq(g, Int32(1), Int32(2), Int32[1])
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_MIP)
    params = params_rrsp(U_CONT_BUDGET, NB_EXCLUSION, 1; gamma = 1.0)
    sol = solve_adv(net, params, x, solver)
    @test sol.status == ST_OK
    @test sol.z ≈ 0.5 atol = 1e-8
    for k in (1, 2, 5)
        sol = solve_adv(net, params_rrsp(U_CONT_BUDGET, NB_EXCLUSION, k; gamma = 1.0), x, solver)
        @test sol.z ≈ 0.5 atol = 1e-8
    end
    # Same instance through RRSP: the commitment is free, so z = 1/2.
    rr = solve_rrsp(net, params, solver)
    @test rr.status == ST_OK
    @test rr.z ≈ 0.5 atol = 1e-8
end

@testset "C2: asp_decompose terminates on cycles" begin
    # Detached directed cycle 3 ⇄ 4 next to the s–t arc.
    g = fixture_cycle_irrelevant()
    @test asp_decompose(g, 1, 2) === nothing
    # Isolated self-loop on vertex 3.
    g2 = graph_new(3, Int32[1, 3], Int32[2, 3])
    @test asp_decompose(g2, 1, 2) === nothing
    # Self-loop on an s–t path.
    g3 = fixture_self_loop()
    @test asp_decompose(g3, 1, 3) === nothing
    # The combinatorial solvers fall back instead of hanging.
    net = network_new(g2, 1, 2, uniform_costs(g2))
    x = path_from_seq(g2, Int32(1), Int32(2), Int32[1])
    sol = solve_inc(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), x, cmax_of(net), AUTO_SOLVER)
    @test sol.status == ST_OK
    # Rec on a cyclic digraph needs the MIP; without one it reports NOT_IMPL.
    p = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    @test solve_rec(net, p, AUTO_SOLVER).status == ST_NOT_IMPL
    @test solve_rec(net, p, TEST_SOLVER).status == ST_OK
end

@testset "C4: library MIPs work after a foreign HiGHS model" begin
    # A model with HiGHS's default thread setting runs first in this process.
    foreign = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(foreign)
    x = JuMP.@variable(foreign, 0 <= x <= 1)
    JuMP.@objective(foreign, Max, x)
    JuMP.optimize!(foreign)
    @test JuMP.termination_status(foreign) == JuMP.MOI.OPTIMAL

    g = fixture_two_paths()
    net = network_new(g, 1, 4, costs_new(ones(4), ones(4), ones(4)))
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 1.0)
    # Default (threads = 0): the library does not touch the thread count.
    sol = solve_rrsp(net, params, Solver(; optimizer = HiGHS.Optimizer, method = METHOD_MIP))
    @test sol.status == ST_OK
end
