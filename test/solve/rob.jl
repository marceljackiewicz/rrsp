# Phase B: classical robust SP under interval uncertainty.
# First-stage cost is omitted from the robust objective.

@testset "ROB interval equals SP under c_hat + d" begin
    g = fixture_two_paths()
    costs = costs_new(
        [1000.0, 1000.0, 1000.0, 1000.0],
        [1.0, 1.0, 10.0, 10.0],
        [0.0, 0.0, 1.0, 1.0],
    )
    net = network_new(g, 1, 4, costs)
    params = params_rob(U_INTERVAL)
    wmax = net.costs.c_hat .+ net.costs.d
    sp = solve_sp(net, wmax, COMB_SOLVER)
    rob = solve_rob(net, params, COMB_SOLVER)
    @test rob.status == ST_OK
    @test rob.z == sp.z == 2.0
    @test rob.first.seq == sp.first.seq == Int32[1, 2]
    @test rob.z_second == rob.z
    @test rob.z_first == path_cost(rob.first, net.costs.C)
    @test rob.method_used == METHOD_COMB
    @test isempty(rob.second.seq)
end

@testset "ROB interval does not put C in the objective" begin
    g = fixture_unique_path(3)
    costs = costs_new([1.0e6, 1.0e6], [3.0, 4.0], [1.0, 1.0])
    net = network_new(g, 1, 3, costs)
    rob = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    test_ok_single_stage(rob; z = 9.0, seq = Int32[1, 2])
    @test rob.z_first == 2.0e6
    @test rob.z_second == 9.0
end

@testset "ROB interval with d = 0 equals SP under c_hat" begin
    g = fixture_dag_not_layered()
    net = fixture_network(g, 1, 3; C = 8.0, c_hat = 2.0, d = 0.0)
    rob = solve_rob(net, params_rob(U_INTERVAL), AUTO_SOLVER)
    sp = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    @test rob.z == sp.z
    @test rob.status == sp.status == ST_OK
end

@testset "ROB interval: t unreachable" begin
    net = fixture_network(fixture_disconnected_t(), 1, 3)
    sol = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    test_infeasible_solution(sol, 1)
    @test sol.method_used == METHOD_COMB
end

@testset "ROB interval on unique path" begin
    net, _ = load_instance("single_path")
    rob = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    # Four arcs, c_hat = 0, d = 2 → z = 8; unique path.
    test_ok_single_stage(rob; z = 8.0, seq = Int32[1, 2, 3, 4])
    @test eval_max(net, rob.first) == rob.z
    @test eval_nominal(net, rob.first) == 0.0
end

@testset "ROB interval on two_beads has a unique optimum" begin
    net, _ = load_instance("two_beads")
    rob = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    test_ok_single_stage(rob; z = 100.0, seq = Int32[1, 4])
    @test eval_max(net, rob.first) == 100.0
end

@testset "ROB interval on parallel arcs with tied max-costs" begin
    net, _ = load_instance("single_arc_paths")
    rob = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    # Every arc has c_hat + d = 6.
    test_ok_single_stage(rob; z = 6.0)
    @test length(rob.first.seq) == 1
    @test sum(Int, rob.first.chi) == 1
end

@testset "ROB nominal is SP under c_hat" begin
    g = fixture_two_beads()
    costs = costs_new(
        zeros(4),
        [50.0, 0.0, 20.0, 40.0],
        [10.0, 100.0, 30.0, 0.0],
    )
    net = network_new(g, 1, 3, costs)
    rob = solve_rob(net, params_rob(U_NOMINAL), COMB_SOLVER)
    sp = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    @test rob.status == ST_OK
    @test rob.z == sp.z == 20.0
    @test rob.first.seq == sp.first.seq
end

@testset "ROB METHOD_MIP with no optimizer is not implemented" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    sol = solve_rob(net, params_rob(U_INTERVAL), MIP_SOLVER_NO_OPT)
    @test sol.status == ST_NOT_IMPL
    @test sol.method_used == METHOD_MIP
end

@testset "ROB interval MIP matches COMB" begin
    net = fixture_network(fixture_two_paths(), 1, 4; c_hat = 1.0, d = 2.0)
    comb = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    mip = solve_rob(net, params_rob(U_INTERVAL), MIP_SOLVER)
    @test mip.status == ST_OK
    @test approxz(mip.z, comb.z)
    @test mip.method_used == METHOD_MIP
end

@testset "ROB rejects k != 0" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    p = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    @test_throws ArgumentError solve_rob(net, p, COMB_SOLVER)
end

@testset "ROB AUTO equals COMB on interval fixtures" begin
    g = fixture_layered(2, 3)
    net = fixture_network(g, 1, 3; c_hat = 1.0, d = 2.0)
    a = solve_rob(net, params_rob(U_INTERVAL), AUTO_SOLVER)
    c = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    @test a.status == c.status == ST_OK
    @test a.z == c.z
    @test a.method_used == c.method_used == METHOD_COMB
end

@testset "ROB continuous: Gamma = 0 equals SP under c_hat" begin
    g = fixture_two_beads()
    costs = costs_new(zeros(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    rob = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = 0.0), COMB_SOLVER)
    sp = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    @test rob.status == ST_OK
    @test rob.z == sp.z
    @test rob.method_used == METHOD_COMB
end

@testset "ROB continuous: large Gamma equals interval ROB" begin
    net, _ = load_instance("two_beads")
    gamma = sum(net.costs.d)
    cont = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = gamma), COMB_SOLVER)
    interval = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    @test cont.status == interval.status == ST_OK
    @test cont.z == interval.z
end

@testset "ROB continuous is nondecreasing in Gamma" begin
    net, _ = load_instance("two_beads")
    prev = -Inf
    for gamma in (0.0, 10.0, 50.0, 100.0, 1000.0)
        z = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = gamma), COMB_SOLVER).z
        @test z >= prev - 1e-12
        prev = z
    end
end

@testset "ROB continuous bottleneck vs spread" begin
    g = fixture_two_paths()
    # Path 1–2–4: cheap nominal, spread deviations. Path 1–3–4: expensive nominal, tiny d.
    costs = costs_new([0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 10.0, 10.0], [10.0, 10.0, 0.0, 1.0])
    net = network_new(g, 1, 4, costs)
    z0 = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = 0.0), COMB_SOLVER)
    z5 = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = 5.0), COMB_SOLVER)
    z50 = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = 50.0), COMB_SOLVER)
    @test z0.z == 2.0
    @test z0.first.seq == Int32[1, 2]
    @test z5.z == 7.0
    @test z5.first.seq == Int32[1, 2]
    @test z50.z == 21.0
    @test z50.first.seq == Int32[3, 4]
end

@testset "ROB continuous matches enumeration" begin
    g = fixture_two_beads()
    costs = costs_new(ones(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    for gamma in (0.0, 15.0, 80.0, 200.0)
        z_or = oracle_rob_cont_z(g, 1, 3, net.costs.c_hat, net.costs.d, gamma)
        sol = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = gamma), COMB_SOLVER)
        @test sol.status == ST_OK
        @test sol.z == z_or
        mip = solve_rob(net, params_rob(U_CONT_BUDGET; gamma = gamma), MIP_SOLVER)
        @test mip.status == ST_OK
        @test approxz(mip.z, z_or)
    end
end

@testset "ROB continuous rejects negative gamma" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    @test_throws ArgumentError solve_rob(net, params_rob(U_CONT_BUDGET; gamma = -1.0), COMB_SOLVER)
end

@testset "ROB interval matches SP under c_hat+d on random DAGs" begin
    rng = MersenneTwister(20260909)
    for trial in 1:25
        n = rand(rng, 2:7)
        m = rand(rng, 1:12)
        g = rand_digraph(rng, n, m; dag = true)
        C = [Float64(rand(rng, 0:5)) for _ in 1:m]
        c_hat = [Float64(rand(rng, 0:5)) for _ in 1:m]
        d = [Float64(rand(rng, 0:5)) for _ in 1:m]
        net = network_new(g, 1, n, costs_new(C, c_hat, d))
        rob = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
        sp = solve_sp(net, c_hat .+ d, COMB_SOLVER)
        @test rob.status == sp.status
        @test rob.z == sp.z
        if rob.status == ST_OK
            @test eval_max(net, rob.first) == rob.z
        end
    end
end

@testset "ROB discrete: Delta = 0 equals SP under c_hat" begin
    g = fixture_two_beads()
    costs = costs_new(zeros(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    rob = solve_rob(net, params_rob(U_DISC_BUDGET; delta = 0), COMB_SOLVER)
    sp = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    @test rob.status == ST_OK
    @test rob.z == sp.z
    @test rob.method_used == METHOD_COMB
end

@testset "ROB discrete: Delta >= m equals interval ROB" begin
    net, _ = load_instance("two_beads")
    disc = solve_rob(net, params_rob(U_DISC_BUDGET; delta = Int(net.graph.m)), COMB_SOLVER)
    interval = solve_rob(net, params_rob(U_INTERVAL), COMB_SOLVER)
    @test disc.status == interval.status == ST_OK
    @test disc.z == interval.z
end

@testset "ROB discrete is nondecreasing in Delta" begin
    net, _ = load_instance("two_beads")
    prev = -Inf
    for delta in 0:5
        z = solve_rob(net, params_rob(U_DISC_BUDGET; delta = delta), COMB_SOLVER).z
        @test z >= prev - 1e-12
        prev = z
    end
end

@testset "ROB discrete on parallel arcs hits the Delta largest d" begin
    net, _ = load_instance("single_arc_paths")
    # Each path is one arc; c_hat + d = 6 on every arc.
    z0 = solve_rob(net, params_rob(U_DISC_BUDGET; delta = 0), COMB_SOLVER)
    z1 = solve_rob(net, params_rob(U_DISC_BUDGET; delta = 1), COMB_SOLVER)
    @test z0.status == ST_OK
    @test z0.z == 1.0
    @test z1.status == ST_OK
    @test z1.z == 6.0
    @test length(z1.first.seq) == 1
end

@testset "ROB discrete matches enumeration" begin
    g = fixture_two_beads()
    costs = costs_new(ones(4), [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, costs)
    for delta in 0:3
        z_or = oracle_rob_disc_z(g, 1, 3, net.costs.c_hat, net.costs.d, delta)
        sol = solve_rob(net, params_rob(U_DISC_BUDGET; delta = delta), COMB_SOLVER)
        @test sol.status == ST_OK
        @test sol.z == z_or
        mip = solve_rob(net, params_rob(U_DISC_BUDGET; delta = delta), MIP_SOLVER)
        @test mip.status == ST_OK
        @test approxz(mip.z, z_or)
    end
end

@testset "ROB discrete MIP with no optimizer is not implemented" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    sol = solve_rob(net, params_rob(U_DISC_BUDGET; delta = 1), MIP_SOLVER_NO_OPT)
    @test sol.status == ST_NOT_IMPL
end

@testset "ROB discrete rejects negative delta" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    @test_throws ArgumentError solve_rob(net, params_rob(U_DISC_BUDGET; delta = -1), COMB_SOLVER)
end

@testset "ROB discrete matches enumeration on random DAGs" begin
    rng = MersenneTwister(20260910)
    for trial in 1:15
        n = rand(rng, 2:5)
        m = rand(rng, 1:7)
        g = rand_digraph(rng, n, m; dag = true)
        C = [Float64(rand(rng, 0:4)) for _ in 1:m]
        c_hat = [Float64(rand(rng, 0:4)) for _ in 1:m]
        d = [Float64(rand(rng, 0:4)) for _ in 1:m]
        net = network_new(g, 1, n, costs_new(C, c_hat, d))
        delta = rand(rng, 0:2)
        rob = solve_rob(net, params_rob(U_DISC_BUDGET; delta = delta), COMB_SOLVER)
        z_or = oracle_rob_disc_z(g, 1, n, c_hat, d, delta)
        @test rob.status == (z_or == Inf ? ST_INFEASIBLE : ST_OK)
        @test rob.z == z_or
    end
end
