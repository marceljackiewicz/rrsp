# Incremental shortest path (Phase C).

function diamond_net()
    g = fixture_two_paths()
    # Arcs: 1: 1→2, 2: 2→4, 3: 1→3, 4: 3→4.
    costs = costs_new(
        [1.0, 1.0, 1.0, 1.0],
        [10.0, 10.0, 1.0, 1.0],
        zeros(4),
    )
    return network_new(g, 1, 4, costs)
end

function long_short_net()
    g = fixture_dag_not_layered()
    # Arcs: 1: 1→3, 2: 1→2, 3: 2→3.
    costs = costs_new([1.0, 1.0, 1.0], [1.0, 10.0, 10.0], zeros(3))
    return network_new(g, 1, 3, costs)
end

@testset "k == 0 recovers x" begin
    net = diamond_net()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    w = net.costs.c_hat
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        params = params_rrsp(U_NOMINAL, nb, 0)
        sol = solve_inc(net, params, x, w, COMB_SOLVER)
        @test sol.status == ST_OK
        @test sol.first.seq == x.seq
        @test sol.second.seq == x.seq
        @test sol.z == path_cost(x, w) == 20.0
        @test sol.method_used == METHOD_COMB
    end
end

@testset "large k equals unconstrained SP" begin
    net = diamond_net()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    w = net.costs.c_hat
    params = params_rrsp(U_NOMINAL, NB_INCLUSION, 4)
    sol = solve_inc(net, params, x, w, COMB_SOLVER)
    sp = solve_sp(net, w, COMB_SOLVER)
    @test sol.status == ST_OK
    @test sol.z == sp.z == 2.0
    @test sol.second.seq == Int32[3, 4]
end

@testset "layered: inclusion k = exclusion k = symdiff 2k" begin
    g = fixture_layered(3, 2)
    m = Int(g.m)
    w = [Float64(a) for a in 1:m]
    net = network_new(g, 1, 4, costs_new(ones(m), w, zeros(m)))
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 3, 5])
    for k in (0, 1, 2, 3)
        zi = solve_inc(net, params_rrsp(U_NOMINAL, NB_INCLUSION, k), x, w, COMB_SOLVER).z
        ze = solve_inc(net, params_rrsp(U_NOMINAL, NB_EXCLUSION, k), x, w, COMB_SOLVER).z
        zs = solve_inc(net, params_rrsp(U_NOMINAL, NB_SYMDIFF, 2 * k), x, w, COMB_SOLVER).z
        @test zi == ze == zs
        @test zi == oracle_inc_z(g, 1, 4, x, w, NB_INCLUSION, k)
    end
end

@testset "non-layered: inclusion and exclusion differ" begin
    net = long_short_net()
    g = net.graph
    x = path_from_seq(g, net.s, net.t, Int32[2, 3])
    w = net.costs.c_hat
    zi = solve_inc(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), x, w, COMB_SOLVER)
    ze = solve_inc(net, params_rrsp(U_NOMINAL, NB_EXCLUSION, 1), x, w, COMB_SOLVER)
    @test zi.status == ST_OK
    @test ze.status == ST_OK
    @test zi.z == 1.0
    @test ze.z == 20.0
    @test zi.z != ze.z
end

@testset "x that is not an s–t path throws" begin
    net = diamond_net()
    bad = path_empty(Int(net.graph.m))
    @test_throws ArgumentError solve_inc(
        net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), bad, net.costs.c_hat, COMB_SOLVER
    )
    wrong = path_empty(3)
    @test_throws ArgumentError solve_inc(
        net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), wrong, net.costs.c_hat, COMB_SOLVER
    )
end

@testset "cycle + exclusion k >= 1: recovery is simple" begin
    g = fixture_cycle_reachable()
    net = fixture_network(g, 1, 4; c_hat = 1.0)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 4])
    w = [1.0, 1.0, 0.0, 1.0]
    sol = solve_inc(net, params_rrsp(U_NOMINAL, NB_EXCLUSION, 1), x, w, MIP_SOLVER)
    @test sol.status == ST_OK
    @test sol.second.seq == path_from_chi(g, net.s, net.t, sol.second.chi).seq
    @test length(sol.second.seq) == length(unique(sol.second.seq))
end

@testset "METHOD_COMB + exclusion on a cyclic graph is not implemented" begin
    g = fixture_cycle_reachable()
    net = fixture_network(g, 1, 4)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 4])
    sol = solve_inc(net, params_rrsp(U_NOMINAL, NB_EXCLUSION, 1), x, fill(1.0, 4), COMB_SOLVER)
    @test sol.status == ST_NOT_IMPL
    sol2 = solve_inc(net, params_rrsp(U_NOMINAL, NB_SYMDIFF, 2), x, fill(1.0, 4), COMB_SOLVER)
    @test sol2.status == ST_NOT_IMPL
end

@testset "METHOD_MIP with no optimizer is not implemented" begin
    net = diamond_net()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    sol = solve_inc(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), x, net.costs.c_hat, MIP_SOLVER_NO_OPT)
    @test sol.status == ST_NOT_IMPL
end

@testset "inclusion COMB vs MIP vs enumeration" begin
    net = diamond_net()
    g = net.graph
    x = path_from_seq(g, net.s, net.t, Int32[1, 2])
    w = net.costs.c_hat
    for k in 0:3
        params = params_rrsp(U_NOMINAL, NB_INCLUSION, k)
        z_or = oracle_inc_z(g, net.s, net.t, x, w, NB_INCLUSION, k)
        comb = solve_inc(net, params, x, w, COMB_SOLVER)
        mip = solve_inc(net, params, x, w, MIP_SOLVER)
        auto = solve_inc(net, params, x, w, TEST_SOLVER)
        @test comb.status == ST_OK
        @test mip.status == ST_OK
        @test auto.status == ST_OK
        @test comb.z == z_or
        @test approxz(mip.z, z_or)
        @test auto.z == z_or
        @test auto.method_used == METHOD_COMB
        @test mip.method_used == METHOD_MIP
        @test eval_recovered(net, params, x, w, COMB_SOLVER) == comb.z
    end
end

@testset "DAG exclusion and symdiff COMB vs MIP vs enumeration" begin
    net = diamond_net()
    g = net.graph
    @test is_dag(g)
    x = path_from_seq(g, net.s, net.t, Int32[1, 2])
    w = net.costs.c_hat
    for (nb, ks) in ((NB_EXCLUSION, 0:2), (NB_SYMDIFF, 0:4))
        for k in ks
            params = params_rrsp(U_NOMINAL, nb, k)
            z_or = oracle_inc_z(g, net.s, net.t, x, w, nb, k)
            comb = solve_inc(net, params, x, w, COMB_SOLVER)
            mip = solve_inc(net, params, x, w, MIP_SOLVER)
            @test comb.status == ST_OK
            @test mip.status == ST_OK
            @test comb.z == z_or
            @test approxz(mip.z, z_or)
            @test comb.method_used == METHOD_COMB
        end
    end
end

@testset "enumeration oracle on random DAGs" begin
    rng = MersenneTwister(20260909)
    for trial in 1:20
        n = rand(rng, 3:5)
        m = rand(rng, n:8)
        g = rand_digraph(rng, n, m; dag = true)
        w = [Float64(rand(rng, 0:5)) for _ in 1:Int(g.m)]
        net = fixture_network(g, 1, n)
        paths = enumerate_simple_st_paths(g, 1, n)
        isempty(paths) && continue
        x = path_from_seq(g, Int32(1), Int32(n), paths[1])
        k = rand(rng, 0:3)
        nb = rand(rng, (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF))
        params = params_rrsp(U_NOMINAL, nb, k)
        z_or = oracle_inc_z(g, 1, n, x, w, nb, k)
        sol = solve_inc(net, params, x, w, COMB_SOLVER)
        if z_or == Inf
            @test sol.status == ST_INFEASIBLE
        else
            @test sol.status == ST_OK
            @test sol.z == z_or
        end
    end
end

@testset "weight length mismatch throws" begin
    net = diamond_net()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    @test_throws ArgumentError solve_inc(
        net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), x, [1.0], COMB_SOLVER
    )
end

@testset "ASP Inc matches oracle and MIP on two_beads and diamond" begin
    net, _ = load_instance("two_beads")
    @test asp_decompose(net.graph, net.s, net.t) !== nothing
    w = net.costs.c_hat
    xseq = solve_sp(net, w, COMB_SOLVER).first.seq
    x = path_from_seq(net.graph, net.s, net.t, xseq)
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        for k in 0:3
            params = params_rrsp(U_NOMINAL, nb, k)
            z_or = oracle_inc_z(net.graph, net.s, net.t, x, w, nb, k)
            comb = solve_inc(net, params, x, w, COMB_SOLVER)
            mip = solve_inc(net, params, x, w, MIP_SOLVER)
            @test comb.status == ST_OK
            @test comb.z == z_or
            @test approxz(mip.z, z_or)
        end
    end
    net2 = diamond_net()
    @test asp_decompose(net2.graph, net2.s, net2.t) !== nothing
    x2 = path_from_seq(net2.graph, net2.s, net2.t, Int32[1, 2])
    w2 = net2.costs.c_hat
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        params = params_rrsp(U_NOMINAL, nb, 2)
        z_or = oracle_inc_z(net2.graph, net2.s, net2.t, x2, w2, nb, 2)
        comb = solve_inc(net2, params, x2, w2, COMB_SOLVER)
        @test comb.status == ST_OK
        @test comb.z == z_or
    end
end

@testset "wheatstone Inc uses DAG comb (not ASP)" begin
    g = fixture_wheatstone()
    @test asp_decompose(g, 1, 4) === nothing
    costs = costs_new(ones(5), [1.0, 2.0, 1.0, 2.0, 0.0], zeros(5))
    net = network_new(g, 1, 4, costs)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 3])
    w = net.costs.c_hat
    params = params_rrsp(U_NOMINAL, NB_SYMDIFF, 2)
    z_or = oracle_inc_z(g, 1, 4, x, w, NB_SYMDIFF, 2)
    sol = solve_inc(net, params, x, w, COMB_SOLVER)
    @test sol.status == ST_OK
    @test sol.z == z_or
end
