# Recoverable shortest path under interval uncertainty (Phase C).

function rec_diamond()
    g = fixture_two_paths()
    # Cheap first-stage on the expensive second-stage path, and conversely.
    costs = costs_new(
        [1.0, 1.0, 20.0, 20.0],
        [10.0, 10.0, 1.0, 1.0],
        zeros(4),
    )
    return network_new(g, 1, 4, costs)
end

@testset "unique path: both stages that path" begin
    g = fixture_unique_path(4)
    net = fixture_network(g, 1, 4; C = 2.0, c_hat = 3.0, d = 1.0)
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    sol = solve_rec(net, params, COMB_SOLVER)
    @test sol.status == ST_OK
    @test sol.first.seq == sol.second.seq == Int32[1, 2, 3]
    @test sol.z == 2.0 * 3 + 4.0 * 3
    @test sol.z_first == 6.0
    @test sol.z_second == 12.0
end

@testset "k == 0 and C == 0 equals ROB interval" begin
    net = rec_diamond()
    costs0 = costs_new(zeros(4), net.costs.c_hat, net.costs.d)
    net0 = network_new(net.graph, net.s, net.t, costs0)
    rec = solve_rec(net0, params_rrsp(U_INTERVAL, NB_INCLUSION, 0), COMB_SOLVER)
    rob = solve_rob(net0, params_rob(U_INTERVAL), COMB_SOLVER)
    @test rec.status == rob.status == ST_OK
    @test rec.z == rob.z
    @test rec.first.seq == rec.second.seq
end

@testset "k == 0 and C != 0 uses the same path both stages" begin
    net = rec_diamond()
    rec = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 0), COMB_SOLVER)
    @test rec.status == ST_OK
    @test rec.first.seq == rec.second.seq
    # min of C(x)+(ĉ+d)(x): path 1-2-4 costs 2+20=22, path 1-3-4 costs 40+2=42.
    @test rec.z == 22.0
    @test rec.first.seq == Int32[1, 2]
end

@testset "k increasing is nonincreasing in z" begin
    net = rec_diamond()
    prev = Inf
    for k in 0:4
        z = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, k), COMB_SOLVER).z
        @test z <= prev
        prev = z
    end
    z0 = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 0), COMB_SOLVER).z
    z2 = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 2), COMB_SOLVER).z
    @test z2 == 4.0  # C on cheap-C path 2, second stage cheap-cmax path 2
    @test z2 < z0
end

@testset "large k: second stage is interval SP" begin
    net = rec_diamond()
    rec = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 4), COMB_SOLVER)
    sp2 = solve_sp(net, cmax_of(net), COMB_SOLVER)
    sp1 = solve_sp(net, net.costs.C, COMB_SOLVER)
    @test rec.status == ST_OK
    @test rec.z_second == sp2.z
    @test rec.z_first == sp1.z
    @test rec.z == sp1.z + sp2.z
end

@testset "layered neighborhood identity on Rec z" begin
    g = fixture_layered(2, 2)
    m = Int(g.m)
    C = ones(m)
    c_hat = [10.0, 1.0, 10.0, 1.0]
    net = network_new(g, 1, 3, costs_new(C, c_hat, zeros(m)))
    for k in 0:2
        zi = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, k), COMB_SOLVER).z
        ze = solve_rec(net, params_rrsp(U_INTERVAL, NB_EXCLUSION, k), COMB_SOLVER).z
        zs = solve_rec(net, params_rrsp(U_INTERVAL, NB_SYMDIFF, 2 * k), COMB_SOLVER).z
        @test zi == ze == zs
    end
end

@testset "infeasible network" begin
    net = fixture_network(fixture_disconnected_t(), 1, 3)
    sol = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), COMB_SOLVER)
    test_infeasible_solution(sol, 1)
end

@testset "METHOD_COMB on a cyclic graph is not implemented" begin
    g = fixture_cycle_reachable()
    net = fixture_network(g, 1, 4)
    sol = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), COMB_SOLVER)
    @test sol.status == ST_NOT_IMPL
end

@testset "METHOD_MIP with no optimizer is not implemented" begin
    net = rec_diamond()
    sol = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), MIP_SOLVER_NO_OPT)
    @test sol.status == ST_NOT_IMPL
end

@testset "continuous and discrete Rec are not implemented" begin
    net = rec_diamond()
    for U in (U_CONT_BUDGET, U_DISC_BUDGET)
        sol = solve_rec(net, params_rrsp(U, NB_INCLUSION, 1), COMB_SOLVER)
        @test sol.status == ST_NOT_IMPL
    end
end

@testset "DAG COMB vs MIP vs enumeration, all neighborhoods" begin
    net = rec_diamond()
    g = net.graph
    C = net.costs.C
    cm = cmax_of(net)
    @test is_dag(g)
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        for k in 0:3
            params = params_rrsp(U_INTERVAL, nb, k)
            z_or = oracle_rec_z(g, net.s, net.t, C, cm, nb, k)
            comb = solve_rec(net, params, COMB_SOLVER)
            mip = solve_rec(net, params, MIP_SOLVER)
            auto = solve_rec(net, params, TEST_SOLVER)
            @test comb.status == ST_OK
            @test mip.status == ST_OK
            @test auto.status == ST_OK
            @test comb.z == z_or
            @test approxz(mip.z, z_or)
            @test auto.z == z_or
            @test auto.method_used == METHOD_COMB
            @test mip.method_used == METHOD_MIP
        end
    end
end

@testset "cyclic Rec MIP vs enumeration" begin
    g = fixture_cycle_reachable()
    costs = costs_new([1.0, 1.0, 5.0, 1.0], [1.0, 1.0, 0.0, 1.0], zeros(4))
    net = network_new(g, 1, 4, costs)
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    z_or = oracle_rec_z(g, 1, 4, net.costs.C, cmax_of(net), NB_INCLUSION, 1)
    mip = solve_rec(net, params, MIP_SOLVER)
    auto = solve_rec(net, params, TEST_SOLVER)
    @test mip.status == ST_OK
    @test auto.status == ST_OK
    @test approxz(mip.z, z_or)
    @test approxz(auto.z, z_or)
    @test auto.method_used == METHOD_MIP
end

@testset "two_beads Rec inclusion vs oracle" begin
    net, _ = load_instance("two_beads")
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    z_or = oracle_rec_z(net.graph, net.s, net.t, net.costs.C, cmax_of(net), NB_INCLUSION, 1)
    comb = solve_rec(net, params, COMB_SOLVER)
    mip = solve_rec(net, params, MIP_SOLVER)
    @test comb.z == z_or
    @test approxz(mip.z, z_or)
end

@testset "enumeration on random DAGs" begin
    rng = MersenneTwister(20260909)
    for trial in 1:12
        n = rand(rng, 3:4)
        m = rand(rng, n:6)
        g = rand_digraph(rng, n, m; dag = true)
        C = [Float64(rand(rng, 0:4)) for _ in 1:Int(g.m)]
        c_hat = [Float64(rand(rng, 0:4)) for _ in 1:Int(g.m)]
        d = [Float64(rand(rng, 0:3)) for _ in 1:Int(g.m)]
        net = network_new(g, 1, n, costs_new(C, c_hat, d))
        isempty(enumerate_simple_st_paths(g, 1, n)) && continue
        nb = rand(rng, (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF))
        k = rand(rng, 0:2)
        params = params_rrsp(U_INTERVAL, nb, k)
        z_or = oracle_rec_z(g, 1, n, C, c_hat .+ d, nb, k)
        sol = solve_rec(net, params, COMB_SOLVER)
        if z_or == Inf
            @test sol.status == ST_INFEASIBLE
        else
            @test sol.status == ST_OK
            @test sol.z == z_or
        end
    end
end

@testset "ASP Rec matches DAG Rec and MIP on two_beads and series" begin
    net, _ = load_instance("two_beads")
    @test asp_decompose(net.graph, net.s, net.t) !== nothing
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        for k in 0:3
            params = params_rrsp(U_INTERVAL, nb, k)
            z_or = oracle_rec_z(net.graph, net.s, net.t, net.costs.C, cmax_of(net), nb, k)
            comb = solve_rec(net, params, COMB_SOLVER)
            mip = solve_rec(net, params, MIP_SOLVER)
            @test comb.status == ST_OK
            @test comb.z == z_or
            @test approxz(mip.z, z_or)
        end
    end
    g = fixture_unique_path(4)
    net2 = fixture_network(g, 1, 4; C = 2.0, c_hat = 3.0, d = 1.0)
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 2)
    sol = solve_rec(net2, params, COMB_SOLVER)
    @test sol.status == ST_OK
    @test sol.z == 2.0 * 3 + 4.0 * 3
    net3 = fixture_network(fixture_two_paths(), 1, 4; C = 1.0, c_hat = 2.0, d = 1.0)
    @test asp_decompose(net3.graph, net3.s, net3.t) !== nothing
    for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF)
        params3 = params_rrsp(U_INTERVAL, nb, 2)
        z_or = oracle_rec_z(net3.graph, 1, 4, net3.costs.C, cmax_of(net3), nb, 2)
        comb3 = solve_rec(net3, params3, COMB_SOLVER)
        @test comb3.status == ST_OK
        @test comb3.z == z_or
    end
end

@testset "wheatstone Rec uses DAG comb (not ASP)" begin
    g = fixture_wheatstone()
    @test asp_decompose(g, 1, 4) === nothing
    costs = costs_new(ones(5), [1.0, 2.0, 1.0, 2.0, 0.0], zeros(5))
    net = network_new(g, 1, 4, costs)
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    z_or = oracle_rec_z(g, 1, 4, net.costs.C, cmax_of(net), NB_INCLUSION, 1)
    sol = solve_rec(net, params, COMB_SOLVER)
    @test sol.status == ST_OK
    @test sol.z == z_or
end
