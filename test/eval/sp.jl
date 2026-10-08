# Phase B–C: nonnegative shortest path. METHOD_MIP with no optimizer must
# return ST_NOT_IMPL, not throw. With an optimizer, MIP matches COMB.

@testset "METHOD_MIP with optimizer === nothing is not implemented" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    sol = solve_sp(net, [1.0, 1.0], MIP_SOLVER_NO_OPT)
    @test sol.status == ST_NOT_IMPL
    @test sol.method_used == METHOD_MIP
    @test sol.z == Inf
    @test isempty(sol.first.seq)
end

@testset "METHOD_COMB is always available for SP" begin
    for (name, g) in all_fixture_graphs()
        @testset "$name" begin
            s = Int32(1)
            t = g.n == 1 ? Int32(1) : g.n
            # Prefer documented terminals when the fixture has them in mind.
            if name == "disconnected_t"
                t = Int32(3)
            elseif name == "isolated_s"
                t = Int32(3)
            elseif name == "empty_st"
                t = Int32(2)
            end
            m = Int(g.m)
            net = fixture_network(g, s, t)
            w = fill(1.0, m)
            sol = solve_sp(net, w, COMB_SOLVER)
            @test sol.status != ST_NOT_IMPL
            @test sol.status == ST_OK || sol.status == ST_INFEASIBLE
            @test sol.method_used == METHOD_COMB
        end
    end
end

@testset "unique path" begin
    g = fixture_unique_path(5)
    net = fixture_network(g, 1, 5; C = 100.0, c_hat = 7.0)
    w = [10.0, 20.0, 30.0, 40.0]
    for slv in (COMB_SOLVER, AUTO_SOLVER)
        sol = solve_sp(net, w, slv)
        test_ok_single_stage(sol; z = 100.0, seq = Int32[1, 2, 3, 4])
        @test sol.first.chi == UInt8[0x01, 0x01, 0x01, 0x01]
        @test sol.z_first == 400.0
        @test sol.z_second == 0.0
        @test path_cost(sol.first, w) == sol.z
    end
end

@testset "parallel arcs: unique minimum" begin
    g = fixture_parallel_arcs()
    net = fixture_network(g, 1, 2)
    w = [5.0, 4.0, 3.0, 2.0, 1.0]
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 1.0, seq = Int32[5])
    @test sum(Int, sol.first.chi) == 1
    @test sol.first.chi[5] == 0x01
end

@testset "parallel arcs: tied weights, unique z, single arc" begin
    g = fixture_parallel_arcs()
    net = fixture_network(g, 1, 2)
    w = fill(3.0, 5)
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 3.0)
    @test length(sol.first.seq) == 1
    @test sum(Int, sol.first.chi) == 1
end

@testset "two paths, different costs" begin
    g = fixture_two_paths()
    net = fixture_network(g, 1, 4)
    # Arcs: 1: 1→2, 2: 2→4, 3: 1→3, 4: 3→4.
    w = [100.0, 100.0, 1.0, 1.0]
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 2.0, seq = Int32[3, 4])
    @test sol.first.chi == UInt8[0x00, 0x00, 0x01, 0x01]
end

@testset "two paths, equal costs: do not assert seq" begin
    g = fixture_two_paths()
    net = fixture_network(g, 1, 4)
    w = fill(1.0, 4)
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 2.0)
    @test length(sol.first.seq) == 2
    @test sum(Int, sol.first.chi) == 2
end

@testset "t unreachable" begin
    g = fixture_disconnected_t()
    net = fixture_network(g, 1, 3)
    sol = solve_sp(net, [1.0], COMB_SOLVER)
    test_infeasible_solution(sol, 1)
    @test sol.method_used == METHOD_COMB
end

@testset "m = 0, s distinct from t" begin
    g = fixture_empty_st()
    net = fixture_network(g, 1, 2)
    sol = solve_sp(net, Float64[], COMB_SOLVER)
    test_infeasible_solution(sol, 0)
end

@testset "s equals t is the trivial path of cost 0" begin
    g = fixture_unique_path(4)
    net = fixture_network(g, 2, 2)
    w = [1.0, 1.0, 1.0]
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 0.0, seq = Int32[])
    @test all(==(0x00), sol.first.chi)
    @test sol.z_first == 0.0
end

@testset "s equals t on an empty one-vertex graph" begin
    g = fixture_graph_n1_empty()
    net = fixture_network(g, 1, 1)
    sol = solve_sp(net, Float64[], COMB_SOLVER)
    test_ok_single_stage(sol; z = 0.0, seq = Int32[])
end

@testset "zero weights: some path of cost 0" begin
    g = fixture_two_paths()
    net = fixture_network(g, 1, 4)
    w = zeros(4)
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 0.0)
    @test !isempty(sol.first.seq)
    @test path_cost(sol.first, w) == 0.0
end

@testset "zero weights on a unique path still reports that path" begin
    g = fixture_unique_path(4)
    net = fixture_network(g, 1, 4)
    w = zeros(3)
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 0.0, seq = Int32[1, 2, 3])
end

@testset "DAG scan and Dijkstra agree on a unique DAG optimum" begin
    g = fixture_unique_path(6)
    @test is_dag(g)
    w = [2.0, 5.0, 1.0, 8.0, 3.0]
    st_d, p_d, z_d = Rrsp._shortest_path_dijkstra(g, Int32(1), Int32(6), w)
    st_g, p_g, z_g = Rrsp._shortest_path_dag(g, Int32(1), Int32(6), w)
    @test st_d == ST_OK
    @test st_g == ST_OK
    @test z_d == z_g == 19.0
    @test p_d.seq == p_g.seq == Int32[1, 2, 3, 4, 5]
end

@testset "DAG vs cyclic graph, same unique s–t optimum" begin
    g_dag = graph_new(2, Int32[1], Int32[2])
    g_cyc = fixture_cycle_irrelevant()
    @test is_dag(g_dag)
    @test !is_dag(g_cyc)
    net_d = fixture_network(g_dag, 1, 2)
    net_c = fixture_network(g_cyc, 1, 2)
    sol_d = solve_sp(net_d, [4.0], COMB_SOLVER)
    sol_c = solve_sp(net_c, [4.0, 99.0, 99.0], COMB_SOLVER)
    test_ok_single_stage(sol_d; z = 4.0, seq = Int32[1])
    test_ok_single_stage(sol_c; z = 4.0, seq = Int32[1])
end

@testset "DAG shortest path rejects a cyclic graph" begin
    g = fixture_cycle_reachable()
    @test !is_dag(g)
    @test_throws ArgumentError Rrsp._shortest_path_dag(g, Int32(1), Int32(4), [1.0, 1.0, 1.0, 1.0])
end

@testset "self-loop is not used on the unique simple path" begin
    g = fixture_self_loop()
    net = fixture_network(g, 1, 3)
    w = [1.0, 0.0, 1.0]
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 2.0, seq = Int32[1, 3])
end

@testset "cycle on the route: unique simple path" begin
    g = fixture_cycle_reachable()
    net = fixture_network(g, 1, 4)
    w = [1.0, 1.0, 0.0, 1.0]
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 3.0, seq = Int32[1, 2, 4])
end

@testset "non-layered DAG prefers the short arc" begin
    g = fixture_dag_not_layered()
    net = fixture_network(g, 1, 3)
    w = fill(1.0, 3)
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 1.0, seq = Int32[1])
end

@testset "diamond plus cheap chord" begin
    g = fixture_diamond_plus_chord()
    net = fixture_network(g, 1, 4)
    w = [10.0, 10.0, 10.0, 10.0, 1.0]
    sol = solve_sp(net, w, COMB_SOLVER)
    test_ok_single_stage(sol; z = 1.0, seq = Int32[5])
end

@testset "isolated origin cannot reach t" begin
    g = fixture_isolated_s()
    net = fixture_network(g, 1, 3)
    sol = solve_sp(net, [1.0], COMB_SOLVER)
    test_infeasible_solution(sol, 1)
end

@testset "isolated pair that is connected is feasible" begin
    g = fixture_isolated_s()
    net = fixture_network(g, 2, 3)
    sol = solve_sp(net, [9.0], COMB_SOLVER)
    test_ok_single_stage(sol; z = 9.0, seq = Int32[1])
end

@testset "weight length mismatch throws" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    @test_throws ArgumentError solve_sp(net, [1.0], COMB_SOLVER)
    @test_throws ArgumentError solve_sp(net, [1.0, 1.0, 1.0], COMB_SOLVER)
end

@testset "negative weight throws" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    @test_throws ArgumentError solve_sp(net, [1.0, -0.1], COMB_SOLVER)
end

@testset "solve_sp does not mutate the weight vector" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    w = [3.0, 4.0]
    w0 = copy(w)
    solve_sp(net, w, COMB_SOLVER)
    @test w == w0
end

@testset "AUTO and COMB agree on fixtures against the enumeration oracle" begin
    rng = MersenneTwister(20260909)
    named = [
        ("unique_path", fixture_unique_path(), Int32(1), Int32(5)),
        ("two_paths", fixture_two_paths(), Int32(1), Int32(4)),
        ("two_beads", fixture_two_beads(), Int32(1), Int32(3)),
        ("dag_not_layered", fixture_dag_not_layered(), Int32(1), Int32(3)),
        ("layered", fixture_layered(3, 2), Int32(1), Int32(4)),
        ("parallel", fixture_parallel_arcs(), Int32(1), Int32(2)),
        ("cycle_reachable", fixture_cycle_reachable(), Int32(1), Int32(4)),
        ("cycle_irrelevant", fixture_cycle_irrelevant(), Int32(1), Int32(2)),
        ("self_loop", fixture_self_loop(), Int32(1), Int32(3)),
        ("disconnected_t", fixture_disconnected_t(), Int32(1), Int32(3)),
        ("empty_st", fixture_empty_st(), Int32(1), Int32(2)),
        ("n1", fixture_graph_n1_empty(), Int32(1), Int32(1)),
    ]
    for (name, g, s, t) in named
        @testset "$name" begin
            m = Int(g.m)
            w = [Float64(rand(rng, 0:7)) for _ in 1:m]
            net = fixture_network(g, s, t)
            z_or = oracle_sp_z(g, s, t, w)
            for slv in (COMB_SOLVER, AUTO_SOLVER)
                sol = solve_sp(net, w, slv)
                if z_or == Inf
                    test_infeasible_solution(sol, m)
                else
                    @test sol.status == ST_OK
                    @test sol.z == z_or
                    @test path_cost(sol.first, w) == sol.z
                    @test sol.method_used == METHOD_COMB
                end
            end
        end
    end
end

@testset "enumeration oracle on random DAGs" begin
    rng = MersenneTwister(20260909)
    for trial in 1:40
        n = rand(rng, 2:8)
        m = rand(rng, 0:14)
        g = rand_digraph(rng, n, m; dag = true)
        w = [Float64(rand(rng, 0:9)) for _ in 1:m]
        net = fixture_network(g, 1, n)
        sol = solve_sp(net, w, AUTO_SOLVER)
        z_or = oracle_sp_z(g, Int32(1), Int32(n), w)
        if z_or == Inf
            test_infeasible_solution(sol, m)
        else
            @test sol.status == ST_OK
            @test sol.z == z_or
            @test is_dag(g)
        end
    end
end

@testset "enumeration oracle on random cyclic digraphs" begin
    rng = MersenneTwister(20260909)
    for trial in 1:40
        n = rand(rng, 2:5)
        m = rand(rng, 0:8)
        g = rand_digraph(rng, n, m; dag = false)
        w = [Float64(rand(rng, 0:6)) for _ in 1:m]
        net = fixture_network(g, 1, n)
        sol = solve_sp(net, w, COMB_SOLVER)
        z_or = oracle_sp_z(g, Int32(1), Int32(n), w)
        if z_or == Inf
            test_infeasible_solution(sol, m)
        else
            @test sol.status == ST_OK
            @test sol.z == z_or
        end
    end
end

@testset "named single_path under first-stage costs" begin
    net, _ = load_instance("single_path")
    sol = solve_sp(net, net.costs.C, COMB_SOLVER)
    test_ok_single_stage(sol; z = 4.0, seq = Int32[1, 2, 3, 4])
end

@testset "named single_arc_paths under nominal second-stage costs" begin
    net, _ = load_instance("single_arc_paths")
    sol = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    test_ok_single_stage(sol; z = 1.0, seq = Int32[5])
end

@testset "named two_beads under nominal second-stage costs" begin
    net, _ = load_instance("two_beads")
    sol = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    test_ok_single_stage(sol; z = 20.0, seq = Int32[2, 3])
end

@testset "MIP SP matches COMB on fixtures" begin
    g = fixture_unique_path(4)
    net = fixture_network(g, 1, 4)
    w = [2.0, 5.0, 1.0]
    comb = solve_sp(net, w, COMB_SOLVER)
    mip = solve_sp(net, w, MIP_SOLVER)
    @test mip.status == ST_OK
    @test approxz(mip.z, comb.z)
    @test mip.method_used == METHOD_MIP
    @test mip.first.seq == comb.first.seq
end

@testset "MIP SP infeasible when t is unreachable" begin
    net = fixture_network(fixture_disconnected_t(), 1, 3)
    sol = solve_sp(net, [1.0], MIP_SOLVER)
    test_infeasible_solution(sol, 1)
    @test sol.method_used == METHOD_MIP
end

@testset "heap Dijkstra agrees with a quadratic scan, ties included" begin
    function scan_reference(g, s, t, w)
        n = Int(g.n)
        dist = fill(Inf, n)
        pred = zeros(Int32, n)
        done = falses(n)
        dist[s] = 0.0
        for _ in 1:n
            best, u = Inf, 0
            for v in 1:n
                if !done[v] && dist[v] < best
                    best, u = dist[v], v
                end
            end
            (u == 0 || best == Inf) && break
            done[u] = true
            u == t && break
            for a in outgoing(g, Int32(u))
                vv = Int(g.head[a])
                nd = dist[u] + w[a]
                if nd < dist[vv]
                    dist[vv], pred[vv] = nd, a
                end
            end
        end
        return dist, pred
    end
    for seed in 1:100
        rng = MersenneTwister(seed)
        n = rand(rng, 5:30)
        m = rand(rng, n:(4n))
        g = graph_new(n, rand(rng, 1:n, m), rand(rng, 1:n, m))
        w = Float64.(rand(rng, 0:3, m))        # small integers: many ties
        status, path, z = Rrsp._shortest_path_dijkstra(g, Int32(1), Int32(n), w)
        dist, pred = scan_reference(g, 1, n, w)
        if isinf(dist[n])
            @test status == ST_INFEASIBLE
        else
            expected = Int32[]
            v = n
            while v != 1
                pushfirst!(expected, pred[v])
                v = Int(g.tail[pred[v]])
            end
            @test status == ST_OK
            @test z == dist[n]
            @test path.seq == expected
        end
    end
end
