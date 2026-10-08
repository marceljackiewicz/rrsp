# Path enumeration and the cut-generation engine (solve_adv_cuts, solve_rrsp_enum)
# against the independent oracles of setup.jl.

@testset "enumerate_st_paths matches the DFS oracle" begin
    for (name, g) in all_fixture_graphs()
        n = Int(g.n)
        n >= 2 || continue
        got = enumerate_st_paths(g, 1, n)
        want = enumerate_simple_st_paths(g, 1, n)
        @test sort([Vector(p.seq) for p in got]) == sort(want)
        for p in got
            @test length(p.chi) == Int(g.m)
            @test sum(p.chi) == length(p.seq)
        end
    end
    g = fixture_layered(3, 2)
    @test length(enumerate_st_paths(g, 1, 4)) == 8
    @test length(enumerate_st_paths(g, 1, 1)) == 1
    @test isempty(enumerate_st_paths(g, 1, 1)[1].seq)
    @test_throws ArgumentError enumerate_st_paths(g, 1, 4; limit = 7)
    @test_throws ArgumentError enumerate_st_paths(g, 1, 9)
end

@testset "solve_adv_cuts continuous matches the exact LP oracle" begin
    rng = MersenneTwister(404)
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO)
    bad = String[]
    n_checked = 0
    for trial in 1:8
        dag = iseven(trial)
        g = rand_digraph(rng, 5, 8; dag = dag)
        n = Int(g.n)
        paths = enumerate_simple_st_paths(g, 1, n)
        isempty(paths) && continue
        m = Int(g.m)
        c = costs_new(rand(rng, 0:3, m), rand(rng, 0:4, m), rand(rng, 0:4, m))
        net = network_new(g, 1, n, c)
        for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF), k in 0:3, gamma in (1.0, 3.0)
            # Exclusion/symmetric difference recoveries on cyclic digraphs need the MIP, which AUTO uses.
            params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = gamma)
            for seq in first(paths, 3)
                x = path_from_seq(g, 1, n, seq)
                z_or = oracle_adv_cont_z(g, 1, n, x, c.c_hat, c.d, nb, k, gamma)
                sol = solve_adv_cuts(net, params, x, solver)
                n_checked += 1
                if !(sol.status == ST_OK && approxz(sol.z, z_or))
                    push!(bad, "trial=$trial nb=$nb k=$k Γ=$gamma x=$seq z=$(sol.z) oracle=$z_or st=$(sol.status)")
                end
            end
        end
    end
    @test n_checked > 100
    isempty(bad) || foreach(println, first(bad, 10))
    @test isempty(bad)
end

@testset "solve_adv_cuts discrete matches subset enumeration" begin
    rng = MersenneTwister(505)
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO)
    bad = String[]
    n_checked = 0
    for trial in 1:6
        g = rand_digraph(rng, 5, 7; dag = true)
        n = Int(g.n)
        paths = enumerate_simple_st_paths(g, 1, n)
        isempty(paths) && continue
        m = Int(g.m)
        c = costs_new(rand(rng, 0:3, m), rand(rng, 0:4, m), rand(rng, 0:4, m))
        net = network_new(g, 1, n, c)
        for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF), k in (0, 1, 2), delta in (1, 2, 3)
            params = params_rrsp(U_DISC_BUDGET, nb, k; delta = delta)
            for seq in first(paths, 2)
                x = path_from_seq(g, 1, n, seq)
                z_or = oracle_adv_disc_z(g, 1, n, x, c.c_hat, c.d, nb, k, delta)
                sol = solve_adv_cuts(net, params, x, solver)
                n_checked += 1
                if !(sol.status == ST_OK && approxz(sol.z, z_or))
                    push!(bad, "trial=$trial nb=$nb k=$k Δ=$delta x=$seq z=$(sol.z) oracle=$z_or")
                end
            end
        end
    end
    @test n_checked > 100
    isempty(bad) || foreach(println, first(bad, 10))
    @test isempty(bad)
end

@testset "solve_adv_cuts agrees with the compact LP" begin
    net = adv_diamond()
    x = path_from_seq(net.graph, net.s, net.t, Int32[1, 2])
    for (nb, k) in ((NB_INCLUSION, 1), (NB_EXCLUSION, 1), (NB_SYMDIFF, 2))
        params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = 5.0)
        a = solve_adv(net, params, x, TEST_SOLVER)
        b = solve_adv_cuts(net, params, x, TEST_SOLVER)
        @test a.status == ST_OK && b.status == ST_OK
        @test approxz(a.z, b.z)
    end
    # Other uncertainty sets and a missing optimizer are not implemented here.
    @test solve_adv_cuts(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), x, TEST_SOLVER).status == ST_NOT_IMPL
    @test solve_adv_cuts(net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 1.0), x, AUTO_SOLVER).status == ST_NOT_IMPL
end

@testset "solve_rrsp_enum matches the exact oracles" begin
    rng = MersenneTwister(606)
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO)
    bad = String[]
    n_checked = 0
    for trial in 1:6
        g = rand_digraph(rng, 5, 8; dag = true)
        n = Int(g.n)
        isempty(enumerate_simple_st_paths(g, 1, n)) && continue
        m = Int(g.m)
        c = costs_new(rand(rng, 0:3, m), rand(rng, 0:4, m), rand(rng, 0:4, m))
        net = network_new(g, 1, n, c)
        for nb in (NB_INCLUSION, NB_EXCLUSION), k in (0, 1, 3)
            for gamma in (1.0, 3.0)
                params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = gamma)
                z_or = oracle_rrsp_cont_z(g, 1, n, c.C, c.c_hat, c.d, nb, k, gamma)
                sol = solve_rrsp_enum(net, params, solver)
                n_checked += 1
                (sol.status == ST_OK && approxz(sol.z, z_or) && approxz(sol.z, sol.z_first + sol.z_second)) ||
                    push!(bad, "cont trial=$trial nb=$nb k=$k Γ=$gamma z=$(sol.z) oracle=$z_or")
            end
            for delta in (1, 2)
                params = params_rrsp(U_DISC_BUDGET, nb, k; delta = delta)
                z_or = oracle_rrsp_disc_z(g, 1, n, c.C, c.c_hat, c.d, nb, k, delta)
                sol = solve_rrsp_enum(net, params, solver)
                n_checked += 1
                (sol.status == ST_OK && approxz(sol.z, z_or)) ||
                    push!(bad, "disc trial=$trial nb=$nb k=$k Δ=$delta z=$(sol.z) oracle=$z_or")
            end
        end
    end
    @test n_checked > 80
    isempty(bad) || foreach(println, first(bad, 10))
    @test isempty(bad)
end

@testset "solve_rrsp_enum options and edge cases" begin
    g = fixture_two_paths()
    net = network_new(g, 1, 4, costs_new([1.0, 0.0, 0.0, 2.0], [1.0, 1.0, 2.0, 2.0], [3.0, 3.0, 1.0, 1.0]))
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO)
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 2.0)
    full = solve_rrsp_enum(net, params, solver)
    @test full.status == ST_OK
    # Restricting the candidates can only make the value larger.
    ps = enumerate_st_paths(g, 1, 4)
    one = solve_rrsp_enum(net, params, solver; paths = ps[1:1])
    @test one.z >= full.z - 1e-9
    # Supplying the lower bounds does not change the result.
    lb = enum_lower_bounds(net, params, ps, solver)
    @test solve_rrsp_enum(net, params, solver; lower_bounds = lb).z ≈ full.z
    @test_throws ArgumentError solve_rrsp_enum(net, params, solver; lower_bounds = [1.0])
    # Interval uncertainty falls back to Rec.
    pint = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    @test solve_rrsp_enum(net, pint, solver).z ≈ solve_rec(net, pint, solver).z
    # No optimizer, no path.
    @test solve_rrsp_enum(net, params, AUTO_SOLVER).status == ST_NOT_IMPL
    g0 = fixture_disconnected_t()
    net0 = network_new(g0, 1, 3, uniform_costs(g0))
    @test solve_rrsp_enum(net0, params, solver).status == ST_INFEASIBLE
end
