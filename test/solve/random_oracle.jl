# Seeded random differential test: continuous-budget ADV and RRSP against the
# exact LP oracle over enumerated simple s–t paths (see `_adv_lp` in setup.jl).
#
# Covers every neighborhood, k = 0..3, two budgets, DAGs and cyclic digraphs.
# The exclusion neighborhood with k ≥ |x| is included on purpose.

@testset "random continuous ADV vs exact LP oracle" begin
    rng = MersenneTwister(20261002)
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO)
    n_checked = 0
    bad = String[]
    for trial in 1:14
        dag = isodd(trial) ? false : true
        g = rand_digraph(rng, 5, 8; dag = dag)
        n = Int(g.n)
        paths = enumerate_simple_st_paths(g, 1, n)
        isempty(paths) && continue
        m = Int(g.m)
        c = costs_new(rand(rng, 0:3, m), rand(rng, 0:4, m), rand(rng, 0:4, m))
        net = network_new(g, 1, n, c)
        for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF), k in 0:3, gamma in (1.0, 3.0)
            params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = gamma)
            for seq in paths
                x = path_from_seq(g, 1, n, seq)
                z_or = oracle_adv_cont_z(g, 1, n, x, c.c_hat, c.d, nb, k, gamma)
                sol = solve_adv(net, params, x, solver)
                n_checked += 1
                if !(sol.status == ST_OK && approxz(sol.z, z_or))
                    push!(bad, "trial=$trial dag=$dag nb=$nb k=$k Γ=$gamma x=$seq z=$(sol.z) oracle=$z_or st=$(sol.status)")
                end
            end
        end
    end
    @test n_checked > 500
    isempty(bad) || foreach(println, first(bad, 10))
    @test isempty(bad)
end

@testset "random continuous RRSP vs exact LP oracle" begin
    rng = MersenneTwister(77)
    solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO)
    bad = String[]
    n_checked = 0
    for trial in 1:8
        dag = isodd(trial) ? false : true
        g = rand_digraph(rng, 5, 8; dag = dag)
        n = Int(g.n)
        isempty(enumerate_simple_st_paths(g, 1, n)) && continue
        m = Int(g.m)
        c = costs_new(rand(rng, 0:3, m), rand(rng, 0:4, m), rand(rng, 0:4, m))
        net = network_new(g, 1, n, c)
        for nb in (NB_INCLUSION, NB_EXCLUSION, NB_SYMDIFF), k in (0, 1, 3), gamma in (1.0, 3.0)
            params = params_rrsp(U_CONT_BUDGET, nb, k; gamma = gamma)
            z_or = oracle_rrsp_cont_z(g, 1, n, c.C, c.c_hat, c.d, nb, k, gamma)
            sol = solve_rrsp(net, params, solver)
            n_checked += 1
            if !(sol.status == ST_OK && approxz(sol.z, z_or))
                push!(bad, "trial=$trial dag=$dag nb=$nb k=$k Γ=$gamma z=$(sol.z) oracle=$z_or st=$(sol.status)")
            end
        end
    end
    @test n_checked > 80
    isempty(bad) || foreach(println, first(bad, 10))
    @test isempty(bad)
end
