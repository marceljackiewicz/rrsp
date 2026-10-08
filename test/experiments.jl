@testset "sample_scenario stays inside U" begin
    g = gen_repsel([2, 2])
    net = network_new(g, 1, 3, overlay_alpha(g; alpha = 0.5, rng = MersenneTwister(7)))
    rng = MersenneTwister(8)
    w = sample_scenario(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), rng)
    @test length(w) == Int(g.m)
    @test all(a -> net.costs.c_hat[a] - 1e-12 <= w[a] <= net.costs.c_hat[a] + net.costs.d[a] + 1e-12, 1:Int(g.m))
    wd = sample_scenario(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 0; delta = 1), MersenneTwister(9))
    n_hi = count(a -> wd[a] > net.costs.c_hat[a] + 1e-12, 1:Int(g.m))
    @test n_hi <= 1
    wc = sample_scenario(net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 0; gamma = 3.0), MersenneTwister(10))
    dumped = sum(wc[a] - net.costs.c_hat[a] for a in 1:Int(g.m))
    @test dumped <= 3.0 + 1e-10
end

@testset "experiment_eval interval Rec row and CSV" begin
    g = gen_repsel([2, 2])
    net = network_new(g, 1, 3, overlay_alpha(g; alpha = 0.5, rng = MersenneTwister(11)))
    params = params_rrsp(U_INTERVAL, NB_INCLUSION, 1)
    row = experiment_eval(
        net, params, COMB_SOLVER;
        problem = :rrsp,
        family = "layered",
        cost_model = "alpha",
        c_mode = "zero",
        seed = 11,
        n_samples = 4,
        rng = MersenneTwister(12),
    )
    @test row.status == "OK"
    @test row.problem == "rrsp"
    @test row.n == 3
    @test row.m == 4
    @test isfinite(row.z)
    @test isfinite(row.z_nom)
    @test isfinite(row.z_wc)
    @test isfinite(row.z_sampled)
    @test row.n_samples == 4
    mktempdir() do dir
        path = joinpath(dir, "rows.csv")
        write_experiment_csv(path, [row])
        lines = readlines(path)
        @test lines[1] == experiment_csv_header()
        @test length(lines) == 2
        write_experiment_csv(path, [row]; append = true)
        @test length(readlines(path)) == 3
    end
end

@testset "experiment_sanity on a layered generator" begin
    g = gen_repsel([2, 2])
    net = network_new(g, 1, 3, overlay_alpha(g; alpha = 0.5, rng = MersenneTwister(13)))
    fails = experiment_sanity(net, COMB_SOLVER)
    @test isempty(fails)
end

@testset "nested models on interval C=0" begin
    g = gen_random_dag(5, 8; rng = MersenneTwister(15))
    net = network_new(g, 1, 5, with_first_stage(overlay_alpha(g; alpha = 0.5, rng = MersenneTwister(16)), :zero))
    pnom = params_rrsp(U_INTERVAL, NB_INCLUSION, 0)
    nom = experiment_eval(net, pnom, COMB_SOLVER; problem = :nom)
    rob = experiment_eval(net, params_rob(U_INTERVAL), COMB_SOLVER; problem = :rob)
    rrsp0 = experiment_eval(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 0), COMB_SOLVER; problem = :rrsp)
    rrsp1 = experiment_eval(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), COMB_SOLVER; problem = :rrsp)
    @test nom.status == "OK"
    @test rob.status == "OK"
    @test rrsp0.status == "OK"
    @test rrsp1.status == "OK"
    @test isapprox(rrsp0.z, rob.z; atol = 1e-8)
    @test rrsp1.z <= rob.z + 1e-8
    @test nom.z <= rob.z + 1e-8
    @test isfinite(nom.z_wc)
end

@testset "experiment_eval ROB discrete and approx" begin
    g = gen_repsel([2, 2])
    net = network_new(g, 1, 3, overlay_alpha(g; alpha = 0.5, rng = MersenneTwister(14)))
    rob = experiment_eval(
        net, params_rob(U_DISC_BUDGET; delta = 1), COMB_SOLVER;
        problem = :rob, family = "repsel",
    )
    @test rob.status == "OK"
    @test rob.k == 0
    ap = experiment_eval(
        net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 1; delta = 1), COMB_SOLVER;
        problem = :approx, family = "repsel",
    )
    @test ap.status == "OK"
    @test isfinite(ap.z)
end
