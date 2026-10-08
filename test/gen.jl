@testset "gen_repsel with equal widths" begin
    g = gen_repsel([2, 2, 2])
    test_star_invariants(g)
    @test g.n == 4
    @test g.m == 6
    @test is_dag(g)
    @test is_st_layered(g, 1, 4)
end

@testset "gen_random_dag" begin
    rng = MersenneTwister(1)
    g = gen_random_dag(7, 12; rng = rng)
    test_star_invariants(g)
    @test is_dag(g)
    net = network_new(g, 1, 7, overlay_uniform(g; rng = MersenneTwister(2)))
    sol = solve_sp(net, net.costs.c_hat, COMB_SOLVER)
    @test sol.status == ST_OK
end

@testset "gen_grid" begin
    g = gen_grid(3, 4)
    test_star_invariants(g)
    @test g.n == 12
    @test is_dag(g)
    @test is_st_layered(g, 1, 12)
    g1 = gen_grid(1, 1)
    @test g1.n == 1
    @test g1.m == 0
end

@testset "gen_repsel and gen_asp" begin
    g = gen_repsel([2, 3, 2])
    test_star_invariants(g)
    @test g.n == 4
    @test g.m == 7
    @test asp_decompose(g, 1, 4) !== nothing
    @test is_st_layered(g, 1, 4)
    rng = MersenneTwister(3)
    ga = gen_asp(8; rng = rng)
    test_star_invariants(ga)
    @test asp_decompose(ga, 1, 2) !== nothing
end

function _arc_count(g::Graph, s::Int, t::Int)::Int
    m = Int(g.m)
    w = ones(Float64, m)
    net = network_new(g, s, t, costs_new(zeros(m), w, zeros(m)))
    sol = solve_sp(net, w, COMB_SOLVER)
    @test sol.status == ST_OK
    return Int(round(sol.z))
end

@testset "wide layered, skips, and ASP detours" begin
    g = gen_wide_layered(12, 4)
    test_star_invariants(g)
    @test g.n == 46
    @test g.m == 168
    @test is_dag(g)
    @test is_st_layered(g, 1, 46)
    @test _arc_count(g, 1, 46) == 12
    z = gen_layered_skips(12, 4; p = 0.0, rng = MersenneTwister(1))
    @test z.n == 46
    @test z.m == 168
    @test is_st_layered(z, 1, 46)
    full = gen_layered_skips(12, 4; p = 1.0, rng = MersenneTwister(1))
    @test full.m == 168 + 152
    @test is_dag(full)
    @test !is_st_layered(full, 1, 46)
    @test _arc_count(full, 1, 46) == 6
    det = gen_asp_detours(12)
    test_star_invariants(det)
    @test det.n == 37
    @test det.m == 48
    @test asp_decompose(det, 1, 37) !== nothing
    @test _arc_count(det, 1, 37) == 12
    @test !is_st_layered(det, 1, 37)
end

@testset "cost overlays" begin
    g = gen_repsel([2, 2])
    rng = MersenneTwister(4)
    ca = overlay_alpha(g; alpha = 0.5, rng = rng)
    @test cost_structure_alpha(ca) !== nothing
    @test isapprox(cost_structure_alpha(ca), 0.5; atol = 1e-12)
    cz = with_first_stage(ca, :zero)
    @test all(==(0.0), cz.C)
    cc = with_first_stage(ca, :c_hat)
    @test cc.C == cc.c_hat
    zn = overlay_zero_nominal(g; n_zero = 1, rng = MersenneTwister(5))
    @test cost_structure_alpha(zn) === nothing
    cb = overlay_designated_bottleneck(g, [1, 3])
    @test count(==(100.0), cb.d) == 2
    @test count(==(1.0), cb.d) == 2
    @test cb.c_hat[1] == 1.0 && cb.c_hat[2] == 8.0
    @test cb.C == cb.c_hat
    @test_throws ArgumentError overlay_designated_bottleneck(g, [5])
end

@testset "overlay_asp_detours draws per alternative" begin
    H = 5
    c = overlay_asp_detours(H; rng = MersenneTwister(21))
    @test length(c.c_hat) == 4H
    @test c.C == c.c_hat
    for gix in 0:(H - 1)
        b = 4 * gix
        # The three detour arcs share one pair.
        @test c.c_hat[b + 2] == c.c_hat[b + 3] == c.c_hat[b + 4]
        @test c.d[b + 2] == c.d[b + 3] == c.d[b + 4]
        @test 1.0 <= c.c_hat[b + 1] < 10.0
        @test 0.0 <= c.d[b + 1] < 10.0
        @test 1.0 / 3 <= c.c_hat[b + 2] < 10.0 / 3
    end
    @test overlay_asp_detours(H; rng = MersenneTwister(21)).c_hat == c.c_hat
    # Same stream as overlay_uniform for a non-gadget graph.
    g = gen_repsel([2, 2, 2])
    u = overlay_uniform(g; rng = MersenneTwister(2))
    @test with_first_stage(u, :c_hat).C == u.c_hat
end
