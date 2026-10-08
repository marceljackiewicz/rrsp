@testset "eval_nominal is c_hat · χ" begin
    g = fixture_unique_path(4)
    c_hat = [2.0, 3.0, 5.0]
    costs = costs_new(fill(9.0, 3), c_hat, fill(1.0, 3))
    net = network_new(g, 1, 4, costs)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 3])
    @test eval_nominal(net, x) == 10.0
    @test eval_nominal(net, x) == path_cost(x, net.costs.c_hat)
    @test @inferred(eval_nominal(net, x)) isa Float64
end

@testset "eval_max is (c_hat + d) · χ" begin
    g = fixture_unique_path(4)
    c_hat = [2.0, 3.0, 5.0]
    d = [1.0, 0.0, 4.0]
    costs = costs_new(fill(9.0, 3), c_hat, d)
    net = network_new(g, 1, 4, costs)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 3])
    @test eval_max(net, x) == 15.0
    @test eval_max(net, x) == path_cost(x, net.costs.c_hat) + path_cost(x, net.costs.d)
    @test @inferred(eval_max(net, x)) isa Float64
end

@testset "eval of the empty path is 0" begin
    g = fixture_two_paths()
    net = fixture_network(g, 1, 4; c_hat = 5.0, d = 2.0)
    x = path_empty(Int(g.m))
    @test eval_nominal(net, x) == 0.0
    @test eval_max(net, x) == 0.0
end

@testset "eval on a proper subpath of two_beads" begin
    g = fixture_two_beads()
    costs = costs_new(
        [0.0, 0.0, 0.0, 0.0],
        [50.0, 0.0, 20.0, 40.0],
        [10.0, 100.0, 30.0, 0.0],
    )
    net = network_new(g, 1, 3, costs)
    y = path_from_seq(g, Int32(1), Int32(3), Int32[2, 3])
    @test eval_nominal(net, y) == 20.0
    @test eval_max(net, y) == 150.0
end

@testset "eval rejects a characteristic vector of the wrong length" begin
    net = fixture_network(fixture_unique_path(3), 1, 3)
    x = path_empty(5)
    @test_throws ArgumentError eval_nominal(net, x)
    @test_throws ArgumentError eval_max(net, x)
end

@testset "eval_max equals eval_nominal when d = 0" begin
    g = fixture_two_paths()
    net = fixture_network(g, 1, 4; c_hat = 3.0, d = 0.0)
    x = path_from_seq(g, Int32(1), Int32(4), Int32[1, 2])
    @test eval_max(net, x) == eval_nominal(net, x) == 6.0
end
