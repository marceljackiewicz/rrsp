@testset "costs_new length mismatch throws" begin
    @test_throws ArgumentError costs_new([1.0], [1.0, 2.0], [0.0])
    @test_throws ArgumentError costs_new(Float64[], [1.0], [0.0])
    @test_throws ArgumentError costs_new([1.0], [1.0], Float64[])
end

@testset "negative deviation throws" begin
    @test_throws ArgumentError costs_new([1.0], [1.0], [-0.1])
    @test_throws ArgumentError costs_new([0.0, 1.0], [0.0, 1.0], [0.0, -1.0])
end

@testset "zero deviation is admissible" begin
    c = costs_new([1.0, 2.0], [3.0, 4.0], [0.0, 0.0])
    @test c.d == [0.0, 0.0]
    @test c.C == [1.0, 2.0]
    @test c.c_hat == [3.0, 4.0]
end

@testset "empty cost vectors match m = 0" begin
    c = costs_new(Float64[], Float64[], Float64[])
    @test isempty(c.C)
end

@testset "costs_new copies the arrays" begin
    C = [1.0]
    c_hat = [2.0]
    d = [3.0]
    c = costs_new(C, c_hat, d)
    C[1] = 9.0
    c_hat[1] = 9.0
    d[1] = 9.0
    @test c.C[1] == 1.0
    @test c.c_hat[1] == 2.0
    @test c.d[1] == 3.0
end

@testset "network_new rejects terminals out of range" begin
    g = fixture_unique_path(3)
    c = uniform_costs(g)
    @test_throws ArgumentError network_new(g, 0, 3, c)
    @test_throws ArgumentError network_new(g, 1, 4, c)
    @test_throws ArgumentError network_new(g, 1, 0, c)
end

@testset "network_new rejects a cost length other than m" begin
    g = fixture_unique_path(3)
    c = costs_new([1.0], [1.0], [0.0])
    @test_throws ArgumentError network_new(g, 1, 3, c)
end

@testset "network with matching costs" begin
    g = fixture_two_beads()
    c = costs_new([0.0, 0.0, 0.0, 0.0], [50.0, 0.0, 20.0, 40.0], [10.0, 100.0, 30.0, 0.0])
    net = network_new(g, 1, 3, c)
    @test net.s == 1
    @test net.t == 3
    @test net.graph.n == 3
    @test net.costs.c_hat[1] == 50.0
    test_star_invariants(net.graph)
end

@testset "network_new copies costs" begin
    g = fixture_empty_st()
    c = costs_new(Float64[], Float64[], Float64[])
    net = network_new(g, 1, 2, c)
    # Mutating a fresh costs object used as input: rebuild with non-empty.
    g2 = fixture_unique_path(2)
    C = [5.0]
    costs = costs_new(C, [1.0], [0.0])
    net2 = network_new(g2, 1, 2, costs)
    costs.C[1] = 0.0
    @test net2.costs.C[1] == 5.0
end

@testset "s equals t is a valid network" begin
    g = fixture_unique_path(3)
    net = network_new(g, 2, 2, uniform_costs(g))
    @test net.s == net.t == 2
end
