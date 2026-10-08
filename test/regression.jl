# Parse every named instance (see `INSTANCES` in setup.jl): topology and costs only (no solver).

@testset "each named instance parses" begin
    for name in INSTANCE_NAMES
        @testset "$name" begin
            net, params = load_instance(name)
            test_star_invariants(net.graph)
            @test 1 <= net.s <= net.graph.n
            @test 1 <= net.t <= net.graph.n
            @test length(net.costs.C) == Int(net.graph.m)
            @test params.k >= 0
            @test params.gamma >= 0.0
            @test params.delta >= 0
        end
    end
end

@testset "named instances round-trip through write_rrsp" begin
    mktempdir() do dir
        for name in INSTANCE_NAMES
            @testset "$name" begin
                net, params = load_instance(name)
                out = joinpath(dir, "$name.rrsp")
                write_rrsp(out, net, params)
                net2, params2 = parse_rrsp(out)
                test_star_invariants(net2.graph)
                @test labeled_endpoints(net2) == labeled_endpoints(net)
                @test params2.neighborhood == params.neighborhood
                @test params2.k == params.k
                @test params2.gamma == params.gamma
            end
        end
    end
end

@testset "single_path geometry" begin
    net, params = load_instance("single_path")
    @test net.graph.n == 5
    @test net.graph.m == 4
    @test is_dag(net.graph)
    @test params.neighborhood == NB_INCLUSION
    @test params.k == 1
    @test params.gamma == 3.0
    @test params.uncertainty == U_CONT_BUDGET
    p = path_from_seq(net.graph, net.s, net.t, Int32[1, 2, 3, 4])
    @test all(==(0x01), p.chi)
end

@testset "single_arc_paths geometry" begin
    net, params = load_instance("single_arc_paths")
    @test net.graph.n == 2
    @test net.graph.m == 5
    @test params.k == 0
    @test all(a -> net.graph.tail[a] == net.s && net.graph.head[a] == net.t, 1:Int(net.graph.m))
end

@testset "two_beads geometry" begin
    net, _ = load_instance("two_beads")
    @test net.graph.n == 3
    @test net.graph.m == 4
    @test out_degree(net.graph, net.s) == 2
    @test in_degree(net.graph, net.t) == 2
end
