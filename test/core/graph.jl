@testset "graph construction" begin
    @testset "n < 1 throws" begin
        @test_throws ArgumentError graph_new(0, Int32[], Int32[])
        @test_throws ArgumentError graph_new(-1, Int32[], Int32[])
    end

    @testset "tail and head length mismatch throws" begin
        @test_throws ArgumentError graph_new(2, Int32[1], Int32[])
        @test_throws ArgumentError graph_new(2, Int32[1, 1], Int32[2])
    end

    @testset "endpoint out of range throws" begin
        @test_throws ArgumentError graph_new(2, Int32[0], Int32[2])
        @test_throws ArgumentError graph_new(2, Int32[1], Int32[0])
        @test_throws ArgumentError graph_new(2, Int32[1], Int32[3])
        @test_throws ArgumentError graph_new(2, Int32[3], Int32[1])
        @test_throws ArgumentError graph_new(2, Int32[-1], Int32[1])
    end

    @testset "label of wrong length throws" begin
        @test_throws ArgumentError graph_new(2, Int32[1], Int32[2]; label = Int32[7])
        @test_throws ArgumentError graph_new(2, Int32[1], Int32[2]; label = Int32[7, 8, 9])
    end

    @testset "m = 0 is valid" begin
        g = fixture_empty_st()
        @test g.n == 2
        @test g.m == 0
        @test g.first_out[1] == 1
        @test g.first_out[Int(g.n) + 1] == 1
        @test g.first_in[Int(g.n) + 1] == 1
        @test isempty(g.out_arc)
        @test isempty(g.in_arc)
        test_star_invariants(g)
        g1 = fixture_graph_n1_empty()
        @test g1.n == 1
        @test g1.m == 0
        test_star_invariants(g1)
    end

    @testset "isolated vertices have empty stars" begin
        g = fixture_isolated_s()
        test_star_invariants(g)
        @test out_degree(g, Int32(1)) == 0
        @test in_degree(g, Int32(1)) == 0
        @test isempty(collect(outgoing(g, Int32(1))))
        @test isempty(collect(incoming(g, Int32(4))))
        @test collect(outgoing(g, Int32(2))) == Int32[1]
        @test collect(incoming(g, Int32(3))) == Int32[1]
    end

    @testset "parallel arcs keep distinct indices" begin
        g = fixture_parallel_arcs()
        test_star_invariants(g)
        @test g.m == 5
        outs = collect(outgoing(g, Int32(1)))
        @test outs == Int32[1, 2, 3, 4, 5]
        @test all(a -> g.tail[a] == 1 && g.head[a] == 2, outs)
        @test collect(incoming(g, Int32(2))) == Int32[1, 2, 3, 4, 5]
        @test length(unique(outs)) == 5
    end

    @testset "self-loop appears in both stars" begin
        g = fixture_self_loop()
        test_star_invariants(g)
        # Arc 2 is 2 → 2.
        @test g.tail[2] == 2 && g.head[2] == 2
        @test 2 in collect(outgoing(g, Int32(2)))
        @test 2 in collect(incoming(g, Int32(2)))
    end

    @testset "graph_new copies endpoint arrays and labels" begin
        tail = Int32[1]
        head = Int32[2]
        label = Int32[10, 20]
        g = graph_new(2, tail, head; label = label)
        tail[1] = 2
        head[1] = 1
        label[1] = 0
        @test g.tail[1] == 1
        @test g.head[1] == 2
        @test g.label[1] == 10
        @test g.label[2] == 20
    end

    @testset "default labels are 1:n" begin
        g = graph_new(3, Int32[1], Int32[2])
        @test g.label == Int32[1, 2, 3]
    end

    @testset "vertex out of range on degree accessors throws" begin
        g = fixture_unique_path(3)
        @test_throws ArgumentError out_degree(g, Int32(0))
        @test_throws ArgumentError out_degree(g, Int32(4))
        @test_throws ArgumentError in_degree(g, 0)
        @test_throws ArgumentError outgoing(g, Int32(0))
        @test_throws ArgumentError incoming(g, Int32(99))
    end
end

@testset "star invariants on every fixture" begin
    for (name, g) in all_fixture_graphs()
        @testset "$name" begin
            test_star_invariants(g)
        end
    end
end

@testset "star invariants on random digraphs" begin
    rng = MersenneTwister(20260908)
    @testset "general" begin
        for trial in 1:40
            n = rand(rng, 1:12)
            m = rand(rng, 0:30)
            g = rand_digraph(rng, n, m; dag = false)
            test_star_invariants(g)
        end
    end
    @testset "acyclic" begin
        for trial in 1:30
            n = rand(rng, 2:12)
            m = rand(rng, 0:25)
            g = rand_digraph(rng, n, m; dag = true)
            test_star_invariants(g)
            @test is_dag(g)
        end
    end
end

@testset "is_dag" begin
    for (name, g) in all_fixture_graphs()
        @testset "$name" begin
            @test is_dag(g) == fixture_is_dag(name)
        end
    end

    @testset "self-loop on a single vertex" begin
        g = graph_new(1, Int32[1], Int32[1])
        test_star_invariants(g)
        @test !is_dag(g)
    end
end

@testset "type stability of star accessors" begin
    g = fixture_unique_path(4)
    @test @inferred(out_degree(g, Int32(1))) isa Int32
    @test @inferred(in_degree(g, Int32(2))) isa Int32
end
