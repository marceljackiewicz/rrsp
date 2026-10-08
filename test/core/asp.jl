function test_asp_tree(g::Graph, s::Integer, t::Integer, tree::AspTree)
    m = Int(g.m)
    @test Int(tree.n_nodes) == 2 * m - 1
    @test length(tree.nodes) == Int(tree.n_nodes)
    @test 1 <= Int(tree.root) <= Int(tree.n_nodes)
    seen_arc = zeros(UInt8, m)
    @inbounds for i in 1:Int(tree.n_nodes)
        nd = tree.nodes[i]
        if nd.op == ASP_LEAF
            @test nd.left == 0 && nd.right == 0
            a = Int(nd.arc)
            @test 1 <= a <= m
            @test seen_arc[a] == 0x00
            seen_arc[a] = 0x01
            @test nd.s == g.tail[a]
            @test nd.t == g.head[a]
        else
            @test nd.arc == 0
            @test 1 <= Int(nd.left) < i
            @test 1 <= Int(nd.right) < i
            L = tree.nodes[nd.left]
            R = tree.nodes[nd.right]
            if nd.op == ASP_SERIES
                @test L.t == R.s
                @test nd.s == L.s
                @test nd.t == R.t
            else
                @test nd.op == ASP_PARALLEL
                @test L.s == R.s == nd.s
                @test L.t == R.t == nd.t
            end
        end
    end
    @test all(==(0x01), seen_arc)
    rt = tree.nodes[tree.root]
    @test rt.s == Int32(s)
    @test rt.t == Int32(t)
    return nothing
end

@testset "single arc" begin
    g = graph_new(2, Int32[1], Int32[2])
    tree = asp_decompose(g, 1, 2)
    @test tree !== nothing
    test_asp_tree(g, 1, 2, tree)
    @test tree.nodes[tree.root].op == ASP_LEAF
end

@testset "two arcs in series" begin
    g = fixture_unique_path(3)
    tree = asp_decompose(g, 1, 3)
    @test tree !== nothing
    test_asp_tree(g, 1, 3, tree)
    @test tree.nodes[tree.root].op == ASP_SERIES
end

@testset "two parallel arcs" begin
    g = graph_new(2, Int32[1, 1], Int32[2, 2])
    tree = asp_decompose(g, 1, 2)
    @test tree !== nothing
    test_asp_tree(g, 1, 2, tree)
    @test tree.nodes[tree.root].op == ASP_PARALLEL
end

@testset "enum integer values" begin
    @test Int(ASP_LEAF) == 0
    @test Int(ASP_SERIES) == 1
    @test Int(ASP_PARALLEL) == 2
end

@testset "unique path, two_beads, diamond, parallel bundle" begin
    for (g, s, t) in (
        (fixture_unique_path(5), 1, 5),
        (fixture_two_beads(), 1, 3),
        (fixture_two_paths(), 1, 4),
        (fixture_parallel_arcs(), 1, 2),
        (fixture_diamond_plus_chord(), 1, 4),
        (fixture_layered(3, 2), 1, 4),
        (fixture_dag_not_layered(), 1, 3),
        (fixture_isolated_s(), 2, 3),
    )
        tail0 = copy(g.tail)
        head0 = copy(g.head)
        tree = asp_decompose(g, s, t)
        @test tree !== nothing
        test_asp_tree(g, s, t, tree)
        @test g.tail == tail0
        @test g.head == head0
    end
end

@testset "named two_beads and single_path" begin
    net, _ = load_instance("two_beads")
    tree = asp_decompose(net.graph, net.s, net.t)
    @test tree !== nothing
    test_asp_tree(net.graph, net.s, net.t, tree)
    net2, _ = load_instance("single_path")
    tree2 = asp_decompose(net2.graph, net2.s, net2.t)
    @test tree2 !== nothing
    test_asp_tree(net2.graph, net2.s, net2.t, tree2)
end

@testset "wheatstone is not ASP" begin
    g = fixture_wheatstone()
    @test is_dag(g)
    @test asp_decompose(g, 1, 4) === nothing
end

@testset "cyclic and disconnected are not ASP" begin
    @test asp_decompose(fixture_cycle_reachable(), 1, 4) === nothing
    @test asp_decompose(fixture_disconnected_t(), 1, 3) === nothing
    @test asp_decompose(fixture_empty_st(), 1, 2) === nothing
    @test asp_decompose(fixture_self_loop(), 1, 3) === nothing
end

@testset "s == t is not a two-terminal ASP graph" begin
    g = fixture_graph_n1_empty()
    @test asp_decompose(g, 1, 1) === nothing
end

@testset "isolated extra vertices still decompose" begin
    g = graph_new(4, Int32[1], Int32[2])
    tree = asp_decompose(g, 1, 2)
    @test tree !== nothing
    test_asp_tree(g, 1, 2, tree)
    @test tree.nodes[tree.root].op == ASP_LEAF
end
