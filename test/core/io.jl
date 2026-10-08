function write_temp_rrsp(dir::AbstractString, name::AbstractString, contents::AbstractString)
    path = joinpath(dir, name)
    open(path, "w") do io
        print(io, contents)
    end
    return path
end

function labeled_endpoints(net::Network)
    g = net.graph
    s_lab = g.label[net.s]
    t_lab = g.label[net.t]
    arcs = [
        (g.label[g.tail[a]], g.label[g.head[a]], net.costs.C[a], net.costs.c_hat[a], net.costs.d[a])
        for a in 1:Int(g.m)
    ]
    return s_lab, t_lab, arcs
end

@testset "exclusion and symmetric difference on the parameter line" begin
    mktempdir() do dir
        p_exc = write_temp_rrsp(
            dir,
            "exc.rrsp",
            """
            1 2
            neighborhood=EXC k=2 uncertainty=NOMINAL gamma=0.0 delta=0
            1 2 1.0 1.0 0.0
            """,
        )
        _, params_exc = parse_rrsp(p_exc)
        @test params_exc.neighborhood == NB_EXCLUSION
        @test params_exc.k == 2
        p_sym = write_temp_rrsp(
            dir,
            "sym.rrsp",
            """
            1 2
            neighborhood=SYM_DIFF k=4 uncertainty=CONT gamma=1.5 delta=0
            1 2 1.0 1.0 0.0
            """,
        )
        _, params_sym = parse_rrsp(p_sym)
        @test params_sym.neighborhood == NB_SYMDIFF
        @test params_sym.k == 4
        @test params_sym.gamma == 1.5
        @test params_sym.uncertainty == U_CONT_BUDGET
    end
end

@testset "file without a parameter line is nominal" begin
    mktempdir() do dir
        path = write_temp_rrsp(
            dir,
            "v2.rrsp",
            """
            1 2
            1 2 1.0 2.0 3.0
            """,
        )
        net, params = parse_rrsp(path)
        test_star_invariants(net.graph)
        @test net.graph.n == 2
        @test net.graph.m == 1
        @test params.uncertainty == U_NOMINAL
        @test params.k == 0
        @test params.gamma == 0.0
        @test params.delta == 0
        @test net.costs.C[1] == 1.0
        @test net.costs.c_hat[1] == 2.0
        @test net.costs.d[1] == 3.0
    end
end

@testset "parameter line" begin
    mktempdir() do dir
        path = write_temp_rrsp(
            dir,
            "v2p.rrsp",
            """
            1 2
            neighborhood=EXC k=3 uncertainty=INTERVAL gamma=0.0 delta=0
            1 2 0.0 1.0 0.0
            """,
        )
        net, params = parse_rrsp(path)
        @test params.neighborhood == NB_EXCLUSION
        @test params.k == 3
        @test params.uncertainty == U_INTERVAL
        @test net.s == 1
        @test net.t == 2
    end
end

@testset "comments and blank lines are skipped" begin
    mktempdir() do dir
        path = write_temp_rrsp(
            dir,
            "comments.rrsp",
            """
            # terminals follow

            1 3

            1 2 0.0 1.0 0.0
            # an arc
            2 3 0.0 1.0 0.0
            """,
        )
        net, params = parse_rrsp(path)
        test_star_invariants(net.graph)
        @test net.graph.m == 2
        @test params.uncertainty == U_NOMINAL
    end
end

@testset "non-dense vertex labels keep appearance order" begin
    mktempdir() do dir
        path = write_temp_rrsp(
            dir,
            "labels.rrsp",
            """
            7 9
            7 3 1.0 2.0 0.5
            3 9 3.0 4.0 0.0
            """,
        )
        net, _ = parse_rrsp(path)
        g = net.graph
        @test g.n == 3
        @test g.label == Int32[7, 9, 3]
        @test net.s == 1
        @test net.t == 2
        @test g.tail == Int32[1, 3]
        @test g.head == Int32[3, 2]
        @test g.label[g.tail[1]] == 7
        @test g.label[g.head[1]] == 3
        @test g.label[g.head[2]] == 9
    end
end

@testset "isolated terminals with no arcs" begin
    mktempdir() do dir
        path = write_temp_rrsp(dir, "empty.rrsp", "1 2\n")
        net, _ = parse_rrsp(path)
        test_star_invariants(net.graph)
        @test net.graph.n == 2
        @test net.graph.m == 0
        @test net.s == 1
        @test net.t == 2
    end
end

@testset "write_rrsp / parse_rrsp round-trip preserves labels, arcs, costs, and params" begin
    g = fixture_two_beads()
    costs = costs_new(
        [0.0, 0.0, 0.0, 0.0],
        [50.0, 0.0, 20.0, 40.0],
        [10.0, 100.0, 30.0, 0.0],
    )
    net = network_new(g, 1, 3, costs)
    params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, 1; gamma = 100.0)
    mktempdir() do dir
        path = joinpath(dir, "out.rrsp")
        write_rrsp(path, net, params)
        net2, params2 = parse_rrsp(path)
        test_star_invariants(net2.graph)
        @test labeled_endpoints(net2) == labeled_endpoints(net)
        @test params2.uncertainty == params.uncertainty
        @test params2.neighborhood == params.neighborhood
        @test params2.k == params.k
        @test params2.gamma == params.gamma
        @test params2.delta == params.delta
    end
end

@testset "round-trip of a network with non-default labels" begin
    g = graph_new(3, Int32[1, 2], Int32[2, 3]; label = Int32[7, 3, 9])
    net = network_new(g, 1, 3, uniform_costs(g; C = 1.0, c_hat = 2.0, d = 0.5))
    mktempdir() do dir
        path = joinpath(dir, "lab.rrsp")
        write_rrsp(path, net, params_nominal())
        net2, params2 = parse_rrsp(path)
        @test net2.graph.label == Int32[7, 9, 3] || net2.graph.label == Int32[7, 3, 9]
        s2, t2, arcs2 = labeled_endpoints(net2)
        s1, t1, arcs1 = labeled_endpoints(net)
        @test (s2, t2) == (s1, t1)
        @test arcs2 == arcs1
        @test params2.uncertainty == U_NOMINAL
    end
end

@testset "missing file throws" begin
    @test_throws ArgumentError parse_rrsp(joinpath(mktempdir(), "no_such_file.rrsp"))
end

@testset "truncated arc line throws" begin
    mktempdir() do dir
        path = write_temp_rrsp(dir, "trunc.rrsp", "1 2\n1 2 0.0 1.0\n")
        @test_throws ArgumentError parse_rrsp(path)
    end
end

@testset "unknown neighborhood token throws" begin
    mktempdir() do dir
        path = write_temp_rrsp(
            dir,
            "badnb.rrsp",
            "1 2\nneighborhood=FOO k=0 uncertainty=NOMINAL\n1 2 0.0 1.0 0.0\n",
        )
        @test_throws ArgumentError parse_rrsp(path)
    end
end

@testset "unknown uncertainty token throws" begin
    mktempdir() do dir
        path = write_temp_rrsp(
            dir,
            "badu.rrsp",
            "1 2\nneighborhood=INC k=0 uncertainty=ELLIPSOID\n1 2 0.0 1.0 0.0\n",
        )
        @test_throws ArgumentError parse_rrsp(path)
    end
end

@testset "empty file throws" begin
    mktempdir() do dir
        path = write_temp_rrsp(dir, "emptyfile.rrsp", "")
        @test_throws ArgumentError parse_rrsp(path)
    end
end

@testset "header with only one terminal throws" begin
    mktempdir() do dir
        path = write_temp_rrsp(dir, "one.rrsp", "1\n")
        @test_throws ArgumentError parse_rrsp(path)
    end
end
