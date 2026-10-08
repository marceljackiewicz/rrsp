@testset "empty path" begin
    p = path_empty(5)
    @test p.seq == Int32[]
    @test length(p.chi) == 5
    @test all(==(0x00), p.chi)
    @test path_cost(p, zeros(5)) == 0.0
    @test path_length(p) == 0
    @test path_cost(p, [1.0, 2.0, 3.0, 4.0, 5.0]) == 0.0
end

@testset "path_empty of m = 0" begin
    p = path_empty(0)
    @test isempty(p.chi)
    @test isempty(p.seq)
end

@testset "s equals t, empty sequence is the trivial path" begin
    g = fixture_unique_path(3)
    p = path_from_seq(g, Int32(2), Int32(2), Int32[])
    @test p.seq == Int32[]
    @test all(==(0x00), p.chi)
end

@testset "empty sequence with s distinct from t throws" begin
    g = fixture_unique_path(3)
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(3), Int32[])
end

@testset "unique path from sequence" begin
    g = fixture_unique_path(5)
    seq = Int32[1, 2, 3, 4]
    p = path_from_seq(g, Int32(1), Int32(5), seq)
    @test p.seq == seq
    @test path_length(p) == 4
    @test p.chi == UInt8[0x01, 0x01, 0x01, 0x01]
    @test path_cost(p, [10.0, 20.0, 30.0, 40.0]) == 100.0
end

@testset "path_from_seq copies the sequence" begin
    g = fixture_unique_path(3)
    seq = Int32[1, 2]
    p = path_from_seq(g, Int32(1), Int32(3), seq)
    seq[1] = 2
    @test p.seq[1] == 1
end

@testset "path_from_seq rejects a walk that does not start at s" begin
    g = fixture_unique_path(4)
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(4), Int32[2, 3])
end

@testset "path_from_seq rejects a skipped vertex" begin
    g = fixture_unique_path(4)
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(4), Int32[1, 3])
end

@testset "path_from_seq rejects a walk that does not end at t" begin
    g = fixture_unique_path(4)
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(4), Int32[1, 2])
end

@testset "path_from_seq rejects an arc id out of range" begin
    g = fixture_unique_path(3)
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(3), Int32[0, 1])
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(3), Int32[1, 3])
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(3), Int32[-1])
end

@testset "path_from_seq rejects a repeated arc" begin
    g = fixture_cycle_reachable()
    # Arcs: 1: 1→2, 2: 2→3, 3: 3→2, 4: 3→4. The cycle 2→3→2 repeats vertices.
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(4), Int32[1, 2, 3, 2, 4])
end

@testset "path_from_seq rejects a self-loop" begin
    g = fixture_self_loop()
    # Arcs: 1: 1→2, 2: 2→2, 3: 2→3.
    @test_throws ArgumentError path_from_seq(g, Int32(1), Int32(3), Int32[1, 2, 3])
end

@testset "path_from_seq rejects s equals t with a nonempty cycle" begin
    g = fixture_cycle_reachable()
    @test_throws ArgumentError path_from_seq(g, Int32(2), Int32(2), Int32[2, 3])
end

@testset "path_cost weight length mismatch throws" begin
    g = fixture_unique_path(3)
    p = path_from_seq(g, Int32(1), Int32(3), Int32[1, 2])
    @test_throws ArgumentError path_cost(p, [1.0])
    @test_throws ArgumentError path_cost(p, Float64[])
end

@testset "path_cost is linear in the characteristic vector" begin
    g = fixture_two_beads()
    p = path_from_seq(g, Int32(1), Int32(3), Int32[1, 4])
    c_hat = [10.0, 100.0, 20.0, 40.0]
    d = [1.0, 2.0, 3.0, 4.0]
    @test path_cost(p, c_hat) == 50.0
    @test path_cost(p, d) == 5.0
    @test path_cost(p, c_hat) + path_cost(p, d) == path_cost(p, c_hat .+ d)
end

@testset "two distinct simple paths on two beads" begin
    g = fixture_two_beads()
    p14 = path_from_seq(g, Int32(1), Int32(3), Int32[1, 4])
    p23 = path_from_seq(g, Int32(1), Int32(3), Int32[2, 3])
    @test p14.chi != p23.chi
    @test p14.seq != p23.seq
    @test sum(p14.chi) == 2
    @test sum(p23.chi) == 2
end

@testset "path_from_chi round-trip on a unique path" begin
    g = fixture_unique_path(5)
    p = path_from_seq(g, Int32(1), Int32(5), Int32[1, 2, 3, 4])
    q = path_from_chi(g, Int32(1), Int32(5), p.chi)
    @test q.seq == p.seq
    @test q.chi == p.chi
end

@testset "path_from_chi round-trip on two beads" begin
    g = fixture_two_beads()
    for seq in (Int32[1, 3], Int32[1, 4], Int32[2, 3], Int32[2, 4])
        p = path_from_seq(g, Int32(1), Int32(3), seq)
        q = path_from_chi(g, Int32(1), Int32(3), p.chi)
        @test q.seq == p.seq
        @test q.chi == p.chi
    end
end

@testset "path_from_chi copies the characteristic vector" begin
    g = fixture_unique_path(3)
    chi = UInt8[0x01, 0x01]
    p = path_from_chi(g, Int32(1), Int32(3), chi)
    chi[1] = 0x00
    @test p.chi[1] == 0x01
end

@testset "path_from_chi rejects a characteristic vector that is not a path" begin
    g = fixture_two_beads()
    # Both arcs leaving 1: branching, not a path.
    @test_throws ArgumentError path_from_chi(g, Int32(1), Int32(3), UInt8[1, 1, 1, 0])
    # No arcs.
    @test_throws ArgumentError path_from_chi(g, Int32(1), Int32(3), UInt8[0, 0, 0, 0])
    # Wrong length.
    @test_throws ArgumentError path_from_chi(g, Int32(1), Int32(3), UInt8[1, 0])
end

@testset "path_from_chi of the empty characteristic vector when s equals t" begin
    g = fixture_two_paths()
    p = path_from_chi(g, Int32(2), Int32(2), zeros(UInt8, Int(g.m)))
    @test isempty(p.seq)
end

@testset "type stability of path_cost" begin
    g = fixture_unique_path(3)
    p = path_from_seq(g, Int32(1), Int32(3), Int32[1, 2])
    w = [1.0, 2.0]
    @test @inferred(path_cost(p, w)) isa Float64
end
