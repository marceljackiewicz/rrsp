# Shared fixtures and invariant checks for Phase A.
# Included from runtests.jl; not a testset of its own.

# ---------------------------------------------------------------------------
# Small named instances, in the `.rrsp` file format
# ---------------------------------------------------------------------------

const INSTANCES = Dict{String,String}(
    "single_path" => """
1 5
neighborhood=INC k=1 uncertainty=CONT gamma=3.0 delta=0
1 2 1.0 0.0 2.0
2 3 1.0 0.0 2.0
3 4 1.0 0.0 2.0
4 5 1.0 0.0 2.0
""",
    "single_arc_paths" => """
1 2
neighborhood=INC k=0 uncertainty=CONT gamma=100.0 delta=0
1 2 1.0 5.0 1.0
1 2 2.0 4.0 2.0
1 2 1.0 3.0 3.0
1 2 2.0 2.0 4.0
1 2 2.0 1.0 5.0
""",
    "two_beads" => """
1 3
neighborhood=INC k=1 uncertainty=CONT gamma=100.0 delta=0
1 2 0.0 50.0 10.0
1 2 0.0 0.0 100.0
2 3 0.0 20.0 30.0
2 3 0.0 40.0 0.0
""",
    "fixed_theta_counterexample" => """
1 4
neighborhood=INC k=1 uncertainty=CONT gamma=20.0 delta=0
1 2 0.0 20.0 5.0
2 4 0.0 20.0 5.0
1 3 0.0 10.0 40.0
3 4 0.0 10.0 40.0
1 4 0.0 0.0 300.0
""",
    "lambda_theta_counterexample" => """
1 5
neighborhood=INC k=2 uncertainty=CONT gamma=1.0 delta=0
1 2 0.0 1.0 100.0
2 3 0.0 1.0 100.0
3 4 0.0 1.0 100.0
4 5 0.0 1.0 100.0
1 6 100.0 0.0 100.0
6 2 100.0 0.0 100.0
2 7 100.0 0.0 100.0
7 3 100.0 0.0 100.0
3 8 100.0 0.0 100.0
8 4 100.0 0.0 100.0
4 9 100.0 0.0 100.0
9 5 100.0 0.0 100.0
""",
    "icalp_counterexample" => """
1 4
neighborhood=INC k=5 uncertainty=CONT gamma=1.0 delta=0
1 2 0.0 1.0 100.0
2 4 0.0 1.0 100.0
2 3 0.0 1.0 1000.0
1 3 0.0 1.0 100.0
3 4 0.0 1.0 1000.0
""",
)

const INSTANCE_NAMES = [
    "single_path",
    "single_arc_paths",
    "two_beads",
    "fixed_theta_counterexample",
    "lambda_theta_counterexample",
    "icalp_counterexample",
]

"""Parse the named instance of `INSTANCES` (through a temporary file)."""
function load_instance(name::AbstractString)
    mktempdir() do dir
        path = joinpath(dir, "$name.rrsp")
        write(path, INSTANCES[name])
        return parse_rrsp(path)
    end
end

# ---------------------------------------------------------------------------
# Graph fixtures
# ---------------------------------------------------------------------------

function fixture_graph_n1_empty()
    return graph_new(1, Int32[], Int32[])
end

function fixture_empty_st()
    return graph_new(2, Int32[], Int32[])
end

function fixture_isolated_s()
    # Vertices 1 and 4 are isolated; the only arc is 2 → 3.
    return graph_new(4, Int32[2], Int32[3])
end

function fixture_self_loop()
    # Path 1 → 2 → 3 together with a self-loop at 2.
    return graph_new(3, Int32[1, 2, 2], Int32[2, 2, 3])
end

function fixture_parallel_arcs()
    # Five parallel arcs 1 → 2, as in the `single_arc_paths` instance.
    return graph_new(2, Int32[1, 1, 1, 1, 1], Int32[2, 2, 2, 2, 2])
end

function fixture_unique_path(n::Int = 5)
    n >= 2 || throw(ArgumentError("unique path requires n >= 2"))
    tail = Int32[i for i in 1:(n - 1)]
    head = Int32[i for i in 2:n]
    return graph_new(n, tail, head)
end

function fixture_two_paths()
    # Diamond: 1 → 2 → 4 and 1 → 3 → 4.
    return graph_new(4, Int32[1, 2, 1, 3], Int32[2, 4, 3, 4])
end

function fixture_diamond_plus_chord()
    # Geometry of the `fixed_theta_counterexample` instance: two paths plus 1 → 4.
    return graph_new(4, Int32[1, 2, 1, 3, 1], Int32[2, 4, 3, 4, 4])
end

function fixture_two_beads()
    # Two layers of parallel arcs, as in the `two_beads` instance.
    return graph_new(3, Int32[1, 1, 2, 2], Int32[2, 2, 3, 3])
end

function fixture_layered(k::Int = 3, width::Int = 2)
    # Vertices 1, …, k+1 with `width` parallel arcs between i and i+1.
    # Every s–t path has k arcs.
    tail = Int32[]
    head = Int32[]
    for i in 1:k
        for _ in 1:width
            push!(tail, Int32(i))
            push!(head, Int32(i + 1))
        end
    end
    return graph_new(k + 1, tail, head)
end

function fixture_dag_not_layered()
    # Short path 1 → 3 and long path 1 → 2 → 3: arc counts differ.
    return graph_new(3, Int32[1, 1, 2], Int32[3, 2, 3])
end

function fixture_cycle_reachable()
    # 1 → 2 → 3 → 4 with a cycle 2 ⇄ 3 on every route to 4 except none: 3 → 2.
    return graph_new(4, Int32[1, 2, 3, 3], Int32[2, 3, 2, 4])
end

function fixture_cycle_irrelevant()
    # s–t arc 1 → 2 together with a disjoint 3 ⇄ 4 cycle.
    return graph_new(4, Int32[1, 3, 4], Int32[2, 4, 3])
end

function fixture_disconnected_t()
    # Arc 1 → 2; vertex 3 has no incident arc (used as an unreachable t).
    return graph_new(3, Int32[1], Int32[2])
end

function fixture_wheatstone()
    # Directed Wheatstone bridge: not two-terminal series-parallel.
    return graph_new(4, Int32[1, 1, 2, 3, 2], Int32[2, 3, 4, 4, 3])
end

function all_fixture_graphs()
    return [
        "n1_empty" => fixture_graph_n1_empty(),
        "empty_st" => fixture_empty_st(),
        "isolated_s" => fixture_isolated_s(),
        "self_loop" => fixture_self_loop(),
        "parallel_arcs" => fixture_parallel_arcs(),
        "unique_path" => fixture_unique_path(),
        "two_paths" => fixture_two_paths(),
        "diamond_plus_chord" => fixture_diamond_plus_chord(),
        "two_beads" => fixture_two_beads(),
        "layered" => fixture_layered(),
        "dag_not_layered" => fixture_dag_not_layered(),
        "cycle_reachable" => fixture_cycle_reachable(),
        "cycle_irrelevant" => fixture_cycle_irrelevant(),
        "disconnected_t" => fixture_disconnected_t(),
        "wheatstone" => fixture_wheatstone(),
    ]
end

function fixture_is_dag(name::AbstractString)::Bool
    cyclic = ("self_loop", "cycle_reachable", "cycle_irrelevant")
    return !(name in cyclic)
end

function rand_digraph(rng::AbstractRNG, n::Int, m::Int; dag::Bool = false)
    n >= 1 || throw(ArgumentError("n < 1"))
    tail = Vector{Int32}(undef, m)
    head = Vector{Int32}(undef, m)
    if dag
        n == 1 && m > 0 && throw(ArgumentError("no acyclic arcs on one vertex"))
        for a in 1:m
            u = rand(rng, 1:(n - 1))
            v = rand(rng, (u + 1):n)
            tail[a] = Int32(u)
            head[a] = Int32(v)
        end
    else
        for a in 1:m
            tail[a] = Int32(rand(rng, 1:n))
            head[a] = Int32(rand(rng, 1:n))
        end
    end
    return graph_new(n, tail, head)
end

# ---------------------------------------------------------------------------
# Star invariants (oracle 1): run on every graph we construct
# ---------------------------------------------------------------------------

function test_star_invariants(g::Graph)
    n = Int(g.n)
    m = Int(g.m)
    @test n >= 1
    @test m >= 0
    @test length(g.tail) == m
    @test length(g.head) == m
    @test length(g.out_arc) == m
    @test length(g.in_arc) == m
    @test length(g.first_out) == n + 1
    @test length(g.first_in) == n + 1
    @test length(g.label) == n
    @test g.first_out[1] == Int32(1)
    @test g.first_out[n + 1] == Int32(m + 1)
    @test g.first_in[1] == Int32(1)
    @test g.first_in[n + 1] == Int32(m + 1)
    for v in 1:n
        @test g.first_out[v] <= g.first_out[v + 1]
        @test g.first_in[v] <= g.first_in[v + 1]
    end
    sum_out = 0
    sum_in = 0
    seen_out = falses(max(m, 0))
    seen_in = falses(max(m, 0))
    for v in 1:n
        vv = Int32(v)
        od = out_degree(g, vv)
        idg = in_degree(g, vv)
        outs = collect(outgoing(g, vv))
        ins = collect(incoming(g, vv))
        @test od == Int32(length(outs))
        @test idg == Int32(length(ins))
        @test od == g.first_out[v + 1] - g.first_out[v]
        @test idg == g.first_in[v + 1] - g.first_in[v]
        sum_out += Int(od)
        sum_in += Int(idg)
        for a in outs
            @test Int32(1) <= a <= g.m
            @test g.tail[a] == vv
            @test !seen_out[Int(a)]
            seen_out[Int(a)] = true
        end
        for a in ins
            @test Int32(1) <= a <= g.m
            @test g.head[a] == vv
            @test !seen_in[Int(a)]
            seen_in[Int(a)] = true
        end
        for i in 2:length(outs)
            @test outs[i - 1] < outs[i]
        end
        for i in 2:length(ins)
            @test ins[i - 1] < ins[i]
        end
    end
    @test sum_out == m
    @test sum_in == m
    if m > 0
        @test all(seen_out)
        @test all(seen_in)
    end
    return nothing
end

function uniform_costs(g::Graph; C::Float64 = 1.0, c_hat::Float64 = 1.0, d::Float64 = 0.0)
    m = Int(g.m)
    return costs_new(fill(C, m), fill(c_hat, m), fill(d, m))
end

function fixture_network(
    g::Graph,
    s::Integer,
    t::Integer;
    C::Float64 = 1.0,
    c_hat::Float64 = 1.0,
    d::Float64 = 0.0,
)
    return network_new(g, s, t, uniform_costs(g; C = C, c_hat = c_hat, d = d))
end

const COMB_SOLVER = Solver(; method = METHOD_COMB)
const AUTO_SOLVER = Solver(; method = METHOD_AUTO)
const MIP_SOLVER_NO_OPT = Solver(; method = METHOD_MIP)

# ---------------------------------------------------------------------------
# Enumeration oracle (simple s–t paths) for nonnegative shortest path
# ---------------------------------------------------------------------------

function enumerate_simple_st_paths(g::Graph, s::Integer, t::Integer)::Vector{Vector{Int32}}
    ss = Int32(s)
    tt = Int32(t)
    n = Int(g.n)
    (1 <= Int(ss) <= n && 1 <= Int(tt) <= n) || throw(ArgumentError("vertex out of range"))
    found = Vector{Vector{Int32}}()
    if ss == tt
        push!(found, Int32[])
        return found
    end
    seq = Int32[]
    visited = falses(n)
    visited[Int(ss)] = true
    function dfs(v::Int32)
        if v == tt
            push!(found, copy(seq))
            return
        end
        for a in outgoing(g, v)
            w = g.head[a]
            if !visited[Int(w)]
                visited[Int(w)] = true
                push!(seq, a)
                dfs(w)
                pop!(seq)
                visited[Int(w)] = false
            end
        end
        return nothing
    end
    dfs(ss)
    return found
end

function oracle_sp_z(g::Graph, s::Integer, t::Integer, w::AbstractVector{<:Real})::Float64
    paths = enumerate_simple_st_paths(g, s, t)
    isempty(paths) && return Inf
    zmin = Inf
    for seq in paths
        p = path_from_seq(g, s, t, seq)
        z = path_cost(p, w)
        if z < zmin
            zmin = z
        end
    end
    return zmin
end

function test_infeasible_solution(sol::Solution, m::Integer)
    @test sol.status == ST_INFEASIBLE
    @test sol.z == Inf
    @test sol.z_first == Inf
    @test sol.z_second == Inf
    @test isempty(sol.first.seq)
    @test isempty(sol.second.seq)
    @test length(sol.first.chi) == Int(m)
    @test length(sol.second.chi) == Int(m)
    @test all(==(0x00), sol.first.chi)
    @test all(==(0x00), sol.second.chi)
    return nothing
end

function test_ok_single_stage(sol::Solution; z::Float64, seq = nothing)
    @test sol.status == ST_OK
    @test sol.z == z
    @test sol.method_used == METHOD_COMB
    @test sol.status != ST_NOT_IMPL
    @test isempty(sol.second.seq)
    @test all(==(0x00), sol.second.chi)
    @test sol.time_sec >= 0.0
    @test sol.mip_gap == 0.0
    if seq !== nothing
        @test sol.first.seq == seq
    end
    return nothing
end

const MIP_SOLVER = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_MIP, silent = true)
const TEST_SOLVER = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO, silent = true)

approxz(a::Float64, b::Float64) = isapprox(a, b; atol = 1e-8, rtol = 1e-8)

function in_neighborhood(x::Path, y::Path, nb::Neighborhood, k::Integer)::Bool
    length(x.chi) == length(y.chi) || throw(ArgumentError("path length mismatch"))
    incl = 0
    excl = 0
    for a in 1:length(x.chi)
        if y.chi[a] != 0x00 && x.chi[a] == 0x00
            incl += 1
        end
        if x.chi[a] != 0x00 && y.chi[a] == 0x00
            excl += 1
        end
    end
    nb == NB_INCLUSION && return incl <= k
    nb == NB_EXCLUSION && return excl <= k
    return (incl + excl) <= k
end

function oracle_inc_z(
    g::Graph,
    s::Integer,
    t::Integer,
    x::Path,
    w::AbstractVector{<:Real},
    nb::Neighborhood,
    k::Integer,
)::Float64
    zmin = Inf
    for seq in enumerate_simple_st_paths(g, s, t)
        y = path_from_seq(g, s, t, seq)
        if in_neighborhood(x, y, nb, k)
            z = path_cost(y, w)
            if z < zmin
                zmin = z
            end
        end
    end
    return zmin
end

function oracle_rec_z(
    g::Graph,
    s::Integer,
    t::Integer,
    C::AbstractVector{<:Real},
    cmax::AbstractVector{<:Real},
    nb::Neighborhood,
    k::Integer,
)::Float64
    zmin = Inf
    for xseq in enumerate_simple_st_paths(g, s, t)
        x = path_from_seq(g, s, t, xseq)
        for yseq in enumerate_simple_st_paths(g, s, t)
            y = path_from_seq(g, s, t, yseq)
            if in_neighborhood(x, y, nb, k)
                z = path_cost(x, C) + path_cost(y, cmax)
                if z < zmin
                    zmin = z
                end
            end
        end
    end
    return zmin
end

function cmax_of(net::Network)::Vector{Float64}
    m = Int(net.graph.m)
    w = Vector{Float64}(undef, m)
    @inbounds for a in 1:m
        w[a] = net.costs.c_hat[a] + net.costs.d[a]
    end
    return w
end

# Reference implementation of the k = 0 continuous adversary (dump Γ onto x).
function dump_gamma_on_path(
    x::Path,
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    gamma::Real,
)::Float64
    rem = Float64(gamma)
    z = 0.0
    arcs = Int[]
    @inbounds for a in 1:length(x.chi)
        x.chi[a] != 0x00 && push!(arcs, a)
    end
    sort!(arcs; by = a -> Float64(d[a]), rev = true)
    @inbounds for a in arcs
        da = Float64(d[a])
        take = da < rem ? da : rem
        rem -= take
        z += Float64(c_hat[a]) + take
    end
    return z
end

function oracle_rob_cont_z(
    g::Graph,
    s::Integer,
    t::Integer,
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    gamma::Real,
)::Float64
    zmin = Inf
    for seq in enumerate_simple_st_paths(g, s, t)
        x = path_from_seq(g, s, t, seq)
        z = dump_gamma_on_path(x, c_hat, d, gamma)
        z < zmin && (zmin = z)
    end
    return zmin
end

# Exact continuous adversary over an explicit set of recoveries Y:
#   max τ  s.t.  τ <= (ĉ + δ)·y  for all y in Y,  0 <= δ <= d,  sum δ <= Γ.
# (An earlier greedy "water-filling" version was only a lower bound.)
function _adv_lp(
    Y::Vector{Path},
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    gamma::Real,
)::Float64
    isempty(Y) && return Inf
    m = length(c_hat)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(model)
    delta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    tau = JuMP.@variable(model)
    for a in 1:m
        JuMP.set_upper_bound(delta[a], Float64(d[a]))
    end
    JuMP.@constraint(model, sum(delta) <= Float64(gamma))
    for y in Y
        JuMP.@constraint(
            model,
            tau <= sum(Float64(c_hat[a]) + delta[a] for a in 1:m if y.chi[a] != 0x00; init = 0.0),
        )
    end
    JuMP.@objective(model, Max, tau)
    JuMP.optimize!(model)
    JuMP.termination_status(model) == JuMP.MOI.OPTIMAL || error("oracle LP not optimal")
    return JuMP.objective_value(model)
end

function oracle_adv_cont_z(
    g::Graph,
    s::Integer,
    t::Integer,
    x::Path,
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    nb::Neighborhood,
    k::Integer,
    gamma::Real,
)::Float64
    Y = Path[]
    for seq in enumerate_simple_st_paths(g, s, t)
        y = path_from_seq(g, s, t, seq)
        in_neighborhood(x, y, nb, k) && push!(Y, y)
    end
    return _adv_lp(Y, c_hat, d, gamma)
end

function oracle_rrsp_cont_z(
    g::Graph,
    s::Integer,
    t::Integer,
    C::AbstractVector{<:Real},
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    nb::Neighborhood,
    k::Integer,
    gamma::Real,
)::Float64
    zmin = Inf
    for xseq in enumerate_simple_st_paths(g, s, t)
        x = path_from_seq(g, s, t, xseq)
        z = path_cost(x, C) + oracle_adv_cont_z(g, s, t, x, c_hat, d, nb, k, gamma)
        z < zmin && (zmin = z)
    end
    return zmin
end

function dump_delta_on_path(
    x::Path,
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    delta::Integer,
)::Float64
    z = 0.0
    devs = Float64[]
    @inbounds for a in 1:length(x.chi)
        x.chi[a] == 0x00 && continue
        z += Float64(c_hat[a])
        push!(devs, Float64(d[a]))
    end
    sort!(devs; rev = true)
    k = Int(delta)
    k > length(devs) && (k = length(devs))
    @inbounds for j in 1:k
        z += devs[j]
    end
    return z
end

function oracle_rob_disc_z(
    g::Graph,
    s::Integer,
    t::Integer,
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    delta::Integer,
)::Float64
    zmin = Inf
    for seq in enumerate_simple_st_paths(g, s, t)
        x = path_from_seq(g, s, t, seq)
        z = dump_delta_on_path(x, c_hat, d, delta)
        z < zmin && (zmin = z)
    end
    return zmin
end

function enumerate_subsets_bounded(m::Int, k::Int)::Vector{Vector{Int}}
    out = Vector{Vector{Int}}()
    buf = Int[]
    function rec(start::Int, left::Int)
        push!(out, copy(buf))
        left == 0 && return
        for i in start:m
            push!(buf, i)
            rec(i + 1, left - 1)
            pop!(buf)
        end
        return nothing
    end
    rec(1, k)
    return out
end

function oracle_adv_disc_z(
    g::Graph,
    s::Integer,
    t::Integer,
    x::Path,
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    nb::Neighborhood,
    k::Integer,
    delta::Integer,
)::Float64
    Y = Path[]
    for seq in enumerate_simple_st_paths(g, s, t)
        y = path_from_seq(g, s, t, seq)
        in_neighborhood(x, y, nb, k) && push!(Y, y)
    end
    isempty(Y) && return Inf
    m = length(c_hat)
    zmax = -Inf
    for S in enumerate_subsets_bounded(m, Int(delta))
        inS = falses(m)
        for a in S
            inS[a] = true
        end
        zmin = Inf
        for y in Y
            z = 0.0
            @inbounds for a in 1:m
                y.chi[a] == 0x00 && continue
                z += Float64(c_hat[a])
                inS[a] && (z += Float64(d[a]))
            end
            z < zmin && (zmin = z)
        end
        zmin > zmax && (zmax = zmin)
    end
    return zmax
end

function oracle_rrsp_disc_z(
    g::Graph,
    s::Integer,
    t::Integer,
    C::AbstractVector{<:Real},
    c_hat::AbstractVector{<:Real},
    d::AbstractVector{<:Real},
    nb::Neighborhood,
    k::Integer,
    delta::Integer,
)::Float64
    zmin = Inf
    for xseq in enumerate_simple_st_paths(g, s, t)
        x = path_from_seq(g, s, t, xseq)
        z = path_cost(x, C) + oracle_adv_disc_z(g, s, t, x, c_hat, d, nb, k, delta)
        z < zmin && (zmin = z)
    end
    return zmin
end
