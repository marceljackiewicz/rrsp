# Discrete-budget experiment: value of recovery under discrete budgeted uncertainty (RRSP, exact).
#
#   julia --project=experiments experiments/disc.jl [--scale=quick|full]
#         [--family=layered,grid,asp,random_dag] [--from=i --to=j] [--force]
#
# Zero first-stage costs, inclusion neighborhood, and smaller digraphs than
# the interval and continuous-budget experiments so that the s–t paths can be enumerated: height 8 and width 2
# for the layered and random-DAG families, a 5×5 grid, and 8 ASP gadgets.
# For each instance the budget Δ runs through 0, 1, ..., Δ_force (see
# lib/forcing.jl) and the recovery budget k through 0, ..., longest.
# The commitment is the best of all s–t paths (`solve_rrsp_enum`).
#
# Output: disc/<family>/g<graph>_c<cost>.csv below `run_root()`.

@isdefined(run_cont) || include(joinpath(@__DIR__, "cont.jl"))

const DISC_HEADER = ["n", "m", "ell", "longest", "npaths", "delta_force", "delta", "k", "z", "z_nom", "z_int"]

function solve_disc_instance(inst::Instance; tol = 1e-4)
    g = inst.g
    solver = enum_solver()
    costs = instance_costs(inst; first_stage = :zero)
    net = network_new(g, inst.s, inst.t, costs)
    paths = enumerate_st_paths(g, inst.s, inst.t)
    pc = path_costs(paths, costs.c_hat, costs.d)
    ell = shortest_arc_count(g, inst.s, inst.t)
    longest = maximum(length(p.seq) for p in paths)

    # Adversarial value when every path is a recovery option; it is the same
    # for all commitments and bounds the value of every neighborhood.
    floors = Dict{Int,Float64}()
    function floor_of(delta::Int)
        return get!(floors, delta) do
            full = params_rrsp(U_DISC_BUDGET, NB_INCLUSION, longest; delta = delta)
            sol = solve_adv_cuts(net, full, paths[1], solver)
            sol.status == ST_OK || error("floor failed: $(sol.status) (Δ = $delta)")
            sol.z
        end
    end
    d_force = discrete_forcing_budget(net, pc, floor_of; tol = tol)
    log_line("  n=$(g.n) m=$(g.m) paths=$(length(paths)) ell=$ell longest=$longest Δforce=$d_force")

    lbs = Dict{Int,Vector{Float64}}()
    lower_bounds(k) = get!(lbs, k) do
        enum_lower_bounds(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, k; delta = 1), paths, solver)
    end

    rows = Tuple[]
    meta = (Int(g.n), Int(g.m), ell, longest, length(paths), d_force)
    for delta in 0:d_force
        z_floor = floor_of(delta)
        z_prev = Inf
        for k in 0:longest
            sol = solve_rrsp_enum(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, k; delta = delta), solver;
                                  paths = paths, lower_bounds = lower_bounds(k))
            sol.status == ST_OK || error("RRSP failed: $(sol.status) (Δ=$delta k=$k)")
            z = sol.z
            z > z_prev + 1e-2 && error("value rose in k (Δ=$delta k=$k: $z_prev -> $z)")
            z_prev = z
            push!(rows, (meta..., delta, k, z, pc.z_nom, pc.z_int))
            if z <= z_floor + tol         # larger k stay at the floor
                for kk in (k + 1):longest
                    push!(rows, (meta..., delta, kk, z_floor, pc.z_nom, pc.z_int))
                end
                break
            end
        end
    end
    return rows
end

run_disc(opts; root = run_root()) = run_instances(:disc, "disc", DISC_HEADER, solve_disc_instance; opts = opts, root = root)

if abspath(PROGRAM_FILE) == @__FILE__
    run_disc(parse_options())
end
