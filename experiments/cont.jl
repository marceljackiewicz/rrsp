# Continuous-budget experiment: value of recovery under continuous budgeted uncertainty (RRSP, exact).
#
#   julia --project=experiments experiments/cont.jl [--scale=quick|full]
#         [--family=layered,grid,asp,random_dag] [--from=i --to=j] [--force]
#
# Zero first-stage costs, inclusion neighborhood. For each instance the
# budget Γ runs through the fractions 0, 0.2, ..., 1 of the forcing budget
# (see lib/forcing.jl) and the recovery budget k through 0, ..., longest.
# On arc-series-parallel digraphs each (Γ, k) is solved by `solve_rrsp`, i.e. by the
# compact MIP. On the other families the commitment is the best of all s–t
# paths (`solve_rrsp_enum`): the general compact MIP (Bold–Goerigk) does not
# finish within minutes on a single instance of these sizes.
#
# Output: cont/<family>/g<graph>_c<cost>.csv below `run_root()`.

using HiGHS

@isdefined(EXPERIMENTS_DIR) || include(joinpath(@__DIR__, "lib", "common.jl"))
@isdefined(EXPERIMENT_SPECS) || include(joinpath(@__DIR__, "lib", "instances.jl"))
@isdefined(run_instances) || include(joinpath(@__DIR__, "lib", "runner.jl"))
@isdefined(PathCosts) || include(joinpath(@__DIR__, "lib", "forcing.jl"))

const CONT_FRACTIONS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
const CONT_HEADER = [
    "n", "m", "ell", "longest", "npaths", "gamma_force", "gamma_sat", "gamma", "gamma_frac", "k", "z", "z_nom", "z_int",
]

enum_solver() = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO, silent = true)
mip_solver() = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_MIP, mip_gap = 1e-9, silent = true)

function solve_cont_instance(inst::Instance; fractions = CONT_FRACTIONS, tol = 1e-5)
    g = inst.g
    solver = enum_solver()
    compact = asp_decompose(g, inst.s, inst.t) !== nothing      # compact MIP applies
    mip = mip_solver()
    costs = instance_costs(inst; first_stage = :zero)
    net = network_new(g, inst.s, inst.t, costs)
    paths = enumerate_st_paths(g, inst.s, inst.t)
    pc = path_costs(paths, costs.c_hat, costs.d)
    ell = shortest_arc_count(g, inst.s, inst.t)
    longest = maximum(length(p.seq) for p in paths)
    g_sat = saturation_budget(pc)
    g_force = continuous_forcing_budget(paths, pc, costs.d, HiGHS.Optimizer)
    g_force + 1e-4 >= g_sat || error("forcing budget $g_force below saturation budget $g_sat")
    log_line("  n=$(g.n) m=$(g.m) paths=$(length(paths)) ell=$ell longest=$longest Γforce=$(round(g_force; digits = 3))")

    # Lower bounds C(x) + min_{y ∈ N(x)} ĉ(y) depend on k but not on Γ.
    lbs = Dict{Int,Vector{Float64}}()
    lower_bounds(k) = get!(lbs, k) do
        enum_lower_bounds(net, params_rrsp(U_CONT_BUDGET, NB_INCLUSION, k; gamma = 1.0), paths, solver)
    end

    rows = Tuple[]
    meta = (Int(g.n), Int(g.m), ell, longest, length(paths), g_force, g_sat)
    for frac in fractions
        gamma = frac * g_force
        full = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, longest; gamma = gamma)
        floor_sol = solve_adv_cuts(net, full, paths[1], solver)
        floor_sol.status == ST_OK || error("floor failed: $(floor_sol.status)")
        z_floor = floor_sol.z
        for k in 0:longest
            if k == longest
                z = z_floor
            else
                params = params_rrsp(U_CONT_BUDGET, NB_INCLUSION, k; gamma = gamma)
                sol = compact ? solve_rrsp(net, params, mip) :
                      solve_rrsp_enum(net, params, solver; paths = paths, lower_bounds = lower_bounds(k))
                sol.status == ST_OK || error("RRSP failed: $(sol.status) (frac=$frac k=$k)")
                z = sol.z
            end
            z + 1e-4 >= z_floor || error("value below the full-recovery value (frac=$frac k=$k)")
            z <= z_floor + tol && (z = z_floor)
            push!(rows, (meta..., gamma, frac, k, z, pc.z_nom, pc.z_int))
            if k < longest && z == z_floor        # larger k stay at the floor
                for kk in (k + 1):longest
                    push!(rows, (meta..., gamma, frac, kk, z_floor, pc.z_nom, pc.z_int))
                end
                break
            end
        end
    end
    return rows
end

run_cont(opts; root = run_root()) = run_instances(:cont, "cont", CONT_HEADER, solve_cont_instance; opts = opts, root = root)

if abspath(PROGRAM_FILE) == @__FILE__
    run_cont(parse_options())
end
