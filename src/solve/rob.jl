"""
    solve_rob(net, params, solver) -> Solution

Classical robust shortest path on `net` under `params`.

Recovery size must be ``k = 0``. First-stage cost ``C`` is omitted from the
objective: `z` is the worst-case second-stage cost of the chosen path.

- `U_INTERVAL`: shortest path under ``\\hat{c}+d``.
- `U_NOMINAL`: shortest path under ``\\hat{c}``.
- `U_CONT_BUDGET`: ``\\min_P \\hat{c}(P) + \\min(Γ, d(P))``, by comparing the
  nominal shortest path (plus ``Γ``) with the interval shortest path.
- `U_DISC_BUDGET`: Bertsimas–Sim, a polynomial number of shortest paths
  (``\\hat{c}(P)`` plus the ``Δ`` largest deviations on ``P``). `METHOD_MIP`
  uses the compact Bertsimas–Sim formulation.

The path is returned in `first`. `z_first` is still ``C(x)`` for reporting;
`z_second` equals `z`. `METHOD_MIP` uses a mixed-integer formulation when
`solver.optimizer` is set.

# Throws
- `ArgumentError`: if `params.k != 0`, if `gamma < 0` under a
  continuous budget, or if `delta < 0` under a discrete budget.
"""
function solve_rob(net::Network, params::Params, solver::Solver)::Solution
    params.k == 0 || throw(ArgumentError("solve_rob requires k == 0"))
    t0 = time_ns()
    m = Int(net.graph.m)
    U = params.uncertainty
    if U == U_DISC_BUDGET
        delta = _require_delta(params.delta)
        return _solve_rob_disc(net, delta, solver, t0)
    end
    if U == U_CONT_BUDGET
        gamma = _require_gamma(params.gamma)
        return _solve_rob_cont(net, gamma, solver, t0)
    end
    w = if U == U_INTERVAL
        ww = Vector{Float64}(undef, m)
        c_hat = net.costs.c_hat
        d = net.costs.d
        @inbounds for a in 1:m
            ww[a] = c_hat[a] + d[a]
        end
        _check_weights(ww, m)
        ww
    else
        _copy_weights(net.costs.c_hat, m)
    end
    if solver.method == METHOD_MIP
        solver.optimizer === nothing && return _not_impl(m, METHOD_MIP, t0)
        return _solve_sp_mip(net, w, solver, t0; rob = true)
    end
    status, path, z = _shortest_path(net.graph, net.s, net.t, w)
    return _finish_sp(net, status, path, z, t0; rob = true)
end
