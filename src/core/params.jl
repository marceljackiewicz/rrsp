"""
    params_nominal() -> Params

Parameters of the deterministic (nominal) problem: no uncertainty, no recovery.
"""
params_nominal() = Params(U_NOMINAL, NB_INCLUSION, Int32(0), 0.0, Int32(0))

"""
    params_rob(U; gamma=0.0, delta=0) -> Params

Parameters of the classical robust problem on uncertainty set `U`: recovery
size ``k = 0``. First-stage cost is omitted from the robust objective by the
solver, not by this constructor.
"""
function params_rob(U::Uncertainty; gamma::Real = 0.0, delta::Integer = 0)
    _check_budgets(gamma, delta)
    return Params(U, NB_INCLUSION, Int32(0), Float64(gamma), Int32(delta))
end

"""
    params_rrsp(U, nb, k; gamma=0.0, delta=0) -> Params

Parameters of a recoverable-robust problem.

# Throws
- `ArgumentError`: if `k < 0`, or if `gamma` is negative or not finite,
  or `delta < 0`.
"""
function params_rrsp(
    U::Uncertainty,
    nb::Neighborhood,
    k::Integer;
    gamma::Real = 0.0,
    delta::Integer = 0,
)
    k < 0 && throw(ArgumentError("k < 0"))
    _check_budgets(gamma, delta)
    return Params(U, nb, Int32(k), Float64(gamma), Int32(delta))
end

function _check_budgets(gamma::Real, delta::Integer)
    (isfinite(gamma) && gamma >= 0) || throw(ArgumentError("gamma must be finite and nonnegative"))
    delta >= 0 || throw(ArgumentError("delta < 0"))
    return nothing
end
