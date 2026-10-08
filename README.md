# Rrsp

Julia library for the recoverable robust shortest path (RRSP) problem and its special cases
(shortest path, incremental, recoverable, robust, adversarial).

The library implements results from:

- M. Jackiewicz, A. Kasperski, and P. Zieliński,
  [Recoverable Robust Shortest Path Problem Under Interval Budgeted Uncertainty Representations](https://onlinelibrary.wiley.com/doi/abs/10.1002/net.22255),
  *Networks* 85 (1), 127–141.
- M. Jackiewicz, A. Kasperski, and P. Zieliński,
  [Computational Complexity of the Recoverable Robust Shortest Path Problem with Discrete Recourse](https://doi.org/10.1016/j.dam.2025.03.004),
  *Discrete Applied Mathematics* 370, 103–110.
- M. Jackiewicz, A. Kasperski and P. Zieliński,
  [Computational complexity of the recoverable robust shortest path problem in acyclic digraphs](https://arxiv.org/abs/2410.09425),
  arXiv:2410.09425.

## Funding

The authors were supported by the National Science Centre, Poland, grant 2022/45/B/HS4/00355.

## Requirements

- [Julia](https://julialang.org/) 1.11 or later
- [JuMP](https://jump.dev/JuMP.jl/stable/) 1.31 or later, for the mixed-integer formulations (installed automatically with the package)

A MIP optimizer is not required to install the package or to run combinatorial algorithms.
Tests use [HiGHS](https://highs.dev/) as a test-only dependency.

## Install

From the directory that **contains** this repository (so that `./Rrsp` is the package root):

```julia
import Pkg
Pkg.develop(path = "./Rrsp")
```

Then, from **inside** the package root (`Rrsp/`):

```bash
julia --project=.
```

```julia
using Rrsp
```

`Pkg.develop` makes the package available in the current Julia environment. `--project=.` uses this repository as its own environment (the default for development and tests).

## Documentation

The package documentation is generated with [Documenter.jl](https://documenter.juliadocs.org/stable/).
It is hosted at [GitHub Pages](https://marceljackiewicz.github.io/rrsp/).

The docs walk through the API (types, I/O, evaluation, solvers).
A quick usage example is in the [working example](https://marceljackiewicz.github.io/rrsp/#Working-example).

To generate the HTML locally, from the package root:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

You can open `docs/build/index.html` in a browser.

## Tests

Test-only dependencies (`Test`, `HiGHS`) live in `test/Project.toml`. They are not dependencies of the library.

From the package root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This loads `Rrsp` into an isolated test environment and runs `test/runtests.jl`.
This runs the unit tests only.

The regression tests of the functions used in the experiments for the thesis run the scripts of `experiments/` as subprocesses and takes longer.
Enable it, if needed, with the environment variable `RRSP_TEST_EXPERIMENTS=1` (it applies to that command only):
```bash
RRSP_TEST_EXPERIMENTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

This is a test, not the way to run the experiments. To reproduce the figures and tables of the thesis, see [Experiments](#experiments).

## Experiments

The figures and tables of the experiments chapter of the thesis (value of recovery under interval, continuous-budget and discrete-budget uncertainty)
are produced by recomputing the experiments, with Julia 1.11 or later.
In the `Rrsp` directory:
```bash
julia --project=experiments -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments experiments/run_all.jl             # recompute on small instances, about a minute
```

The step-by-step guide (recomputation, thesis scale, comparison with the thesis sources)
is in the documentation page [Reproducing the experiments](https://marceljackiewicz.github.io/rrsp/experiments/).

## Input and output conventions

### File format

```
s t
neighborhood=INC k=1 uncertainty=CONT gamma=100.0 delta=0
u v C c_hat d
...
```

- `s` and `t` are vertex **labels** (file identifiers).
- The parameter line is optional. If it is omitted, `params_nominal()` is used (nominal uncertainty, inclusion neighborhood, k = Γ = Δ = 0).
- Blank lines and comments (`# …`) are skipped.
- Each remaining line is one arc: tail label, head label, then first-stage cost `C`, nominal-cost `c_hat`, maximal cost deviation `d`.

Parameter keys and values:

| Key | Meaning | Values |
|---|---|---|
| `neighborhood` | path neighborhood type | `INC`, `EXC`, `SYM_DIFF` |
| `k` | neighborhood size parameter | nonnegative integer |
| `uncertainty` | uncertainty set type (or none) | `NOMINAL`, `INTERVAL`, `CONT`, `DISC` |
| `gamma` | continuous budget | floating-point; used when `uncertainty=CONT` |
| `delta` | discrete budget | nonnegative integer; used when `uncertainty=DISC` |

`CONT_BUDGET` and `DISC_BUDGET` are accepted as aliases of `CONT` and `DISC`. `SYMDIFF` is accepted as an alias of `SYM_DIFF`.

### Reading and writing

An instance is two objects: a `Network` (digraph, terminals, costs) and `Params` (uncertainty set, neighborhood, recovery size k, budgets Γ (continuous) and Δ (discrete)).
Files are read and written with `parse_rrsp` and `write_rrsp`.

```julia
# a small instance: two layers of two parallel arcs
write("two_beads.rrsp", """
1 3
neighborhood=INC k=1 uncertainty=CONT gamma=100.0 delta=0
1 2 0.0 50.0 10.0
1 2 0.0 0.0 100.0
2 3 0.0 20.0 30.0
2 3 0.0 40.0 0.0
""")

net, params = parse_rrsp("two_beads.rrsp")
write_rrsp("out.rrsp", net, params)
```

### Indexing

- Vertices are the dense (the arrays have no gaps) integers `1, …, n`; arcs are `1, …, m`.
- An arc is identified by its index, not by the pair of endpoints (parallel arcs are allowed).
- Arc `a` leaves `graph.tail[a]` and enters `graph.head[a]`.
- Identifiers given in instance files need not be dense. They are put into dense arrays upon deserialization.

### Costs

Each arc cost is a triple and is stored in `Costs`.
In every arc line of a file:

| Field | File column | Meaning |
|---|---|---|
| `C` | 3 | first-stage cost `C_a` |
| `c_hat` | 4 | nominal second-stage cost `ĉ_a` |
| `d` | 5 | maximum second-stage deviation `d_a ≥ 0` |

### Solving

Once input digraph is parsed, various functions can be used to solve different problems on it.
For example:
```julia
net, params = parse_rrsp("two_beads.rrsp")   # the file written in the example above
sp = solve_sp(net, net.costs.c_hat, Solver(; method = METHOD_COMB))
rob = solve_rob(net, params_rob(U_INTERVAL), Solver(; method = METHOD_COMB))
rec = solve_rec(net, params_rrsp(U_INTERVAL, NB_INCLUSION, 1), Solver(; method = METHOD_COMB))
eval_nominal(net, sp.first)   # ĉ(x)
eval_max(net, rob.first)      # (ĉ + d)(x)
eval_recovered(net, params_rrsp(U_NOMINAL, NB_INCLUSION, 1), rec.first, net.costs.c_hat, Solver(; method = METHOD_COMB))
```

The list of functions:
| Function | Problem | Notes |
|---|---|---|
| `solve_sp` | nonnegative shortest path | `METHOD_AUTO` / `METHOD_COMB`: combinatorial. `METHOD_MIP`: compact path formulation, needs `Solver.optimizer` (returns `ST_NOT_IMPL` if it is `nothing`). If `s = t`, the result is the trivial empty path of cost 0. Unreachable `t` is `ST_INFEASIBLE`. |
| `solve_rob` | robust shortest path (`k = 0`); first-stage cost is omitted from the robust objective | Interval uncertainty: the same algorithm under `ĉ+d`. Continuous budget: `min_P ĉ(P) + min(Γ,d(P))`. Discrete budget: Bertsimas–Sim algorithm (a polynomial number of shortest paths). |
| `solve_inc` | recovery from a given first-stage path | |
| `solve_adv` | worst-case recovered cost of a given first-stage path | |
| `solve_rec` | interval recoverable problem `min C(x)+(ĉ+d)(y)` with `y` in the neighborhood of `x` | Combinatorial on DAGs and on ASP graphs. |
| `solve_rrsp` | two-stage problem (RRSP) | Continuous budget: compact MIP (or an ASP MIP on inclusion). Discrete budget with `k > 0` has no compact formulation and returns `ST_NOT_IMPL`; use `solve_rrsp_enum` or `approx_rrsp`. |
| `solve_rrsp_enum` | RRSP, solved exactly by enumerating the commitment over the s–t paths | For digraphs with a moderate number of paths; needs an optimizer. |
| `approx_rrsp` | RRSP, `1/α` heuristic | DAGs only. |

**Choosing an optimizer.** The optimizer is an argument of `Solver`; there is no global default.
Install any [JuMP-compatible](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) MIP solver (for example HiGHS, Gurobi, CPLEX or SCIP) and pass its constructor:

```julia
using Rrsp, HiGHS
solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO, silent = true)
```

The combinatorial methods need no optimizer (`optimizer = nothing`, the default). `METHOD_MIP`, `solve_rrsp` with a continuous budget and `solve_rrsp_enum` do, and return `ST_NOT_IMPL` without one. The formulations are mixed-integer linear, so any solver with binary variables works; only HiGHS has been tested.

### Paths

A `Path` has a characteristic vector `chi` of length `m` (`chi[a] ∈ {0,1}`) and a sequence `seq` of the same arcs in traversal order from `s` to `t`.
The empty sequence is the trivial path when `s = t`, and the absence of a path otherwise.
