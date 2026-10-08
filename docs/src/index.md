# RRSP Solver Documentation

Julia library for the recoverable robust shortest path (RRSP) problem and its special cases
(shortest path, incremental, recoverable, robust, adversarial).
This document is the user API of the Robust Recoverable Shortest Path Solver.
The documentation adheres to the `0.1.0` version of the `Rrsp` package.

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

The authors were supported by the National Science Centre, Poland, grant 2022/45/B/HS4/00355.

## How to install

### Requirements

- [Julia](https://julialang.org/) 1.11 or later
- [JuMP](https://jump.dev/JuMP.jl/stable/) 1.31 or later, for the mixed-integer formulations (installed automatically with the package)

A MIP optimizer is not required to install the package or to run combinatorial algorithms.
Tests use [HiGHS](https://highs.dev/) as a test-only dependency.
The optimizer is set as an argument of each solving function: [`Solver`](@ref), see [Choosing an optimizer](@ref).

### Installation

Use `Pkg` to make the package available in a session.

1. Clone the repository.
2. From the directory that **contains** the package root (so that `./Rrsp` contains this repository):

```julia
import Pkg
Pkg.develop(path = "./Rrsp")
```

3. Load the package:

```julia
using Rrsp
```

From **inside** the package root, `julia --project=.` uses this repository as its environment (the default for development and tests).

### Tests

From the package root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This runs the unit tests only.
Test-only dependencies (`Test`, `HiGHS`) live in `test/Project.toml`. They are not dependencies of the library.

The regression tests of the functions used in the experiments for the thesis run the scripts of `experiments/` as subprocesses and takes longer.
Enable it, if needed, with the environment variable `RRSP_TEST_EXPERIMENTS=1` (it applies to that command only):
```bash
RRSP_TEST_EXPERIMENTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

This is a test, not the way to run the experiments. To reproduce the figures and tables of the thesis, see [Reproducing the experiments](experiments.md).

### Building this documentation locally

From the package root:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` in a browser.

## Indexing conventions

- Vertices are the dense (the arrays have no gaps) integers `1, …, n`; arcs are `1, …, m`.
- An arc is identified by its index, not by the pair of endpoints (parallel arcs are allowed).
  Arc `a` leaves `graph.tail[a]` and enters `graph.head[a]`.
- File identifiers need not be dense. They are stored in `graph.label`.
  Dense indices are assigned in first-appearance order: origin `s`, destination `t` if distinct, then arc endpoints in file order.

## Enumerations

```@docs
Neighborhood
Uncertainty
Rrsp.Method
Status
```

## Data types

```@docs
Graph
Costs
Network
Params
Path
Solution
Solver
```

[`AspOp`](@ref), [`AspNode`](@ref), and [`AspTree`](@ref) describe a binary
decomposition of a two-terminal arc-series-parallel (ASP) digraph. Incremental
and recoverable combinatorial solvers call [`asp_decompose`](@ref) and, when a
tree exists, run an ``O(k m)`` DP on it instead of the DAG dynamic network.

```@docs
AspOp
AspNode
AspTree
```

## Construction

```@docs
graph_new
costs_new
network_new
params_nominal
params_rob
params_rrsp
path_empty
path_from_seq
path_from_chi
solution_empty
```

The keyword constructor of [`Solver`](@ref) is documented on that type.

## Graph queries

```@docs
out_degree
in_degree
outgoing
incoming
is_dag
asp_decompose
```

## Generators and experiment utilities

Structured topologies and cost overlays for the experiment chapter.

```@docs
gen_wide_layered
gen_layered_skips
gen_random_dag
gen_grid
gen_repsel
gen_asp
gen_asp_detours
is_st_layered
layer_vertex_ids
overlay_uniform
overlay_alpha
overlay_zero_nominal
overlay_designated_bottleneck
overlay_asp_detours
with_first_stage
sample_scenario
experiment_eval
write_experiment_csv
experiment_csv_header
experiment_sanity
```

The experiments of the thesis chapter (interval, continuous-budget and discrete-budget uncertainty) are built from these
generators; how to run them, from `git clone` to the plots, is on the page
[Reproducing the experiments](experiments.md).

## Paths and evaluation

[`path_cost`](@ref) is the inner product of a weight vector with ``χ``.
[`eval_nominal`](@ref) and [`eval_max`](@ref) evaluate a committed path;
they do not choose a first-stage path. [`eval_recovered`](@ref) is the
incremental cost of recovering from a committed path.
[`eval_worstcase`](@ref) is the adversarial cost of a committed path.

```@docs
path_cost
path_length
eval_nominal
eval_max
eval_recovered
eval_worstcase
```

## Parsing input

Parsing is separate from solving. One [`Network`](@ref) can be reused with
many [`Params`](@ref) objects (for example a sweep on ``k``).

```@docs
parse_rrsp
write_rrsp
```

`parse_rrsp` and `write_rrsp` use the same format.

```
s t
neighborhood=INC k=1 uncertainty=CONT gamma=100.0 delta=0
u v C c_hat d
...
```

- `s` and `t` are vertex **labels** (file identifiers).
- The parameter line is optional. If it is omitted, [`params_nominal`](@ref) is
  used (nominal uncertainty, inclusion neighborhood, ``k = Γ = Δ = 0``).
- Blank lines and comments (`# …`) are skipped.
- Each remaining line is one arc: tail label, head label, then first-stage
  cost `C`, nominal second-stage cost `c_hat` (``\hat{c}``), and deviation `d`
  with ``d \ge 0``.

| Key | Meaning | Values |
|---|---|---|
| `neighborhood` | path neighborhood type | `INC`, `EXC`, `SYM_DIFF` |
| `k` | neighborhood size parameter ``k`` | nonnegative integer |
| `uncertainty` | uncertainty set type (or none) | `NOMINAL`, `INTERVAL`, `CONT`, `DISC` |
| `gamma` | continuous budget ``Γ`` | floating-point; used when `uncertainty=CONT` |
| `delta` | discrete budget ``Δ`` | nonnegative integer; used when `uncertainty=DISC` |

`CONT_BUDGET` and `DISC_BUDGET` are accepted as aliases of `CONT` and `DISC`.
`SYMDIFF` is accepted as an alias of `SYM_DIFF`.

## Solving

Every solving function takes a [`Solver`](@ref).
The path of a single-stage problem is returned in `first`; `second` is empty.
If `status != ST_OK`, both paths are empty and the cost fields are `Inf`.

| Function | Problem | Notes |
|---|---|---|
| [`solve_sp`](@ref) | nonnegative shortest path | `METHOD_AUTO` / `METHOD_COMB`: combinatorial. `METHOD_MIP`: compact path formulation, needs `Solver.optimizer` (returns `ST_NOT_IMPL` if it is `nothing`). If `s = t`, the result is the trivial empty path of cost 0. Unreachable `t` is `ST_INFEASIBLE`. |
| [`solve_rob`](@ref) | robust shortest path (`k = 0`); first-stage cost is omitted from the robust objective | Interval uncertainty: the same algorithm under `ĉ+d`. Continuous budget: `min_P ĉ(P) + min(Γ,d(P))`. Discrete budget: Bertsimas–Sim algorithm (a polynomial number of shortest paths; `METHOD_MIP` uses the compact formulation). |
| [`solve_inc`](@ref) | recovery from a given first-stage path | |
| [`solve_adv`](@ref) | worst-case recovered cost of a given first-stage path | |
| [`solve_rec`](@ref) | interval recoverable problem `min C(x)+(ĉ+d)(y)` with `y` in the neighborhood of `x` | Combinatorial on DAGs and on ASP graphs. |
| [`solve_rrsp`](@ref) | two-stage problem (RRSP) | Continuous budget: compact MIP (or an ASP MIP on inclusion). Discrete budget with `k > 0` has no compact formulation and returns `ST_NOT_IMPL`; use `solve_rrsp_enum` or `approx_rrsp`. |
| [`solve_rrsp_enum`](@ref) | RRSP, solved exactly by enumerating the commitment over the s–t paths | For digraphs with a moderate number of paths; needs an optimizer. |
| [`approx_rrsp`](@ref) | RRSP, `1/α` heuristic | DAGs only. |

### Choosing an optimizer

The optimizer is an argument of [`Solver`](@ref); there is no global default.
Install any [JuMP-compatible](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
MIP solver (for example HiGHS, Gurobi, CPLEX or SCIP) and pass its constructor:

```julia
using Rrsp, HiGHS
solver = Solver(; optimizer = HiGHS.Optimizer, method = METHOD_AUTO, silent = true)
```

The combinatorial methods need no optimizer (`optimizer = nothing`, the
default). `METHOD_MIP`, [`solve_rrsp`](@ref) with a continuous budget and
[`solve_rrsp_enum`](@ref) do, and return `ST_NOT_IMPL` without one. The
formulations are mixed-integer linear, so any solver with binary variables
works; only HiGHS has been tested.

### Shortest path

```@docs
solve_sp
```

### Classical robust shortest path

```@docs
solve_rob
```

### Incremental shortest path

```@docs
solve_inc
```

### Recoverable shortest path

```@docs
solve_rec
```

### Adversarial evaluation

```@docs
solve_adv
```

### Recoverable-robust shortest path

```@docs
solve_rrsp
```

### Enumeration of the commitment

When the number of ``s``–``t`` paths is moderate, the commitment ``x`` can be
chosen by enumeration. For every candidate the adversarial problem is solved
by cut generation: a master problem over the attack (an LP for a continuous
budget, a MIP for a discrete budget) and, as separation oracle,
[`solve_inc`](@ref). This covers the discrete budget at ``k > 0`` and serves as
an independent check of the compact formulations. All of it needs
`Solver.optimizer`. The continuous-budget and discrete-budget experiments use it.

```@docs
enumerate_st_paths
solve_adv_cuts
enum_lower_bounds
solve_rrsp_enum
```

### Approximation

```@docs
cost_structure_alpha
cost_structure_kappa
approx_rrsp
approx_bound_factors
```

## Working example

The example uses a small instance, written to `two_beads.rrsp` below, and the
combinatorial solver (no MIP optimizer).

Make sure `Rrsp` is visible. Use `import` to avoid name collisions; if that
is not an issue, `using` is fine.

```julia
julia> import Rrsp
```

Write the instance file (two layers of two parallel arcs):

```julia
julia> write("two_beads.rrsp", """
       1 3
       neighborhood=INC k=1 uncertainty=CONT gamma=100.0 delta=0
       1 2 0.0 50.0 10.0
       1 2 0.0 0.0 100.0
       2 3 0.0 20.0 30.0
       2 3 0.0 40.0 0.0
       """)
```

Parse the file. The result is a network and a parameter record; interval
robustness uses [`params_rob`](@ref), not the file's recovery size.

```julia
julia> net, params = Rrsp.parse_rrsp("two_beads.rrsp")

julia> solver = Rrsp.Solver(; method = Rrsp.METHOD_COMB)
```

Nominal shortest path under ``\hat{c}``, and classical robust shortest path
under interval uncertainty (shortest path under ``\hat{c}+d``):

```julia
julia> sp = Rrsp.solve_sp(net, net.costs.c_hat, solver)

julia> rob = Rrsp.solve_rob(net, Rrsp.params_rob(Rrsp.U_INTERVAL), solver)
```

Recoverable shortest path under interval uncertainty with inclusion neighborhood
size ``k = 1``:

```julia
julia> rec = Rrsp.solve_rec(net, Rrsp.params_rrsp(Rrsp.U_INTERVAL, Rrsp.NB_INCLUSION, 1), solver)
```

Check `sp.status == Rrsp.ST_OK`. The objective is `sp.z`; the path is
`sp.first.seq` (arc indices in traversal order) and `sp.first.chi`. Evaluate
a committed path without resolving:

```julia
julia> Rrsp.eval_nominal(net, sp.first)   # ĉ(x)

julia> Rrsp.eval_max(net, rob.first)      # (ĉ + d)(x)

julia> Rrsp.eval_recovered(net, Rrsp.params_rrsp(Rrsp.U_NOMINAL, Rrsp.NB_INCLUSION, 1), rec.first, net.costs.c_hat, solver)
```

The instance can be written back with [`write_rrsp`](@ref):

```julia
julia> Rrsp.write_rrsp("out.rrsp", net, params)
```

If `s = t`, [`solve_sp`](@ref) returns the trivial empty path of cost `0`.
If `t` is unreachable, the status is `ST_INFEASIBLE`.

## Index

```@index
```




