# Reproducing the experiments

This page takes you from `git clone` to the figures and tables of the
experiments chapter of the thesis, by recomputing the experiments. Two levels,
the second a superset of the first:

1. **A quick recomputation** (about a minute). Every experiment is run on
   reduced digraphs and the plots and tables are produced from the fresh data.
2. **The thesis-scale recomputation** (long-running; the discrete-budget
   experiment is by far the slowest).
   Run in slices on several cores if you need it.

## What is reproduced

| Content | Experiment | Result data | Generated file |
|---|---|---|---|
| value of recovery, independent costs | interval | `data/interval/replicates.csv` | `figures/fig_scale_vor.tex` |
| the same, as a table | interval | `data/interval/replicates.csv` | `tables/tab_interval_vor.tex` |
| digraph sizes | interval | `data/interval/replicates.csv` | `tables/tab_interval_instances.tex` |
| designated-path costs | interval | `data/interval/bottleneck.csv` | `figures/fig_scale_vor_bottleneck.tex` |
| the same, as a table | interval | `data/interval/bottleneck.csv` | `tables/tab_bottleneck_vor.tex` |
| continuous budget | continuous budget | `data/cont/all.csv` | `figures/fig_scale_vor_cont.tex` |
| tables of the continuous-budget experiment | continuous budget | `data/cont/all.csv` | `tables/tab_cont_vor.tex`, `tables/tab_cont_instances.tex` |
| discrete budget | discrete budget | `data/disc/all.csv` | `figures/fig_scale_vor_disc.tex` |
| tables of the discrete-budget experiment | discrete budget | `data/disc/all.csv` | `tables/tab_disc_vor.tex`, `tables/tab_disc_instances.tex` |

All paths are relative to `experiments/`. Each figure also gets an `.svg`
preview next to its `.tex`. The value of recovery of an instance is
``\mathrm{VoR}(k) = 100\,(Z(0) - Z(k))/Z(0)`` for the recovery budget ``k`` of
the inclusion neighborhood, where ``Z`` is the optimal value of the
recoverable robust problem (interval uncertainty, continuous budget or
discrete budget, according to the experiment).

## Step 1: get the code

You need [Git](https://git-scm.com/) and [Julia](https://julialang.org/downloads/)
1.11 or newer. LaTeX is optional, only for the PDF
preview.

```bash
git clone https://github.com/marceljackiewicz/rrsp.git
cd rrsp
```

## Step 2: install the dependencies

The experiments have their own environment, `experiments/`, which contains
this package, [HiGHS](https://highs.dev/) and [JuMP](https://jump.dev/). The
first call downloads and precompiles them (a few minutes).

```bash
julia --project=experiments -e 'using Pkg; Pkg.instantiate()'
```

All commands below are run from the repository root.

## Step 3: recompute everything on small instances

```bash
julia --project=experiments experiments/run_all.jl --scale=quick --pdf
```

This runs the interval, continuous-budget and discrete-budget experiments, then builds the figures and tables from the
fresh results. It takes about a minute. Nothing tracked is overwritten: raw
results go to `experiments/out/data/` and the TeX/SVG files to
`experiments/out/artifacts/`:

```
experiments/out/artifacts/
├── figures/   fig_scale_vor.tex  fig_scale_vor_bottleneck.tex
│              fig_scale_vor_cont.tex  fig_scale_vor_disc.tex  (and *.svg)
└── tables/    tab_{interval,bottleneck,cont,disc}_vor.tex
               tab_{interval,cont,disc}_instances.tex
```

`--pdf` compiles `preview.pdf` with `pdflatex`, with the thesis macros
(`\NeighborhoodSize`, `\STPath`, ...) replaced by plain stand-ins.

The files are the `figure` and `table` environments of the chapter. To use
them in the thesis, `\input` them in a document that defines the thesis
macros and loads `tikz`, `xcolor` and `amsmath`.

The `quick` scale uses fewer cost draws (2 per fixed digraph, 2 random DAGs)
and, for the continuous-budget and discrete-budget experiments, smaller digraphs,
so the curves differ from the thesis. The interval experiment is the exception: its digraphs are cheap, so `quick` only
reduces the number of draws.

## Step 4: the thesis scale

`--scale=full` uses the instances of the chapter:

| | Interval | Continuous budget | Discrete budget |
|---|---|---|---|
| Layered | height 12, width 4 | height 12, width 2 | height 8, width 2 |
| Grid | 7×7 | 7×7 | 5×5 |
| Arc-series-parallel | 12 gadgets | 12 gadgets | 8 gadgets |
| Random DAG | 10 digraphs × 3 cost draws; layered, arcs skipping layers with probability 0.5 | same, width 2 | same, height 8, width 2 |
| Fixed families | 30 cost draws on one digraph | same | same |
| Commitment | exact (`solve_rec`) | compact MIP (`solve_rrsp`) on arc-series-parallel digraphs, otherwise best of all ``s``-``t`` paths | best of all ``s``-``t`` paths |
| Relative running time | shortest (combinatorial algorithm) | much longer than the interval experiment (compact MIP on arc-series-parallel digraphs; elsewhere exhaustive search over paths, one adversary MIP per candidate, because the general compact MIP is far slower) | much longer than the continuous-budget experiment (exhaustive search over paths; no compact MIP exists for ``k > 0``) |

Seeds are fixed (`experiments/lib/instances.jl`), so every run builds the same
instances.

```bash
julia --project=experiments experiments/interval.jl --scale=full        # shortest
julia --project=experiments experiments/cont.jl --scale=full        # long
julia --project=experiments experiments/disc.jl --scale=full        # longest
julia --project=experiments experiments/make_artifacts.jl --data=experiments/out/data
```

The continuous-budget and discrete-budget experiments write one file per instance,
`experiments/out/data/cont/<family>/g<graph>_c<cost>.csv`, and skip files that
already exist, so an interrupted run resumes where it stopped (`--force`
recomputes). Slices of a run are independent processes; for example, one per
family and half of the cost draws:

```bash
for f in layered grid asp; do
  julia --project=experiments experiments/cont.jl --scale=full --family=$f --from=0 --to=14 &
  julia --project=experiments experiments/cont.jl --scale=full --family=$f --from=15 --to=29 &
done
julia --project=experiments experiments/cont.jl --scale=full --family=random_dag --from=0 --to=9 &
wait
```

For the fixed families `--from/--to` select cost draws `0:29`, for
`random_dag` they select digraphs `0:9`. The HiGHS thread pool is
process-global, so run several single-threaded processes rather than one
multi-threaded one.

The exact search of the continuous-budget and discrete-budget experiments enumerates the commitments
([`solve_rrsp_enum`](@ref)). Its cost depends on the number of paths and on the
budget; the discrete-budget runs with a large ``Δ`` on the 5×5
grid are the slowest. Instead of waiting for a full run you can inspect a
single point with the library:

```julia
using Rrsp, HiGHS, Random
g = gen_grid(5, 5)
net = network_new(g, 1, 25, with_first_stage(overlay_uniform(g; rng = MersenneTwister(1)), :zero))
paths = enumerate_st_paths(g, 1, 25)
solver = Solver(; optimizer = HiGHS.Optimizer, silent = true)
sol = solve_rrsp_enum(net, params_rrsp(U_DISC_BUDGET, NB_INCLUSION, 2; delta = 3), solver; paths = paths)
```

### Replacing the stored results

After a complete `--scale=full` run:

```bash
julia --project=experiments experiments/promote_data.jl
julia --project=experiments experiments/make_artifacts.jl
```

`promote_data.jl` merges the per-instance files into `experiments/data/` and
refuses to do so unless every family has all of its instances. Where the
stored continuous-budget and discrete-budget data came from is described in `experiments/data/README.md`.

## Files

| Path | Role |
|---|---|
| `experiments/interval.jl`, `cont.jl`, `disc.jl` | the three experiments; options `--scale`, `--family`, `--from`, `--to`, `--force` |
| `experiments/run_all.jl` | all experiments, then the artifacts |
| `experiments/make_artifacts.jl` | CSV → figures and tables |
| `experiments/compare_thesis.jl` | diff against the thesis sources |
| `experiments/promote_data.jl` | fresh run → tracked `data/` |
| `experiments/data/` | stored result CSVs (tracked) |
| `experiments/out/` | fresh runs, figures and tables (not tracked) |
| `experiments/lib/` | instances, aggregation, TikZ and table writers |

Environment variables: `RRSP_DATA_DIR`, `RRSP_RUN_DIR` and `RRSP_ARTIFACT_DIR`
move the three directories, and `RRSP_SCALE` sets the default scale.
