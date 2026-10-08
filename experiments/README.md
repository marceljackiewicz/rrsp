# Experiments

The experiments chapter: value of recovery under interval, continuous-budget and discrete-budget uncertainty (`interval.jl`, `cont.jl`, `disc.jl`). Full guide, from `git clone` to the plots: [Reproducing the experiments](https://marceljackiewicz.github.io/rrsp/experiments/) (source: `docs/src/experiments.md`).

From the package root, after `julia --project=experiments -e 'using Pkg; Pkg.instantiate()'`:

```bash
julia --project=experiments experiments/run_all.jl [--scale=quick|full]  # recompute -> out/
julia --project=experiments experiments/compare_thesis.jl --thesis=DIR   # diff with the thesis sources
```

| Path | Role |
|---|---|
| `interval.jl`, `cont.jl`, `disc.jl` | the experiments (`--scale`, `--family`, `--from`, `--to`, `--force`) |
| `run_all.jl` | all experiments, then the artifacts, below `out/` |
| `make_artifacts.jl` | result CSVs → TikZ figures, SVG previews, LaTeX tables |
| `compare_thesis.jl` | compare the artifacts with the thesis sources (read-only) |
| `promote_data.jl` | copy a complete full-scale run into `data/` |
| `data/` | stored result CSVs (tracked; provenance in `data/README.md`) |
| `out/` | fresh runs, figures and tables (not tracked) |
| `lib/` | instances and seeds, forcing budgets, aggregation, TikZ and table writers |

`quick` (default) runs reduced digraphs in about a minute. `full` is the thesis scale: the interval experiment is the shortest, the continuous-budget one takes much longer, and the discrete-budget one longer still; runs resume and can be sliced over processes.
