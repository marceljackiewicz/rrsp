# Tracked result data

Raw results behind the figures and tables of the experiments chapter. From
these files `experiments/make_artifacts.jl` regenerates, in seconds, every
figure and table (`experiments/out/artifacts/`).

| File | Experiment | Rows | Produced by |
|---|---|---|---|
| `interval/replicates.csv` | interval uncertainty, independent costs | 2103 | `interval.jl --scale=full` |
| `interval/bottleneck.csv` | interval uncertainty, designated-path costs | 70 | `interval.jl --scale=full` |
| `cont/all.csv` | continuous budget | 13680 | earlier drivers, see below |
| `disc/all.csv` | discrete budget | 34422 | earlier drivers, see below |

## Columns

* All files: `family` (`layered`, `random_dag`, `grid`, `asp`), `graph` (index of
  the random DAG, `0` for the fixed families), `cost` (index of the cost draw),
  the seeds, `n`, `m` (digraph size), `ell` (fewest arcs of an s–t path),
  `longest` (most arcs), `k` (recovery budget), `z` (optimal value).
* `interval/bottleneck.csv` has no `cost` and `cost_seed`: one instance per family.
* `cont`, `disc`: `npaths` (number of s–t paths), `z_nom = min ĉ(p)`,
  `z_int = min (ĉ+d)(p)`. `cont` adds the forcing budget `gamma_force`, the
  saturation budget `gamma_sat` and the budget `gamma = gamma_frac * gamma_force`;
  `disc` adds `delta_force` and the budget `delta`.
* Once the value reaches the full-recovery value, larger `k` repeat it.

## Provenance

* **`interval`** was recomputed with the current code. It agrees exactly with the
  data of the earlier pipeline (all 2103 replicate rows and the 70 bottleneck
  rows).
* **`cont` and `disc`** were computed by the earlier drivers (a hand-written
  HiGHS C-API implementation of the same enumeration). The columns
  `commit_chat` and `time_sec` were dropped, and the per-instance files were
  concatenated and sorted. The current drivers (`cont.jl`, `disc.jl`, on
  `solve_rrsp_enum`) reproduce these numbers: all 78 rows of one `cont`
  instance (max difference 1e-12) and 30 sampled `disc` points on
  six instances (max difference 1.4e-14). A full `disc` run with the current
  drivers is much slower than the earlier code and has not been repeated.

To replace a file with a fresh full run use `experiments/promote_data.jl`.
