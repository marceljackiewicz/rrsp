# Shared helpers for the experiment drivers: directories, command-line
# options, and small CSV readers/writers (no external dependencies).

using Printf
using Random
using Statistics
using Rrsp

const EXPERIMENTS_DIR = normpath(joinpath(@__DIR__, ".."))

"""Directory with the raw result CSVs that the artifacts are built from.
Default: the tracked `experiments/data`. Override with `RRSP_DATA_DIR`."""
data_root() = get(ENV, "RRSP_DATA_DIR", joinpath(EXPERIMENTS_DIR, "data"))

"""Where a fresh run writes raw CSVs (untracked): `experiments/out/data`.
Override with `RRSP_RUN_DIR`."""
run_root() = get(ENV, "RRSP_RUN_DIR", joinpath(EXPERIMENTS_DIR, "out", "data"))

"""Directory for the generated TeX/SVG artifacts. Default: `experiments/out/artifacts` (untracked).
Override with `RRSP_ARTIFACT_DIR`."""
artifact_root() = get(ENV, "RRSP_ARTIFACT_DIR", joinpath(EXPERIMENTS_DIR, "out", "artifacts"))

# --- command line: `--key=value` and bare `--flag` --------------------------

function parse_options(args::AbstractVector{<:AbstractString} = ARGS)::Dict{String,String}
    opts = Dict{String,String}()
    for a in args
        startswith(a, "--") || error("unrecognized argument $a (expected --key=value)")
        kv = split(a[3:end], '='; limit = 2)
        opts[kv[1]] = length(kv) == 2 ? kv[2] : "true"
    end
    return opts
end

opt(opts, key, default) = get(opts, key, default)
opt_int(opts, key, default) = parse(Int, get(opts, key, string(default)))

"""`quick` (minutes, few draws) or `full` (thesis scale: hours)."""
function opt_scale(opts)::String
    scale = get(opts, "scale", get(ENV, "RRSP_SCALE", "quick"))
    scale in ("quick", "full") || error("--scale must be quick or full, got $scale")
    return scale
end

# --- CSV --------------------------------------------------------------------

function _parse_cell(s::AbstractString)
    v = tryparse(Int, s)
    v !== nothing && return v
    v = tryparse(Float64, s)
    v !== nothing && return v
    return String(s)
end

"""Read a CSV with a header line into a vector of `NamedTuple`s."""
function read_csv(path::AbstractString)::Vector{NamedTuple}
    rows = NamedTuple[]
    open(path) do io
        header = Tuple(Symbol.(split(strip(readline(io)), ',')))
        for line in eachline(io)
            isempty(strip(line)) && continue
            cells = split(line, ',')
            length(cells) == length(header) || error("$path: malformed line: $line")
            push!(rows, NamedTuple{header}(Tuple(_parse_cell(c) for c in cells)))
        end
    end
    return rows
end

"""Read every `*.csv` below `dir` (sorted by path); a missing directory gives no rows."""
function read_csv_tree(dir::AbstractString)::Vector{NamedTuple}
    rows = NamedTuple[]
    isdir(dir) || return rows
    files = String[]
    for (root, _, names) in walkdir(dir)
        for n in names
            endswith(n, ".csv") && push!(files, joinpath(root, n))
        end
    end
    for f in sort(files)
        append!(rows, read_csv(f))
    end
    return rows
end

"""Write `rows` (tuples, in `header` order) atomically."""
function write_csv(path::AbstractString, header::AbstractVector{<:AbstractString}, rows)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    open(tmp, "w") do io
        println(io, join(header, ','))
        for r in rows
            println(io, join(r, ','))
        end
    end
    mv(tmp, path; force = true)
    return path
end

function log_line(args...)
    println(args...)
    flush(stdout)
end
