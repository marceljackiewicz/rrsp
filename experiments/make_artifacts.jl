# Figures and tables of the experiments chapter from the result CSVs.
#
#   julia --project=experiments experiments/make_artifacts.jl
#         [--data=DIR] [--out=DIR] [--no-svg] [--pdf]
#
# --data   directory with interval/, cont/, disc/ (default: the tracked
#          experiments/data; use experiments/out/data for a fresh run)
# --out    output directory (default: experiments/out/artifacts)
# --pdf    also compile preview.pdf (needs pdflatex; the thesis macros are
#          replaced by plain stand-ins)
#
# Writes below --out:
#   figures/fig_scale_vor{,_bottleneck,_cont,_disc}.tex   (+ .svg previews)
#   tables/tab_{interval,bottleneck,cont,disc}_vor.tex
#   tables/tab_{interval,cont,disc}_instances.tex

@isdefined(EXPERIMENTS_DIR) || include(joinpath(@__DIR__, "lib", "common.jl"))
@isdefined(GRID_STEPS) || include(joinpath(@__DIR__, "lib", "aggregate.jl"))
@isdefined(single_panel_figure) || include(joinpath(@__DIR__, "lib", "tikz.jl"))
@isdefined(vor_table) || include(joinpath(@__DIR__, "lib", "tables.jl"))

function write_text(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    write(path, text)
    log_line("wrote $path")
    return path
end

function require_csv(path::AbstractString)
    isfile(path) || error("missing $path (run the experiment or point --data to the result directory)")
    return read_csv(path)
end

function require_tree(dir::AbstractString)
    rows = read_csv_tree(dir)
    isempty(rows) && error("no CSV files below $dir (run the experiment or point --data to the result directory)")
    return rows
end

const PREVIEW_PREAMBLE = raw"""
\documentclass{article}
\usepackage[margin=1.6cm]{geometry}
\usepackage{amsmath,amssymb,xcolor,tikz}
\providecommand{\NeighborhoodSize}{k}
\providecommand{\STPath}{s\text{--}t}
\providecommand{\FirstStageCostVector}{C}
\providecommand{\NominalSecondStageCostVector}{\hat c}
\providecommand{\VerticesSet}{V}
\providecommand{\ArcsSet}{A}
\providecommand{\RobProblemShort}{\mathrm{Rob}}
\providecommand{\IntervalUncertaintySet}{\mathcal{U}}
\begin{document}
"""

"""Compile `preview.pdf` from the generated files; returns its path or `nothing`."""
function make_preview(out::AbstractString, files::Vector{String})
    Sys.which("pdflatex") === nothing && (@warn "pdflatex not found; skipping --pdf"; return nothing)
    tex = joinpath(out, "preview.tex")
    open(tex, "w") do io
        print(io, PREVIEW_PREAMBLE)
        for f in files
            endswith(f, ".tex") || continue
            println(io, "\\input{", relpath(f, out), "}")
            println(io, "\\clearpage")
        end
        println(io, "\\end{document}")
    end
    for _ in 1:2
        run(pipeline(Cmd(`pdflatex -interaction=nonstopmode -halt-on-error preview.tex`; dir = out); stdout = devnull))
    end
    for ext in ("aux", "log")
        rm(joinpath(out, "preview." * ext); force = true)
    end
    log_line("wrote ", joinpath(out, "preview.pdf"))
    return joinpath(out, "preview.pdf")
end

"""Generate every figure and table. Returns the list of written files."""
function make_artifacts(; data::AbstractString = data_root(), out::AbstractString = artifact_root(), svg::Bool = true)
    files = String[]
    emit(rel, text) = push!(files, write_text(joinpath(out, rel), text))

    # Interval uncertainty, independent costs
    replicates = require_csv(joinpath(data, "interval", "replicates.csv"))
    interval = interval_instance_curves(replicates)
    interval_avg = Dict(f => family_curve(filter(ic -> ic.family == f, interval)) for f in FAMILY_ORDER)
    interval_grid = Dict(f => family_curve(filter(ic -> ic.family == f, interval); grid = true) for f in FAMILY_ORDER)
    ymax, ystep = nice_axis(max_y(interval_avg))
    emit("figures/fig_scale_vor.tex", single_panel_figure(interval_avg, FIGURE_INTERVAL_TAIL; ymax, ystep, digits = 4))
    svg && emit("figures/fig_scale_vor.svg", single_panel_svg(interval_avg; ymax, ystep))
    emit("tables/tab_interval_vor.tex", vor_table(TABLE_INTERVAL_VOR_HEAD, TABLE_INTERVAL_VOR_TAIL, interval_grid))
    emit("tables/tab_interval_instances.tex",
         instance_table(TABLE_INTERVAL_INSTANCES_HEAD, TABLE_INTERVAL_INSTANCES_TAIL, instance_sizes(replicates), [:n, :m, :ell, :longest]))

    # Interval uncertainty, designated-path (bottleneck) costs
    bottleneck = bottleneck_curves(require_csv(joinpath(data, "interval", "bottleneck.csv")))
    ymax, ystep = nice_axis(max_y(bottleneck))
    emit("figures/fig_scale_vor_bottleneck.tex", single_panel_figure(bottleneck, FIGURE_BOTTLENECK_TAIL; ymax, ystep, digits = 3))
    svg && emit("figures/fig_scale_vor_bottleneck.svg", single_panel_svg(bottleneck; ymax, ystep))
    emit("tables/tab_bottleneck_vor.tex", vor_table(TABLE_BOTTLENECK_VOR_HEAD, TABLE_BOTTLENECK_VOR_TAIL, bottleneck))

    # Continuous budget
    cont_rows = require_tree(joinpath(data, "cont"))
    cont = budget_family_curves(continuous_instance_curves, cont_rows)
    ymax, ystep = nice_axis(max_y(cont))
    emit("figures/fig_scale_vor_cont.tex", budget_figure(cont, FIGURE_CONT_TAIL; ymax, ystep))
    svg && emit("figures/fig_scale_vor_cont.svg", budget_svg(cont; ymax, ystep))
    emit("tables/tab_cont_vor.tex", budget_vor_table(TABLE_CONT_VOR_HEAD, TABLE_CONT_VOR_TAIL, cont; zero = "0"))
    emit("tables/tab_cont_instances.tex",
         instance_table(TABLE_CONT_INSTANCES_HEAD, TABLE_CONT_INSTANCES_TAIL, instance_sizes(cont_rows), [:n, :m, :npaths, :ell, :longest]))

    # Discrete budget
    disc_rows = require_tree(joinpath(data, "disc"))
    disc = budget_family_curves(discrete_instance_curves, disc_rows)
    ymax, ystep = nice_axis(max_y(disc))
    emit("figures/fig_scale_vor_disc.tex", budget_figure(disc, FIGURE_DISC_TAIL; ymax, ystep))
    svg && emit("figures/fig_scale_vor_disc.svg", budget_svg(disc; ymax, ystep))
    emit("tables/tab_disc_vor.tex", budget_vor_table(TABLE_DISC_VOR_HEAD, TABLE_DISC_VOR_TAIL, disc; zero = "0.00"))
    emit("tables/tab_disc_instances.tex",
         instance_table(TABLE_DISC_INSTANCES_HEAD, TABLE_DISC_INSTANCES_TAIL, instance_sizes(disc_rows), [:n, :m, :npaths, :ell, :longest]))
    return files
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_options()
    out = opt(opts, "out", artifact_root())
    files = make_artifacts(; data = opt(opts, "data", data_root()), out = out, svg = !haskey(opts, "no-svg"))
    haskey(opts, "pdf") && make_preview(out, files)
end
