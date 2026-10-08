# Compare the generated figures and tables with the ones in the thesis sources.
#
#   julia --project=experiments experiments/compare_thesis.jl --thesis=DIR [--artifacts=DIR]
#
# DIR is the thesis directory with experiments.tex and figures/experiments/.
# Nothing in DIR is modified. For every file the report says "identical" or
# lists the differing lines; the exit code is 1 if any file differs.

include(joinpath(@__DIR__, "lib", "common.jl"))

const FIGURE_FILES = ["fig_scale_vor", "fig_scale_vor_bottleneck", "fig_scale_vor_cont", "fig_scale_vor_disc"]
const TABLE_LABELS = Dict(
    "interval_vor" => "tab:experiments-interval-vor",
    "bottleneck_vor" => "tab:experiments-bottleneck-vor",
    "cont_vor" => "tab:experiments-cont-vor",
    "disc_vor" => "tab:experiments-disc-vor",
    "interval_instances" => "tab:experiments-interval-instances",
    "cont_instances" => "tab:experiments-cont-instances",
    "disc_instances" => "tab:experiments-disc-instances",
)

"""The `table` environment of `src` that contains `label`."""
function table_block(src::AbstractString, label::AbstractString)::Union{Nothing,String}
    for m in eachmatch(r"\\begin\{table\}.*?\\end\{table\}"s, src)
        occursin("\\label{" * label * "}", m.match) && return String(m.match)
    end
    return nothing
end

function compare_lines(generated::AbstractString, thesis::AbstractString)
    a = split(strip(generated, '\n'), '\n')
    b = split(strip(thesis, '\n'), '\n')
    diffs = String[]
    for i in 1:max(length(a), length(b))
        x = i <= length(a) ? a[i] : "<missing>"
        y = i <= length(b) ? b[i] : "<missing>"
        x == y || push!(diffs, "  line $i\n    generated: $x\n    thesis:    $y")
    end
    return diffs
end

function main(opts)
    thesis = get(opts, "thesis", "")
    isempty(thesis) && error("give the thesis directory with --thesis=DIR")
    art = opt(opts, "artifacts", artifact_root())
    src = read(joinpath(thesis, "experiments.tex"), String)
    failures = 0
    report(name, diffs) = begin
        if isempty(diffs)
            log_line("identical   ", name)
        else
            failures += 1
            log_line("DIFFERENT   ", name, " ($(length(diffs)) lines)")
            foreach(d -> log_line(d), first(diffs, 5))
        end
    end
    for f in FIGURE_FILES
        report("figures/$f.tex", compare_lines(read(joinpath(art, "figures", f * ".tex"), String),
                                              read(joinpath(thesis, "figures", "experiments", f * ".tex"), String)))
    end
    for (name, label) in sort(collect(TABLE_LABELS))
        block = table_block(src, label)
        block === nothing && (failures += 1; log_line("MISSING     $label in experiments.tex"); continue)
        report("tables/tab_$name.tex", compare_lines(read(joinpath(art, "tables", "tab_$name.tex"), String), block))
    end
    return failures
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(parse_options()) == 0 ? 0 : 1)
end
