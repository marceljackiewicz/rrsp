# LaTeX tables in the thesis layout. The head and tail text of each table comes
# from lib/tex_templates.jl; only the numbers are generated.

"""Table entry: two decimals, with `zero` for a value that rounds to zero."""
function entry(x::Real; zero::AbstractString = "0")::String
    s = @sprintf("%.2f", x)
    return "\$" * (s == "0.00" || s == "-0.00" ? zero : s) * "\$"
end

grid_row(label::AbstractString, c::Curve; zero = "0") =
    label * " & " * join((entry(interpolate(c, t); zero = zero) for t in GRID), " & ") * " \\\\"

"""Table of the mean VoR on the grid `j/12` (interval or bottleneck experiment)."""
function vor_table(head::AbstractString, tail::AbstractString, curves::AbstractDict{String,Curve})::String
    io = IOBuffer()
    print(io, head)
    for f in FAMILY_ORDER
        haskey(curves, f) && println(io, grid_row(FAMILY_ROW[f], curves[f]))
    end
    print(io, tail)
    return String(take!(io))
end

"""Table with a block of four family rows per budget fraction."""
function budget_vor_table(head::AbstractString, tail::AbstractString, curves::AbstractDict; zero::AbstractString)::String
    io = IOBuffer()
    print(io, head)
    for (b, frac) in enumerate(BUDGET_FRACTIONS)
        b > 1 && println(io, "\\hline")
        for (i, f) in enumerate(FAMILY_ORDER)
            lead = i == 1 ? "\$$(frac)\$" : ""
            println(io, grid_row("$lead & $(FAMILY_ROW[f])", curves[frac][f]; zero = zero))
        end
    end
    print(io, tail)
    return String(take!(io))
end

const INSTANCE_ROW_ORDER = (("layered", "Layered"), ("grid", "Grid"), ("asp", "Arc-series-parallel"), ("random_dag", "Random DAG"))

range_cell(r::SizeRange) = r.lo == r.hi ? "\$$(r.lo)\$" : "\$$(r.lo)\$--\$$(r.hi)\$"

"""Table of digraph sizes; `fields` lists the columns after the family."""
function instance_table(head::AbstractString, tail::AbstractString, sizes, fields::Vector{Symbol})::String
    io = IOBuffer()
    print(io, head)
    for (f, label) in INSTANCE_ROW_ORDER
        haskey(sizes, f) || continue
        println(io, "        ", label, " & ", join((range_cell(sizes[f][c]) for c in fields), " & "), " \\\\")
    end
    print(io, tail)
    return String(take!(io))
end
