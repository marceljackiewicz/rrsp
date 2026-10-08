# TikZ figures in the style of the thesis (macros such as \NeighborhoodSize come
# from the thesis preamble) and plain SVG previews of the same curves.

using Printf

include(joinpath(@__DIR__, "tex_templates.jl"))

function f4(x)
    s = @sprintf("%.4f", x)
    return s == "-0.0000" ? "0.0000" : s      # rounding noise around zero
end
f2(x) = @sprintf("%0.2f", x)
f3(x) = @sprintf("%0.3f", x)
fg(x) = @sprintf("%g", x)

"""Break points as `(x,y)` pairs, 4 decimals (`digits = 4`) or shortest 3-decimal form."""
function coordinates(c::Curve; digits::Int = 4)::String
    if digits == 4
        return join(("($(f4(c.x[i])),$(f4(c.y[i])))" for i in eachindex(c.x)), " ")
    end
    return join(("($(round(c.x[i]; digits = digits)),$(round(c.y[i]; digits = digits)))" for i in eachindex(c.x)), " ")
end

# --- y axis ranges -------------------------------------------------------------------

"""Largest y value of a `Curve`, or of all curves in a (possibly nested) dictionary."""
max_y(c::Curve) = maximum(c.y)
max_y(d::AbstractDict) = maximum(max_y, values(d))

"""
    nice_axis(ymax_data) -> (ymax, ystep)

Smallest integer axis range `0:ystep:ymax` with at most five ticks above zero
that contains `ymax_data`; `ystep` is 1, 2 or 5 times a power of ten.
"""
function nice_axis(ymax_data::Real)
    for e in 0:6, m in (1, 2, 5)
        step = m * 10^e
        n = ceil(Int, ymax_data / step)
        n <= 5 && return (max(n, 1) * step, step)
    end
    error("y range too large: $ymax_data")
end

# --- single-panel figures: interval and bottleneck ----------------------------------

const FAMILY_STYLE = Dict(
    "layered" => ("cLayered", "line width=1.4pt"),
    "random_dag" => ("cRandomDag", "line width=1.4pt, dashed"),
    "grid" => ("cGrid", "line width=2.1pt, line cap=round, dash pattern=on 0pt off 5.5pt"),
    "asp" => ("cAsp", "line width=1.6pt, dash pattern=on 7pt off 3pt on 1pt off 3pt"),
)
const FAMILY_LEGEND = Dict("layered" => "layered", "random_dag" => "random DAG", "grid" => "grid", "asp" => "ASP, \$\\NeighborhoodSize/(3\\ell)\$")

"""
    single_panel_figure(curves, tail; ymax, ystep, digits) -> String

One axes with a curve per family (`curves::Dict{String,Curve}`), a legend to
the right, and the caption/label text `tail`.
"""
function single_panel_figure(curves::AbstractDict{String,Curve}, tail::AbstractString; ymax::Int, ystep::Int, digits::Int)::String
    ycm = 4.4
    io = IOBuffer()
    println(io, "\\begin{figure}[h!]")
    println(io, "\\begin{center}")
    println(io, "\\definecolor{cLayered}{RGB}{0,0,0}")
    println(io, "\\definecolor{cRandomDag}{RGB}{0,114,178}")
    println(io, "\\definecolor{cGrid}{RGB}{0,130,60}")
    println(io, "\\definecolor{cAsp}{RGB}{180,20,30}")
    println(io, "\\begin{tikzpicture}[x=8.0000cm, y=$(f4(ycm / ymax))cm]")
    println(io, "\\draw (0,0) rectangle (1,$ymax);")
    println(io, "\\foreach \\y in {0,$ystep,...,$ymax} {")
    println(io, "    \\draw ([xshift=-1.2mm]0,\\y) -- (0,\\y);")
    println(io, "    \\node[left] at ([xshift=-1.6mm]0,\\y) {\\y};")
    println(io, "}")
    println(io, "\\foreach \\x in {0,0.25,0.5,0.75,1} {")
    println(io, "    \\draw (\\x,0) -- ([yshift=-1.2mm]\\x,0);")
    println(io, "    \\node[below] at ([yshift=-1.5mm]\\x,0) {\\x};")
    println(io, "}")
    println(io, "\\node[rotate=90, anchor=south] at ([xshift=-11mm]0,$(ymax / 2)) {value of recovery (\\%)};")
    println(io, "\\node[below] at ([yshift=-8mm]0.5,0) {\$\\NeighborhoodSize/\\ell\$};")
    families = [f for f in FAMILY_ORDER if haskey(curves, f)]
    for f in families
        col, style = FAMILY_STYLE[f]
        println(io, "\\draw[$col, $style] plot coordinates {$(coordinates(curves[f]; digits = digits))};")
    end
    println(io, "\\begin{scope}[x=1cm, y=1cm, xshift=8.45cm, yshift=$(f3(ycm - 0.15))cm]")
    for (i, f) in enumerate(families)
        col, style = FAMILY_STYLE[f]
        yl = f2(-0.48 * (i - 1))
        println(io, "\\draw[$col, $style] (0,$yl) -- (1.25,$yl);")
        println(io, "\\node[anchor=west] at (1.40,$yl) {$(FAMILY_LEGEND[f])};")
    end
    println(io, "\\end{scope}")
    println(io, "\\end{tikzpicture}")
    print(io, tail)
    return String(take!(io))
end

# --- four-panel figures: continuous and discrete budgets -------------------------------

const FRACTION_STYLE = Dict(
    0.2 => ("cF2", "line width=1.2pt, dashed"),
    0.4 => ("cF4", "line width=1.4pt"),
    0.6 => ("cF6", "line width=1.3pt, dash pattern=on 6pt off 2.5pt"),
    0.8 => ("cF8", "line width=1.2pt, dash pattern=on 1.2pt off 2pt"),
)
const PANEL_POSITION = Dict("layered" => (0.0, 0.0), "random_dag" => (5.15, 0.0), "grid" => (0.0, -5.85), "asp" => (5.15, -5.85))

"""
    budget_figure(curves, tail; ymax, ystep) -> String

Four panels (one per family) with a curve per budget fraction;
`curves[frac][family]` is on the grid `j/12`.
"""
function budget_figure(curves::AbstractDict, tail::AbstractString; ymax::Int, ystep::Int)::String
    io = IOBuffer()
    println(io, "\\begin{figure}[h!]")
    println(io, "\\begin{center}")
    println(io, "\\definecolor{cF2}{RGB}{0,114,178}")
    println(io, "\\definecolor{cF4}{RGB}{0,0,0}")
    println(io, "\\definecolor{cF6}{RGB}{0,130,60}")
    println(io, "\\definecolor{cF8}{RGB}{180,20,30}")
    println(io, "\\begin{tikzpicture}[x=1cm, y=1cm]")
    yscale = 3.1 / ymax
    yticks = join(0:ystep:ymax, ",")
    for f in FAMILY_ORDER
        px, py = PANEL_POSITION[f]
        left = px == 0.0
        println(io, "\\begin{scope}[xshift=$(f2(px))cm, yshift=$(f2(py))cm, x=4.00cm, y=$(f3(yscale))cm]")
        println(io, "\\draw (0,0) rectangle (1,$ymax);")
        println(io, "\\foreach \\y in {$yticks} {")
        println(io, "    \\draw ([xshift=-1.1mm]0,\\y) -- (0,\\y);")
        println(io, "    \\node[left] at ([xshift=-1.4mm]0,\\y) {\\y};")
        println(io, "}")
        println(io, "\\foreach \\x in {0,0.5,1} {")
        println(io, "    \\draw (\\x,0) -- ([yshift=-1.1mm]\\x,0);")
        println(io, "    \\node[below] at ([yshift=-1.3mm]\\x,0) {\\x};")
        println(io, "}")
        println(io, "\\node[anchor=south] at (0.5,$(fg(ymax + 0.6))) {$(FAMILY_TITLE[f])};")
        left && println(io, "\\node[rotate=90, anchor=south] at ([xshift=-9mm]0,$(fg(ymax / 2))) {value of recovery (\\%)};")
        xlabel = f == "asp" ? "\$\\NeighborhoodSize/(3\\ell)\$" : "\$\\NeighborhoodSize/\\ell\$"
        println(io, "\\node[below] at ([yshift=-6.5mm]0.5,0) {$xlabel};")
        for frac in BUDGET_FRACTIONS
            col, style = FRACTION_STYLE[frac]
            println(io, "\\draw[$col, $style] plot coordinates {$(coordinates(curves[frac][f]))};")
        end
        println(io, "\\end{scope}")
    end
    println(io, "\\begin{scope}[xshift=10.5cm, yshift=-3.6cm]")
    for (i, frac) in enumerate(BUDGET_FRACTIONS)
        col, style = FRACTION_STYLE[frac]
        yl = i == 1 ? "0" : f2(-0.45 * (i - 1))
        println(io, "\\draw[$col, $style] (0,$yl) -- (1.15,$yl);")
        println(io, "\\node[anchor=west] at (1.30,$yl) {\$$(frac)\$};")
    end
    println(io, "\\end{scope}")
    println(io, "\\end{tikzpicture}")
    print(io, tail)
    return String(take!(io))
end

# --- SVG previews --------------------------------------------------------------------------

const SVG_COLOR = Dict("layered" => "#000000", "random_dag" => "#0072b2", "grid" => "#00823c", "asp" => "#b41428")
const FRACTION_COLOR = Dict(0.2 => "#0072b2", 0.4 => "#000000", 0.6 => "#00823c", 0.8 => "#b41428")
const SVG_DASH = Dict("layered" => "", "random_dag" => "6 4", "grid" => "2 4", "asp" => "10 4 2 4",
                      0.2 => "6 4", 0.4 => "", 0.6 => "8 3", 0.8 => "2 3")

"""One SVG axes (origin at `(ox, oy)`, width `w`, height `h`) with its curves."""
function svg_axes(io, ox, oy, w, h, ymax, ystep, title, xlabel, curves)
    println(io, "<rect x=\"$ox\" y=\"$oy\" width=\"$w\" height=\"$h\" fill=\"none\" stroke=\"#444\"/>")
    for y in 0:ystep:ymax
        py = oy + h - h * y / ymax
        println(io, "<line x1=\"$(ox - 4)\" y1=\"$py\" x2=\"$ox\" y2=\"$py\" stroke=\"#444\"/>")
        println(io, "<text x=\"$(ox - 7)\" y=\"$(py + 4)\" text-anchor=\"end\" font-size=\"11\">$y</text>")
    end
    for x in (0, 0.25, 0.5, 0.75, 1)
        px = ox + w * x
        println(io, "<line x1=\"$px\" y1=\"$(oy + h)\" x2=\"$px\" y2=\"$(oy + h + 4)\" stroke=\"#444\"/>")
        println(io, "<text x=\"$px\" y=\"$(oy + h + 16)\" text-anchor=\"middle\" font-size=\"11\">$x</text>")
    end
    println(io, "<text x=\"$(ox + w / 2)\" y=\"$(oy + h + 32)\" text-anchor=\"middle\" font-size=\"12\">$xlabel</text>")
    isempty(title) || println(io, "<text x=\"$(ox + w / 2)\" y=\"$(oy - 6)\" text-anchor=\"middle\" font-size=\"13\">$title</text>")
    for (c, color, dash) in curves
        pts = join(("$(round(ox + w * c.x[i]; digits = 2)),$(round(oy + h - h * min(c.y[i], ymax) / ymax; digits = 2))" for i in eachindex(c.x)), " ")
        da = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
        println(io, "<polyline points=\"$pts\" fill=\"none\" stroke=\"$color\" stroke-width=\"2\"$da/>")
    end
end

function svg_legend(io, x, y, entries)
    for (i, (label, color, dash)) in enumerate(entries)
        py = y + 20 * (i - 1)
        da = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
        println(io, "<line x1=\"$x\" y1=\"$py\" x2=\"$(x + 36)\" y2=\"$py\" stroke=\"$color\" stroke-width=\"2\"$da/>")
        println(io, "<text x=\"$(x + 44)\" y=\"$(py + 4)\" font-size=\"12\">$label</text>")
    end
end

svg_wrap(w, h, body) = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$w\" height=\"$h\" viewBox=\"0 0 $w $h\" font-family=\"sans-serif\">\n" *
                       "<rect width=\"$w\" height=\"$h\" fill=\"white\"/>\n" * body * "</svg>\n"

function single_panel_svg(curves::AbstractDict{String,Curve}; ymax::Int, ystep::Int)::String
    io = IOBuffer()
    fams = [f for f in FAMILY_ORDER if haskey(curves, f)]
    svg_axes(io, 60, 20, 420, 240, ymax, ystep, "", "k / ℓ  (ASP: k / 3ℓ)", [(curves[f], SVG_COLOR[f], SVG_DASH[f]) for f in fams])
    svg_legend(io, 510, 40, [(FAMILY_TITLE[f], SVG_COLOR[f], SVG_DASH[f]) for f in fams])
    println(io, "<text transform=\"translate(16,140) rotate(-90)\" text-anchor=\"middle\" font-size=\"12\">value of recovery (%)</text>")
    return svg_wrap(640, 310, String(take!(io)))
end

function budget_svg(curves::AbstractDict; ymax::Int, ystep::Int)::String
    io = IOBuffer()
    pos = Dict("layered" => (60, 30), "random_dag" => (330, 30), "grid" => (60, 330), "asp" => (330, 330))
    for f in FAMILY_ORDER
        ox, oy = pos[f]
        xlabel = f == "asp" ? "k / 3ℓ" : "k / ℓ"
        svg_axes(io, ox, oy, 220, 220, ymax, ystep, FAMILY_TITLE[f], xlabel,
                 [(curves[frac][f], FRACTION_COLOR[frac], SVG_DASH[frac]) for frac in BUDGET_FRACTIONS])
    end
    svg_legend(io, 600, 60, [("fraction $frac", FRACTION_COLOR[frac], SVG_DASH[frac]) for frac in BUDGET_FRACTIONS])
    println(io, "<text transform=\"translate(14,140) rotate(-90)\" text-anchor=\"middle\" font-size=\"12\">value of recovery (%)</text>")
    println(io, "<text transform=\"translate(14,440) rotate(-90)\" text-anchor=\"middle\" font-size=\"12\">value of recovery (%)</text>")
    return svg_wrap(760, 620, String(take!(io)))
end
