#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

module graph_db

include("../source/main.jl")

function getEmptyGraph()::rrsp.Graph
    return rrsp.Graph(rrsp.Node[], rrsp.Arc[])
end

function getSingleArcGraph(c::rrsp.Cost)::rrsp.Graph
    n1::rrsp.Node = rrsp.Node(0);
    n2::rrsp.Node = rrsp.Node(1);
    nodes::Vector{rrsp.Node} = [n1, n2]

    arc::rrsp.Arc = rrsp.Arc(n1, n2, c);
    arcs::Vector{rrsp.Arc} = [arc]

    return rrsp.Graph(nodes, arcs)
end

end  # module graph_db
