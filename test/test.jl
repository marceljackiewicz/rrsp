#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

include("../source/main.jl")

import Test

function testShortestPathSinglePath()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    p::rrsp.Path = rrsp.getShortestPath(instance.graph, instance.s, instance.t)
    Test.@test p.arcs == instance.graph.arcs
end

function testShortestPathSingleArcPaths()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    p::rrsp.Path = rrsp.getShortestPath(instance.graph, instance.s, instance.t)
    Test.@test p.arcs == [argmin(arc -> arc.cost.first, instance.graph.arcs)]
end

function testRrspContBudgetSinglePathUsingDagModel()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    solution::rrsp.RrspSolution = rrsp.getRrspContBudgetDag(instance)
    Test.@test solution.first_stage_path == solution.first_stage_path  # fix the test
end

testShortestPathSinglePath()
testShortestPathSingleArcPaths()
testRrspContBudgetSinglePathUsingDagModel()
