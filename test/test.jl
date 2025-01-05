#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

include("../source/main.jl")

import Test

function testShortestPathSinglePath()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    p::rrsp.Path = rrsp.getShortestPath(instance.graph, instance.s, instance.t)
    Test.@test p.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testShortestPathSingleArcPaths()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    p::rrsp.Path = rrsp.getShortestPath(instance.graph, instance.s, instance.t)
    min_arc::rrsp.Arc = argmin(arc -> arc.cost.first, instance.graph.arcs)
    Test.@test p.arcs == [arc == min_arc for arc in instance.graph.arcs]
end

function testRrspContBudgetSinglePathUsingDagModel()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    solution::rrsp.RrspSolution = rrsp.getRrspContBudgetDag(instance)
    Test.@test solution.first_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
    Test.@test solution.second_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testAspTreeDecomposition()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    tree::rrsp.AspTree = rrsp.getAspTreeDecomposition(instance.graph)

    println("-----------", tree.root_idx)
    for n in tree.nodes
        println("\t", n)
    end

    instance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    tree = rrsp.getAspTreeDecomposition(instance.graph)
    println("-----------", tree.root_idx)
    for n in tree.nodes
        println("\t", n)
    end

    instance = rrsp.parseInstanceFromFile("data/two_beads.rrsp")
    tree = rrsp.getAspTreeDecomposition(instance.graph)
    println("-----------", tree.root_idx)
    for n in tree.nodes
        println("\t", n)
    end
end

# testShortestPathSinglePath()
# testShortestPathSingleArcPaths()
# testRrspContBudgetSinglePathUsingDagModel()
testAspTreeDecomposition()
