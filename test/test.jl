#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

include("../source/main.jl")

import Test

function testShortestPathSinglePath()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    p::rrsp.Path = rrsp.solveShortestPath(instance.graph, instance.s_idx, instance.t_idx)
    Test.@test p.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testShortestPathSingleArcPaths()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    p::rrsp.Path = rrsp.solveShortestPath(instance.graph, instance.s_idx, instance.t_idx)
    min_arc::rrsp.Arc = argmin(arc -> arc.cost.first, instance.graph.arcs)
    Test.@test p.arcs == [arc == min_arc for arc in instance.graph.arcs]
end

function testIncrementalShortestPathSingleArcPaths()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    first_stage_path::rrsp.Path = rrsp.solveShortestPath(instance.graph, instance.s_idx, instance.t_idx)  # any first stage path as x
    inc_path::rrsp.Path = rrsp.solveIncrementalShortestPath(instance, first_stage_path)
    Test.@test first_stage_path.arcs == inc_path.arcs  # just one path to choose from
end

function testRecoverableShortestPathSingleArcPaths()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    solution::rrsp.RrspSolution = rrsp.solveRecoverableShortestPath(instance)
    Test.@test solution.first_stage_path.arcs == solution.second_stage_path.arcs  # just one path to choose from
end


function testRecoverableShortestPathFixedThetaCounterexample()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/fixed_theta_counterexample.rrsp")
    solution::rrsp.RrspSolution = rrsp.solveRecoverableShortestPath(instance)
    Test.@test solution.second_stage_path.arcs == [1, 1, 0, 0, 0]
end

function testRrspDagContBudgetSinglePathUsingModel()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    solution::rrsp.RrspSolution = rrsp.solveRrspContBudgetDag(instance)
    Test.@test solution.first_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
    Test.@test solution.second_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testRrspGeneralContBudgetSinglePathUsingModel()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_path.rrsp")
    solution::rrsp.RrspSolution = rrsp.solveRrspContBudget(instance)
    Test.@test solution.first_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
    Test.@test solution.second_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testRrspApproxRatioAcyclic()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    ratio_discrete::Float64 = rrsp.getRecSpToRrspAcyclicDiscreteBudgetApproxRatio(instance)
    ratio_cont::Float64 = rrsp.getRecSpToRrspAcyclicContBudgetApproxRatio(instance)
end

function areTwoSolutionsEqual(s1::rrsp.RrspSolution, s2::rrsp.RrspSolution)
    return (s1.first_stage_path.arcs == s2.first_stage_path.arcs) && (s1.second_stage_path.arcs == s2.second_stage_path.arcs) && (s1.value == s2.value)
end

function testRecSpAsp()
    instance::rrsp.RrspInstance = rrsp.parseInstanceFromFile("data/fixed_theta_counterexample.rrsp")
    solution_combinatorial::rrsp.RrspSolution = rrsp.solveRecSpInAsp(instance)
    solution_model::rrsp.RrspSolution = rrsp.solveRecoverableShortestPath(instance)

    Test.@test areTwoSolutionsEqual(solution_combinatorial, solution_model)
end

testShortestPathSinglePath()
testShortestPathSingleArcPaths()
testIncrementalShortestPathSingleArcPaths()
testRecoverableShortestPathSingleArcPaths()
testRecoverableShortestPathFixedThetaCounterexample()
testRrspDagContBudgetSinglePathUsingModel()
testRrspGeneralContBudgetSinglePathUsingModel()
testRrspApproxRatioAcyclic()
testRecSpAsp()
