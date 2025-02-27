#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

include("../src/Rrsp.jl")

import Test

function testShortestPathSinglePath()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_path.rrsp")
    p::Rrsp.Path = Rrsp.solveShortestPath(instance.graph, instance.s_idx, instance.t_idx)
    Test.@test p.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testShortestPathSingleArcPaths()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    p::Rrsp.Path = Rrsp.solveShortestPath(instance.graph, instance.s_idx, instance.t_idx)
    min_arc::Rrsp.Arc = argmin(arc -> arc.cost.first, instance.graph.arcs)
    Test.@test p.arcs == [arc == min_arc for arc in instance.graph.arcs]
end

function testIncrementalShortestPathSingleArcPaths()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    first_stage_path::Rrsp.Path = Rrsp.solveShortestPath(instance.graph, instance.s_idx, instance.t_idx)  # any first stage path as x
    inc_path::Rrsp.Path = Rrsp.solveIncrementalShortestPath(instance, first_stage_path)
    Test.@test first_stage_path.arcs == inc_path.arcs  # just one path to choose from
end

function testRecoverableShortestPathSingleArcPaths()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    solution::Rrsp.RrspSolution = Rrsp.solveRecoverableShortestPath(instance)
    Test.@test solution.first_stage_path.arcs == solution.second_stage_path.arcs  # just one path to choose from
end


function testRecoverableShortestPathFixedThetaCounterexample()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/fixed_theta_counterexample.rrsp")
    solution::Rrsp.RrspSolution = Rrsp.solveRecoverableShortestPath(instance)
    Test.@test solution.second_stage_path.arcs == [1, 1, 0, 0, 0]
end

function testRrspDagContBudgetSinglePathUsingModel()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_path.rrsp")
    solution::Rrsp.RrspSolution = Rrsp.solveRrspContBudgetDag(instance)
    Test.@test solution.first_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
    Test.@test solution.second_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testRrspGeneralContBudgetSinglePathUsingModel()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_path.rrsp")
    solution::Rrsp.RrspSolution = Rrsp.solveRrspContBudget(instance)
    Test.@test solution.first_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
    Test.@test solution.second_stage_path.arcs == [1 for _ in 1:length(instance.graph.arcs)]
end

function testRrspApproxRatioAcyclic()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/single_arc_paths.rrsp")
    ratio_discrete::Float64 = Rrsp.getRecSpToRrspAcyclicDiscreteBudgetApproxRatio(instance)
    ratio_cont::Float64 = Rrsp.getRecSpToRrspAcyclicContBudgetApproxRatio(instance)
end

function areTwoSolutionsEqual(s1::Rrsp.RrspSolution, s2::Rrsp.RrspSolution)
    return (s1.first_stage_path.arcs == s2.first_stage_path.arcs) && (s1.second_stage_path.arcs == s2.second_stage_path.arcs) && (s1.value == s2.value)
end

function testRecSpAsp()
    instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/fixed_theta_counterexample.rrsp")
    solution_combinatorial::Rrsp.RrspSolution = Rrsp.solveRecSpInAsp(instance)
    solution_model::Rrsp.RrspSolution = Rrsp.solveRecoverableShortestPath(instance)

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
