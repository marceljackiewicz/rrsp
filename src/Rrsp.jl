#= Robust Recoverable Shortest Path Solver
#
#  Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński
=#

module Rrsp

import JuMP

optimizer::Any = nothing  # should be set before using the solver

include("types.jl")
export Path
export RrspInstance
export RrspSolution

include("parsing.jl")
export parseInstanceFromFile

include("deterministic.jl")
export solveDeterministicShortestPath

include("incremental.jl")
export solveIncrementalShortestPath

include("recoverable.jl")
export solveRecoverableShortestPath
export solveRecoverableShortestPathInAsp

include("recoverable_robust.jl")
export solveRrspContBudget
export solveRrspContBudgetDag

end  # module Rrsp
