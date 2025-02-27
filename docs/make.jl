push!(LOAD_PATH, "../src/")

using Documenter, Rrsp

makedocs(sitename="Robust Recoverable Shortest Path Solver", remotes=nothing)
