# RRSP Solver Documentation

The purpose of this document is to illustrate the user API of Robust Recoverable Shortest Path Solver.
Technical documentation of the solver is left in the source code in the form of docstrings and inline comments.

The documentation adheres to the `0.1.0` version of `Rrsp` package.

## Contents
```@contents
```

## How to install
Use `Pkg` package manager to install this package for a session.

1. Clone the repository.
2. Add the package to local registry by
    - `julia`
    - `import Pkg; Pkg.develop(path="./Rrsp")`
3. Build the package by importing it for the first time:
    - `import Rrsp`

### Dependencies
The only pacakge dependency of `Rrsp` is [JuMP](https://jump.dev/JuMP.jl/stable/).
Users are responsible for installing MIP solvers they want to use with `Rrsp`.

## API Data Types
For simplicity, the API uses only one object type for input and output, respectively.
```@docs
RrspInstance
RrspSolution
```

The `Path` structure, the building block of the output object type, is left transparent:
```@docs
Path
```
To fetch the arcs of which the solution path consist you can use the `Path.arcs` vector.

## Setting the Solver

For now choosing the solver is not allowed.

## Parsing Input
The instance of RRSP is assumed to be given in a text file.
Parsing the instance is separate from solving the problem.
Hence, you can keep multiple `RrspInstance` objects at runtime at once.

```@docs
parseInstanceFromFile
```

## Solving Shortest Path Problems
Once the instance is deserialized you can choose which problem should be solved with it.
Solving the instance doesn't change its data so one instance object can be resused
to obtain solutions to different variants of the problem.

Functions solving RRSP and the subproblems are grouped here according to the problem type:

### Deterministic SP
```@docs
solveDeterministicShortestPath
```

### Incremental SP
```@docs
solveIncrementalShortestPath
```

### Recoverable SP
```@docs
solveRecoverableShortestPath
solveRecoverableShortestPathInAsp
```

### Recoverable Robust SP
```@docs
solveRrspContBudget
solveRrspContBudgetDag
```

## Working Example
To illustrate the use of API we present an example of solving a small problem instance.

Assume that there exists a file under `data/instance.rrsp` path containing a problem instance written down as shown in Parsing Input section.\
Make sure `Rrsp` definitions are visible. Use `import` to avoid name collisions in your project. If not an issue, you can use `using`.

    julia> import Rrsp

Parse the file to create an instance object:

    julia> instance::Rrsp.RrspInstance = Rrsp.parseInstanceFromFile("data/instance.rrsp")

Use one of the solver functions to solve specific variant of the RRSP problem.
Here, we use function for solving RRSP with continuous budget in general digraphs:

    solution::Rrsp.RrspSolution = Rrsp.solveRrspContBudget(instance)

You can check if solving was succesful with `solution.value != Inf`.
To access the optimal paths use `solution.first_stage_path` and `solution.second_stage_path` objects.

## Index
```@index
```
