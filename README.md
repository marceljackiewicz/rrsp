# Robust Recoverable Shortest Path Solver
Authors: Marcel Jackiewicz, Adam Kasperski, Paweł Zieliński

This repository contains implementation of algorithms, in the form of [Julia](https://julialang.org/) package, for solving RRSP problem described in the following articles

- Recoverable Robust Shortest Path Problem Under Interval Budgeted Uncertainty Representations, M. Jackiewicz, A. Kasperski and P. Zieliński, Networks 85 (1), 127-141, [Networks](https://onlinelibrary.wiley.com/doi/abs/10.1002/net.22255)

- Computational Complexity of the Recoverable Robust Shortest Path Problem with Discrete Recourse, M. Jackiewicz, A. Kasperski and P. Zieliński, arXiv preprint arXiv:2403.20000, [arxiv](https://arxiv.org/abs/2403.20000)

- Computational complexity of the recoverable robust shortest path problem in acyclic digraphs, A. Kasperski and P. Zieliński, arXiv preprint arXiv:2410.09425, [arxiv](https://arxiv.org/abs/2410.09425)

## Documentation

The package documentation is generated using [Documenter.jl](https://documenter.juliadocs.org/stable/) package.
We opted for not keeping script generated content in the repository (aside from package configuration) and since bitbucket doesn't allow hosting web pages currently, you can generate the documentation locally.
To do so, clone the repository, enter `Rrsp/docs` directory and use

    $ julia --project make.jl

from your command line. Make sure `Documenter` package is installed. (You can obtain it via `Pkg` package manager.)
The documentation is available as a web page at `Rrsp/docs/build/index.html`. Use any web browser to open it.

## Project Funding Notice

The authors were supported by the National Science Centre, Poland, grant 2022/45/B/HS4/00355.
