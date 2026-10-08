using Documenter
using Rrsp

makedocs(;
    sitename = "Robust Recoverable Shortest Path Solver",
    modules = [Rrsp],
    authors = "Marcel Jackiewicz",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://marceljackiewicz.github.io/rrsp/",
        repolink = "https://github.com/marceljackiewicz/rrsp",
        edit_link = "master",
        size_threshold_warn = 150 * 1024,
    ),
    pages = [
        "Home" => "index.md",
        "Reproducing the experiments" => "experiments.md",
    ],
    checkdocs = :exports,
    remotes = nothing,
)

# GitHub Actions deploys `docs/build` with actions/deploy-pages (see
# `.github/workflows/documentation.yml`). Local builds stop here; open
# `docs/build/index.html` in a browser.
