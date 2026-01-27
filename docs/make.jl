using SortedIntervals
using Documenter

DocMeta.setdocmeta!(SortedIntervals, :DocTestSetup, :(using SortedIntervals); recursive=true)

makedocs(;
    modules=[SortedIntervals],
    authors="Galen Lynch <galen@galenlynch.com>",
    sitename="SortedIntervals.jl",
    format=Documenter.HTML(;
        canonical="https://galenlynch.github.io/SortedIntervals.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/galenlynch/SortedIntervals.jl",
    devbranch="main",
)
