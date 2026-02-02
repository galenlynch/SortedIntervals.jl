# SortedIntervals [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://galenlynch.github.io/SortedIntervals.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://galenlynch.github.io/SortedIntervals.jl/dev/) [![Build Status](https://github.com/galenlynch/SortedIntervals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/galenlynch/SortedIntervals.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/galenlynch/SortedIntervals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/galenlynch/SortedIntervals.jl)

Fast set operations on sorted, non-overlapping intervals represented as plain tuples. Exploits sorted order for O(n) performance on intersection, difference, complement, and merge operations.

## Quick Example

```julia
using SortedIntervals

trials = [(0.0, 5.0), (10.0, 15.0), (20.0, 25.0)]
artifacts = [(3.0, 12.0), (22.0, 23.0)]

# Find clean portions of each trial
clean = intervals_diff(trials, artifacts)
# [(0.0, 3.0), (12.0, 15.0), (20.0, 22.0), (23.0, 25.0)]

# Merge intervals separated by small gaps
joined = join_intervals(clean, 1.0)

# Find which trials overlap with artifacts
affected = find_all_overlapping(trials, artifacts)
# BitVector([1, 1, 1])
```

## When to Use What

| Package | Representation | Best for |
|---|---|---|
| **SortedIntervals.jl** | `(start, stop)` tuples | Batch operations on sorted interval lists |
| [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) | Typed intervals (`ClosedInterval`, etc.) | Type-safe single-interval operations |
| [IntervalTrees.jl](https://github.com/BioJulia/IntervalTrees.jl) | Tree data structure | Repeated queries on static interval sets |
| [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) | Rigorous interval arithmetic | Numerical error propagation with guaranteed bounds |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/galenlynch/SortedIntervals.jl")
```

## Documentation

See the [full documentation](https://galenlynch.github.io/SortedIntervals.jl/dev/) for a usage guide and API reference.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
