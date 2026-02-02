```@meta
CurrentModule = SortedIntervals
```

# SortedIntervals.jl

SortedIntervals provides fast set operations on sorted, non-overlapping intervals represented as plain `NTuple{2}` tuples. By assuming intervals are sorted by start time, most operations run in O(n) time using linear scans rather than O(n^2) pairwise comparisons.

## Design Principles

**Tuples, not types.** Intervals are `(start, stop)` tuples. There is no `Interval` struct. This means zero overhead, easy interop with existing code, and no need to convert between representations.

**Sorted-input assumption.** Functions that operate on collections of intervals assume they are sorted by start time. This is what makes O(n) algorithms possible. Use [`intervals_are_ordered`](@ref) to validate, or [`intervals_are_partially_ordered`](@ref) if overlaps are expected.

**In-place variants.** Functions that resize their output (like [`join_intervals!`](@ref) and [`expand_intervals!`](@ref)) have mutating versions that modify arrays in-place to reduce allocations. Non-mutating versions copy first.

**Accessor functions.** Many functions accept an optional function `f` that extracts `(start, stop)` from each element, so you can operate on collections of structs or named tuples without converting them first:
```julia
# Works on any collection, as long as f returns (start, stop)
intervals_are_ordered(x -> x.time_range, my_trials)
```

## Example Workflow

A typical workflow for temporal data analysis (e.g., neuroscience trial data):

```julia
using SortedIntervals

# Define trial windows and artifact periods
trials = [(0.0, 5.0), (10.0, 15.0), (20.0, 25.0)]
artifacts = [(3.0, 12.0), (22.0, 23.0)]

# Validate inputs
intervals_are_ordered(trials)     # true
intervals_are_ordered(artifacts)  # true

# Remove artifact periods from trials
clean = intervals_diff(trials, artifacts)
# [(0.0, 3.0), (12.0, 15.0), (20.0, 22.0), (23.0, 25.0)]

# Merge fragments separated by small gaps
joined = join_intervals(clean, 1.5)
# [(0.0, 3.0), (10.0, 15.0), (20.0, 25.0)]

# Find gaps in the clean data within a recording window
gaps = interval_complements(0.0, 30.0, joined)

# Select events that fall within clean periods
event_times = [1.0, 4.0, 11.0, 21.0, 28.0]
events_in_first = mask_events(event_times, joined[1]...)
```

See the [Usage Guide](@ref guide) for detailed examples of each function category, or the [API Reference](@ref api) for complete function signatures.
