# SortedIntervals [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://galenlynch.github.io/SortedIntervals.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://galenlynch.github.io/SortedIntervals.jl/dev/) [![Build Status](https://github.com/galenlynch/SortedIntervals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/galenlynch/SortedIntervals.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/galenlynch/SortedIntervals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/galenlynch/SortedIntervals.jl)

Efficient operations on sorted, non-overlapping intervals represented as tuples. Designed for temporal data analysis where intervals represent time ranges (trials, epochs, events).

## Overview

This package provides algorithms for interval arithmetic optimized for the common case where intervals are sorted by start time and mostly non-overlapping. Operations exploit sorted order for O(n) performance rather than O(n^2) naive approaches.

Intervals are represented as `NTuple{2, T}` (2-tuples) where `T` is typically a numeric type. The first element is the start, the second is the end (inclusive).

## Features

### Overlap Detection

```julia
using SortedIntervals

# Check if two intervals overlap
check_overlap(1.0, 3.0, 2.0, 4.0)  # true - intervals (1,3) and (2,4) overlap
check_overlap((1.0, 3.0), (5.0, 7.0))  # false

# Check if one interval contains another
is_subinterval((2.0, 3.0), (1.0, 5.0))  # true - (2,3) is inside (1,5)

# Find which intervals in list A overlap with any in list B
overlapping = find_all_overlapping(intervals_a, intervals_b)
```

### Intersection Operations

```julia
# Intersect two intervals
interval_intersect(1.0, 5.0, 3.0, 7.0)  # (3.0, 5.0)
interval_intersect((1.0, 5.0), (3.0, 7.0))  # (3.0, 5.0)

# Measure of intersection (overlap duration)
interval_intersect_measure((0.0, 10.0), (5.0, 15.0))  # 5.0

# Intersect two sorted lists of intervals
valid_epochs = interval_intersections(trial_times, artifact_free_times)
```

### Set Operations

```julia
# Set difference: intervals in A not covered by B
remaining = intervals_diff(all_times, excluded_times)

# Complement: gaps between intervals
gaps = interval_complements(0.0, 100.0, busy_intervals)

# Union of overlapping intervals
merged = overlap_interval_union((1.0, 5.0), (3.0, 7.0))  # (1.0, 7.0)
```

### Merging and Expansion

```julia
# Join intervals separated by small gaps
merged = join_intervals!(intervals, min_gap=0.1)

# Expand intervals and merge overlaps
expanded = expand_intervals(intervals, expansion=1.0)

# Convert points to intervals by joining nearby points
intervals = throttle([1.0, 1.1, 1.2, 5.0, 5.1], min_gap=0.5)
# Returns [(1.0, 1.2), (5.0, 5.1)]
```

### Validation

```julia
# Check intervals are sorted and non-overlapping
intervals_are_ordered([(1.0, 2.0), (3.0, 4.0)])  # true
intervals_are_ordered([(1.0, 5.0), (3.0, 4.0)])  # false - overlap

# Check intervals are sorted (overlaps allowed)
intervals_are_partially_ordered([(1.0, 5.0), (3.0, 7.0)])  # true
```

### Utilities

```julia
# Interval properties
measure((1.0, 5.0))   # 4.0 (duration/length)
midpoint((1.0, 5.0))  # 3.0

# Clipping
clip_int((0.0, 10.0), (2.0, 8.0))  # (2.0, 8.0)

# Find indices within an interval
ib, ie = interval_indices(sorted_times, start, stop)

# Mask events to interval
events_in_range = mask_events(event_times, start, stop)

# Parse range strings (useful for user input)
parse_ranges_str("1-5, 8, 10-12")  # [1, 2, 3, 4, 5, 8, 10, 11, 12]
```

## Relation to Other Packages

### vs IntervalSets.jl

[IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) provides rich interval types (`ClosedInterval`, `OpenInterval`, etc.) with a focus on type safety and abstract interval representations. SortedIntervals.jl differs:

- Uses plain tuples `(start, stop)` instead of special types
- Optimized for sorted collections of many intervals
- Focus on batch operations (intersect lists, find overlaps in bulk)
- Less type ceremony, more raw performance

Use IntervalSets.jl when you need: type-safe interval representations, dispatch on interval types, integration with other JuliaMath packages.

Use SortedIntervals.jl when you need: fast operations on large sorted interval lists, simple tuple-based representations, neuroscience/temporal analysis workflows.

### vs IntervalArithmetic.jl

[IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) implements rigorous interval arithmetic for numerical analysis with guaranteed bounds. SortedIntervals.jl is for discrete interval collections (time ranges, epochs) not numerical error propagation.

### vs IntervalTrees.jl

[IntervalTrees.jl](https://github.com/BioJulia/IntervalTrees.jl) provides tree-based data structures for interval queries with O(log n) lookup. SortedIntervals.jl assumes pre-sorted intervals and uses linear scans, which is often faster for:
- One-time operations on already-sorted data
- Sequential access patterns
- Avoiding tree construction overhead

Use IntervalTrees.jl for repeated queries on static interval sets.

## Design Philosophy

1. **Tuples over types**: Intervals are `NTuple{2}` for zero-overhead representation and easy interop with existing code.

2. **Sorted assumption**: Most operations assume intervals are sorted by start time. This enables O(n) algorithms instead of O(n^2).

3. **In-place when possible**: Functions like `join_intervals!` and `expand_intervals!` modify arrays in-place to reduce allocations.

4. **Neuroscience heritage**: Designed for trial-based experimental data where intervals represent stimulus presentations, behavioral epochs, or artifact periods.

## Installation

```julia
# Not yet registered - install from path or URL
using Pkg
Pkg.develop(path="/path/to/SortedIntervals")
```

## API Reference

### Overlap Detection
- `check_overlap(s1, e1, s2, e2)` / `check_overlap(a, b)` - Test interval overlap
- `is_subinterval(child, parent)` - Test containment
- `find_overlaps(intervals)` - Find all pairwise overlaps
- `find_all_overlapping(as, bs)` - Boolean mask of intervals in A overlapping B

### Intersection
- `interval_intersect(a, b)` - Intersection of two intervals (or nothing)
- `interval_intersect_measure(a, b)` - Overlap measure (0 if none)
- `interval_intersections(as, bs)` - All intersections between sorted lists
- `interval_intersections_overlapping(as, bs)` - Same, allowing overlapping inputs

### Set Operations
- `intervals_diff(as, bs)` - Set difference of interval lists
- `interval_complements(start, stop, intervals)` - Gaps between intervals
- `overlap_interval_union(a, b)` - Union of overlapping intervals

### Merging
- `join_intervals!(intervals, min_gap)` - Merge nearby intervals in-place
- `join_intervals(intervals, min_gap)` - Non-mutating version
- `expand_intervals!(intervals, expansion)` - Expand and merge
- `expand_intervals(intervals, expansion)` - Non-mutating version
- `throttle(points, min_gap)` - Convert points to intervals

### Validation
- `intervals_are_ordered(intervals)` - Check sorted and non-overlapping
- `intervals_are_partially_ordered(intervals)` - Check sorted (overlaps OK)

### Utilities
- `measure(interval)` - Length/duration of interval
- `midpoint(interval)` - Center point
- `clip_int(interval, bounds)` - Clip to bounds
- `interval_indices(basis, start, stop)` - Index range within interval
- `mask_events(times, start, stop)` - View of times within interval
- `maximum_interval_overlap(intervals, target)` - Best-matching interval
- `parse_ranges_str(s)` - Parse "1-5, 8" format
- `measure_to_bounds(start, duration)` - Convert (start, duration) to (start, end)
- `clip_interval_duration(req, bounds)` - Clip while preserving duration

### Extrema
- `reduce_extrema(a, b)` - Combined min/max of two intervals
- `extrema_red(intervals)` - Overall min/max across interval collection
- `clipsize!(vec, n)` - Resize vector and shrink capacity

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
