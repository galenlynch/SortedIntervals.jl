```@meta
CurrentModule = SortedIntervals
```

# [Usage Guide](@id guide)

Intervals in SortedIntervals are plain `(start, stop)` tuples. Most functions that operate on collections of intervals assume the collection is sorted by start time. Functions that require non-overlapping input will validate this and throw an error if violated.

## Overlap Detection

Test whether intervals overlap, or find overlapping pairs within a collection.

[`check_overlap`](@ref) tests whether two intervals share any points (closed-interval semantics, so touching endpoints count as overlapping):

```julia
check_overlap((1, 5), (3, 8))   # true
check_overlap((1, 5), (5, 8))   # true  (touching)
check_overlap((1, 5), (6, 8))   # false
```

It also works on a vector, checking all pairs:

```julia
check_overlap([(1, 3), (2, 5), (7, 9)])  # true — (1,3) and (2,5) overlap
```

[`is_subinterval`](@ref) tests containment:

```julia
is_subinterval((2, 4), (1, 5))  # true — (2,4) is inside (1,5)
is_subinterval((2, 6), (1, 5))  # false
```

[`find_overlaps`](@ref) returns, for each interval, the indices of all other intervals it overlaps with. Assumes sorted input so it can stop early:

```julia
find_overlaps([(1, 5), (3, 7), (10, 12)])
# [[2], [1], []] — interval 1 overlaps interval 2 and vice versa
```

[`find_all_overlapping`](@ref) takes two sorted, non-overlapping interval lists and returns a `BitVector` marking which intervals in the first list overlap with any interval in the second:

```julia
find_all_overlapping([(1, 3), (5, 7), (10, 12)], [(2, 6)])
# BitVector([1, 1, 0])
```

## Intersections

Compute the overlapping region between intervals.

[`interval_intersect`](@ref) returns the intersection of two intervals, or `nothing` if they don't overlap:

```julia
interval_intersect((1, 5), (3, 8))  # (3, 5)
interval_intersect((1, 2), (3, 4))  # nothing
```

[`interval_intersect_measure`](@ref) returns just the length of the overlap:

```julia
interval_intersect_measure((0, 10), (5, 15))  # 5
interval_intersect_measure((1, 2), (3, 4))    # 0
```

[`interval_intersections`](@ref) computes all pairwise intersections between two sorted, non-overlapping lists:

```julia
interval_intersections([(1, 5), (7, 10)], [(3, 8)])
# [(3, 5), (7, 8)]
```

[`interval_intersections_overlapping`](@ref) does the same but allows the input lists to contain overlapping intervals (they must still be sorted by start time). Overlapping results are merged:

```julia
interval_intersections_overlapping([(1, 5), (3, 7)], [(2, 10)])
# [(2, 7)]
```

## Set Operations

[`intervals_diff`](@ref) computes the set difference — the portions of the first list not covered by the second:

```julia
intervals_diff([(1, 10)], [(3, 5), (7, 8)])
# [(1, 3), (5, 7), (8, 10)]
```

[`interval_complements`](@ref) finds the gaps between intervals within a bounding range. An optional `contraction` parameter shrinks each gap boundary inward:

```julia
interval_complements(0, 10, [(2, 4), (6, 8)])
# [(0, 2), (4, 6), (8, 10)]

interval_complements(0.0, 10.0, [(4.0, 6.0)], 1.0)
# [(1.0, 3.0), (7.0, 9.0)]
```

[`overlap_interval_union`](@ref) returns the bounding interval that covers both inputs:

```julia
overlap_interval_union((1, 5), (3, 8))  # (1, 8)
```

## Merging and Expansion

[`join_intervals`](@ref) merges successive intervals when the gap between them is at most `min_gap`:

```julia
join_intervals([(1, 3), (4, 6), (10, 12)], 1)
# [(1, 6), (10, 12)]
```

[`join_intervals!`](@ref) does the same in-place, resizing the input vector.

[`expand_intervals`](@ref) widens each interval by `expand / 2` on each side, then merges any resulting overlaps:

```julia
expand_intervals([(2, 4), (8, 10)], 2)
# [(1, 5), (7, 11)]
```

[`expand_intervals!`](@ref) does the same in-place.

[`throttle`](@ref) converts a sorted vector of points into intervals by grouping consecutive points that are within `min_gap` of each other:

```julia
throttle([1, 2, 3, 10, 11, 20], 2)
# [(1, 3), (10, 11), (20, 20)]
```

## Validation

[`intervals_are_ordered`](@ref) checks that intervals are well-formed (`start <= stop`), sorted by start time, and non-overlapping:

```julia
intervals_are_ordered([(1, 3), (4, 6)])  # true
intervals_are_ordered([(1, 5), (3, 7)])  # false — overlapping
intervals_are_ordered([(5, 1)])          # false — malformed
```

[`intervals_are_partially_ordered`](@ref) is the same but allows overlaps:

```julia
intervals_are_partially_ordered([(1, 5), (3, 7)])  # true
```

Both accept an accessor function for working with non-tuple collections:

```julia
intervals_are_ordered(x -> x.time, my_data)
```

## Utilities

[`measure`](@ref) and [`midpoint`](@ref) compute basic interval properties:

```julia
measure((3, 7))   # 4
midpoint((2, 6))  # 4.0
measure(nothing)  # 0 — useful after interval_intersect
```

[`clip_int`](@ref) clamps an interval to lie within bounds (each endpoint clamped independently):

```julia
clip_int((1, 10), (3, 8))  # (3, 8)
```

[`clip_interval_duration`](@ref) clips an interval to bounds while attempting to preserve its duration by shifting:

```julia
clip_interval_duration(8, 12, 0, 10)  # (6, 10) — shifted left to fit
clip_interval_duration(-2, 5, 0, 10)  # (0, 7)  — shifted right to fit
```

[`interval_indices`](@ref) finds the index range in a sorted vector that falls within an interval, using binary search:

```julia
interval_indices([10, 20, 30, 40, 50], 15, 35)  # (2, 3)
```

[`mask_events`](@ref) returns a view of a sorted vector containing only elements within a range:

```julia
mask_events([1, 3, 5, 7, 9], 2, 6)  # [3, 5]
```

[`maximum_interval_overlap`](@ref) finds which interval in a collection has the greatest overlap with a target:

```julia
maximum_interval_overlap([(1, 5), (4, 10), (12, 15)], (3, 8))
# (2, 4) — index 2, overlap of 4
```

[`measure_to_bounds`](@ref) converts `(start, duration)` to `(start, stop)`:

```julia
measure_to_bounds(5, 3)  # (5, 8)
```

[`reduce_extrema`](@ref) and [`extrema_red`](@ref) compute bounding ranges:

```julia
reduce_extrema((1, 5), (3, 8))              # (1, 8)
extrema_red([(1, 5), (3, 8), (6, 7)])       # (1, 8)
```

[`parse_ranges_str`](@ref) parses human-readable range strings:

```julia
parse_ranges_str("1-3, 7, 10-12")  # [1, 2, 3, 7, 10, 11, 12]
```

[`clipsize!`](@ref) resizes a vector and releases excess memory capacity:

```julia
v = Vector{Int}(undef, 100)
clipsize!(v, 3)  # now length 3 with no excess allocation
```
