module SortedIntervals

export
    # Overlap detection
    check_overlap,
    is_subinterval,
    find_overlaps,
    find_all_overlapping,
    # Intersection
    interval_intersect,
    interval_intersect_measure,
    interval_intersections,
    interval_intersections_overlapping,
    # Set operations
    intervals_diff,
    interval_complements,
    overlap_interval_union,
    # Merge/expand
    join_intervals!,
    join_intervals,
    expand_intervals!,
    expand_intervals,
    # Utilities
    measure,
    midpoint,
    clip_int,
    relative_interval,
    throttle,
    mask_events,
    interval_indices,
    maximum_interval_overlap,
    # Validation
    intervals_are_ordered,
    intervals_are_partially_ordered,
    # Conversion
    parse_ranges_str,
    measure_to_bounds,
    clip_interval_duration,
    # Extrema
    reduce_extrema,
    extrema_red,
    # Internal (but useful)
    clipsize!

# =============================================================================
# Utility function
# =============================================================================

"""
    clipsize!(a::AbstractVector, n::Integer) -> AbstractVector

Resize `a` to length `n` and release excess memory capacity.
"""
clipsize!(a::AbstractVector, n::Integer) = sizehint!(resize!(a, n), n)

# =============================================================================
# Intervals
# =============================================================================

"""
    check_overlap(start1, stop1, start2, stop2) -> Bool
    check_overlap(a::NTuple{2,<:Number}, b::NTuple{2,<:Number}) -> Bool
    check_overlap(intervals::AbstractVector{<:NTuple{2}}) -> Bool

Test whether two intervals overlap, or whether any pair of intervals in a vector overlaps.

Intervals are treated as closed (endpoints are included).
"""
check_overlap(start1, stop1, start2, stop2) = (start1 <= stop2) & (start2 <= stop1)

function check_overlap(tupa::NTuple{2,<:Number}, tupb::NTuple{2,<:Number})
    check_overlap(tupa[1], tupa[2], tupb[1], tupb[2])
end

"""
    is_subinterval(startchild, stopchild, startparent, stopparent) -> Bool
    is_subinterval(child::NTuple{2,<:Number}, parent::NTuple{2,<:Number}) -> Bool

Test whether the child interval is entirely contained within the parent interval.
"""
function is_subinterval(startchild, stopchild, startparent, stopparent)
    (startchild >= startparent) & (stopchild <= stopparent)
end

function is_subinterval(tupa::NTuple{2,<:Number}, tupb::NTuple{2,<:Number})
    is_subinterval(tupa[1], tupa[2], tupb[1], tupb[2])
end

function check_overlap(a::AbstractVector{<:NTuple{2}})
    na = length(a)
    for i = 1:na, j = (i+1):na
        if check_overlap(a[i][1], a[i][2], a[j][1], a[j][2])
            return true
        end
    end
    false
end

"""
    find_overlaps(intervals::AbstractVector{<:Tuple{<:Any,<:Any}}) -> Vector{Vector{Int}}

Find all pairwise overlaps in a vector of intervals. Returns a vector of index lists, where
`result[i]` contains the indices of all intervals that overlap with interval `i`.

Assumes `intervals` is sorted by start time.

# Examples
```jldoctest
julia> find_overlaps([(1, 5), (3, 7), (8, 10)])
3-element Vector{Vector{Int64}}:
 [2]
 [1]
 []
```
"""
function find_overlaps(a::AbstractVector{<:Tuple{<:Any,<:Any}})
    na = length(a)
    overlap_idx = Vector{Vector{Int}}(undef, na)
    @inbounds for i = 1:na
        overlap_idx[i] = Vector{Int}()
    end
    @inbounds for i = 1:na
        thisstop = a[i][2]
        for j = (i+1):na
            a[j][1] > thisstop && break
            push!(overlap_idx[i], j)
            push!(overlap_idx[j], i)
        end
    end
    overlap_idx
end

"""
    find_all_overlapping(intsa, intsb) -> BitVector
    find_all_overlapping(fa, fb, intsa, intsb) -> BitVector

Return a `BitVector` of length `length(intsa)` indicating which intervals in `intsa` overlap
with any interval in `intsb`. Optional accessor functions `fa` and `fb` extract `(start, stop)`
tuples from each element.

Both inputs must be sorted and non-overlapping (see [`intervals_are_ordered`](@ref)).

# Examples
```jldoctest
julia> find_all_overlapping([(1, 3), (5, 7), (10, 12)], [(2, 6)])
3-element BitVector:
 1
 1
 0
```
"""
function find_all_overlapping(fa, fb, intsa, intsb)
    na = length(intsa)
    nb = length(intsb)
    outs = falses(na)
    ib = 1
    @inbounds for ia in eachindex(intsa)
        ab, ae = fa(intsa[ia])
        while ib <= nb && fb(intsb[ib])[2] <= ab
            ib += 1
        end
        ib > nb && break
        any_overlap = false
        for icheck = ib:nb
            bb, be = fb(intsb[icheck])
            bb < ae || break
            any_overlap = true
        end
        outs[ia] = any_overlap
    end
    outs
end

find_all_overlapping(intsa, intsb) = find_all_overlapping(identity, identity, intsa, intsb)

"""
    interval_intersect(start1, stop1, start2, stop2) -> NTuple{2} | Nothing
    interval_intersect(a::NTuple{2}, b::NTuple{2}) -> NTuple{2} | Nothing

Return the intersection of two intervals, or `nothing` if they do not overlap.
"""
function interval_intersect(start1::T, stop1::T, start2::T, stop2::T) where {T}
    ifelse(
        check_overlap(start1, stop1, start2, stop2),
        (max(start1, start2), min(stop1, stop2)),
        nothing,
    )
end

interval_intersect(b1, e1, b2, e2) = interval_intersect(promote(b1, e1, b2, e2)...)

function interval_intersect(a::NTuple{2}, b::NTuple{2})
    interval_intersect(a[1], a[2], b[1], b[2])
end

"""
    interval_intersect_measure(start1, stop1, start2, stop2) -> Number
    interval_intersect_measure(a::NTuple{2}, b::NTuple{2}) -> Number

Return the length of the intersection of two intervals, or zero if they do not overlap.
"""
function interval_intersect_measure(start1::T, stop1::T, start2::T, stop2::T) where {T}
    ifelse(
        check_overlap(start1, stop1, start2, stop2),
        min(stop1, stop2) - max(start1, start2),
        zero(T),
    )
end

interval_intersect_measure(b1, e1, b2, e2) =
    interval_intersect_measure(promote(b1, e1, b2, e2)...)

function interval_intersect_measure(a::NTuple{2}, b::NTuple{2})
    interval_intersect_measure(a[1], a[2], b[1], b[2])
end

@inline partially_ordered_crit(a, prevstart, _) = a >= prevstart
@inline well_ordered_crit(a, _, prevend) = a >= prevend

"""
Checks that each interval is well ordered, the set of intervals is sorted,
and the intervals are non-overlapping.
Returns a boolean.
"""
function _intervals_are_ordered(f, crit, ints)
    iter_result = iterate(ints)
    isnothing(iter_result) && return true
    element, state = iter_result
    prev_start, prev_end = f(element)
    ok = prev_start <= prev_end
    iter_result = iterate(ints, state)
    while ok & !isnothing(iter_result)
        (element, state) = iter_result
        a, b = f(element)
        ok &= a <= b
        ok &= crit(a, prev_start, prev_end)
        prev_start = a
        prev_end = b
        iter_result = iterate(ints, state)
    end
    ok
end

"""
    intervals_are_ordered(ints) -> Bool
    intervals_are_ordered(f, ints) -> Bool

Check that each interval is well-formed (`start <= stop`), intervals are sorted by start
time, and no intervals overlap. An optional accessor function `f` extracts `(start, stop)`
from each element.

See also [`intervals_are_partially_ordered`](@ref).

# Examples
```jldoctest
julia> intervals_are_ordered([(1, 3), (4, 6), (7, 9)])
true

julia> intervals_are_ordered([(1, 5), (3, 7)])
false
```
"""
intervals_are_ordered(f, ints) = _intervals_are_ordered(f, well_ordered_crit, ints)
intervals_are_ordered(ints) = intervals_are_ordered(identity, ints)

"""
    intervals_are_partially_ordered(ints) -> Bool
    intervals_are_partially_ordered(f, ints) -> Bool

Check that each interval is well-formed and intervals are sorted by start time. Unlike
[`intervals_are_ordered`](@ref), overlapping intervals are permitted.

# Examples
```jldoctest
julia> intervals_are_partially_ordered([(1, 5), (3, 7)])
true

julia> intervals_are_partially_ordered([(3, 7), (1, 5)])
false
```
"""
intervals_are_partially_ordered(f, ints) =
    _intervals_are_ordered(f, partially_ordered_crit, ints)
intervals_are_partially_ordered(ints) = intervals_are_partially_ordered(identity, ints)

"""
    interval_intersections(intsa, intsb) -> Vector{NTuple{2}}

Compute all pairwise intersections between two collections of intervals. Both inputs must
be sorted and non-overlapping (see [`intervals_are_ordered`](@ref)). Returns a sorted,
non-overlapping vector of intersection intervals.

See also [`interval_intersections_overlapping`](@ref) for inputs that may contain overlaps.

# Examples
```jldoctest
julia> interval_intersections([(1, 5), (7, 10)], [(3, 8)])
2-element Vector{Tuple{Int64, Int64}}:
 (3, 5)
 (7, 8)
```
"""
function interval_intersections(intsa, intsb)
    na = length(intsa)
    nb = length(intsb)
    outs = similar(intsa, na + nb)
    nout = 0
    ib = 1
    for (ab, ae) in intsa
        while ib <= nb && intsb[ib][2] <= ab
            ib += 1
        end
        ib > nb && break
        icheck = ib
        while icheck <= nb && intsb[icheck][1] < ae
            nout += 1
            outs[nout] = interval_intersect(ab, ae, intsb[icheck]...)
            icheck += 1
        end
    end
    clipsize!(outs, nout)
    outs
end

function flush_growing_intersect_intervals!(outs, nout, working, newminimum)
    i = 1
    new_nout = nout
    working_len = length(working)
    while i <= working_len
        if working[i][2] < newminimum
            new_nout += 1
            outs[new_nout] = working[i]
            deleteat!(working, i)
            working_len -= 1
        else
            i += 1
        end
    end
    new_nout
end

"""
    overlap_interval_union(a::NTuple{2}, b::NTuple{2}) -> NTuple{2}
    overlap_interval_union(ab, ae, bb, be) -> NTuple{2}

Return the bounding interval that covers both `a` and `b`. The inputs are assumed to
overlap; if they don't, the result spans the gap between them.
"""
overlap_interval_union(ab, ae, bb, be) = (min(ab, bb), max(ae, be))
overlap_interval_union(inta, intb) =
    overlap_interval_union(inta[1], inta[2], intb[1], intb[2])

function push_growing_intersect_interval!(working, newint)
    mergeno = 0
    i = 1
    working_int = newint
    working_len = length(working)
    while i <= working_len
        if check_overlap(working[i], working_int)
            working_int = overlap_interval_union(working[i], working_int)
            if mergeno == 0
                mergeno = i
                i += 1
            else
                deleteat!(working, i)
                working_len -= 1
            end
        else
            i += 1
        end
    end
    if mergeno > 0
        working[mergeno] = working_int
    else
        push!(working, working_int)
    end
end

"""
    interval_intersections_overlapping(intsa, intsb) -> Vector{NTuple{2}}

Like [`interval_intersections`](@ref), but allows the input interval lists to contain
overlapping intervals. Both inputs must still be sorted by start time
(see [`intervals_are_partially_ordered`](@ref)). Overlapping results are merged.
"""
function interval_intersections_overlapping(intsa, intsb)
    na = length(intsa)
    nb = length(intsb)
    outs = similar(intsa, na + nb)
    nout = 0
    ib = 1
    working_intersects = similar(intsa, 0)
    for (ab, ae) in intsa
        nout = flush_growing_intersect_intervals!(outs, nout, working_intersects, ab)
        while ib <= nb && intsb[ib][2] <= ab
            ib += 1
        end
        ib > nb && break
        icheck = ib
        while icheck <= nb && intsb[icheck][1] < ae
            bb, be = intsb[icheck]
            newint = interval_intersect(ab, ae, bb, be)
            push_growing_intersect_interval!(working_intersects, newint)
            icheck += 1
        end
    end
    for i in eachindex(working_intersects)
        nout += 1
        outs[nout] = working_intersects[i]
    end
    clipsize!(outs, nout)
    outs
end

"""
    measure(a::NTuple{2,<:Number}) -> Number
    measure(::Nothing) -> Int

Return the length (duration) of an interval `a`, i.e. `a[2] - a[1]`. Returns `0` for
`nothing`, which is useful when chained with [`interval_intersect`](@ref).
"""
@inline measure(a::NTuple{2,<:Number}) = a[2] - a[1]
@inline measure(::Nothing) = 0

"""
    midpoint(a::NTuple{2,<:Number}) -> Number

Return the midpoint of an interval `a`.
"""
@inline midpoint(a::NTuple{2,<:Number}) = (a[1] + a[2]) / 2

"""
    reduce_extrema(a::NTuple{2,T}, b::NTuple{2,T}) -> NTuple{2,T}
    reduce_extrema(s1, s2, t1, t2) -> NTuple{2}

Return `(min(a[1], b[1]), max(a[2], b[2]))` — the bounding interval that contains both
input intervals. Useful as a reduction operator over a collection of intervals.
"""
function reduce_extrema(s1::T, s2::T, t1::T, t2::T) where {T<:Number}
    return (min(s1, t1), max(s2, t2))
end
function reduce_extrema(s::NTuple{2,T}, t::NTuple{2,T}) where {T<:Number}
    return reduce_extrema(s..., t...)
end

"""
    extrema_red(a::AbstractVector{<:Number}) -> NTuple{2}
    extrema_red(a::AbstractArray{<:Number, 2}) -> NTuple{2}
    extrema_red(a::AbstractVector{<:NTuple{2,Number}}) -> NTuple{2}

Return the overall `(minimum, maximum)` across a collection of intervals or numbers.

For a vector of numbers, equivalent to `extrema`. For a `2×N` matrix, treats each column as
an interval. For a vector of 2-tuples, finds the global minimum start and maximum stop.
"""
extrema_red(a::AbstractVector{<:Number}) = extrema(a)

function extrema_red(a::AbstractArray{<:Number,2})
    na = size(a, 2)
    na > 0 || throw(ArgumentError("Collection must not be empty"))
    size(a, 1) == 2 || throw(ArgumentError("First dimension must be size 2"))
    cmin = a[1, 1]
    cmax = a[2, 1]
    for i = 2:na
        cmin = min(cmin, a[1, i])
        cmax = max(cmax, a[2, i])
    end
    return (cmin, cmax)
end

function extrema_red(a::A) where {T<:NTuple{2,Number},A<:AbstractVector{T}}
    na = length(a)
    na > 0 || throw(ArgumentError("Collection must not be empty"))
    cmin = a[1][1]
    cmax = a[1][2]
    for i = 2:na
        cmin = min(cmin, a[i][1])
        cmax = max(cmax, a[i][2])
    end
    return (cmin, cmax)
end

"""
    clip_int(int_begin, int_end, bound_begin, bound_end) -> NTuple{2}
    clip_int(input::NTuple{2,<:Number}, bounds::NTuple{2,<:Number}) -> NTuple{2}

Clamp an interval to lie within the given bounds. Each endpoint is clamped independently.
"""
function clip_int(
    int_begin::Number,
    int_end::Number,
    bound_begin::Number,
    bound_end::Number,
)
    (clamp(int_begin, bound_begin, bound_end), clamp(int_end, bound_begin, bound_end))
end
function clip_int(input::NTuple{2,<:Number}, bounds::NTuple{2,<:Number})
    clip_int(input..., bounds...)
end

"""
    relative_interval(interval, reference) -> NTuple{2}
    relative_interval(int_begin, int_end, ref_begin, ref_end) -> NTuple{2}

Express `interval` in coordinates relative to `reference`, clamped to
`(0, measure(reference))`. Equivalent to shifting `interval` so that
`ref_begin` maps to zero, then clipping to the reference duration.
"""
function relative_interval(int_begin::Number, int_end::Number, ref_begin::Number, ref_end::Number)
    clip_int(int_begin - ref_begin, int_end - ref_begin, 0, ref_end - ref_begin)
end
function relative_interval(interval::NTuple{2,<:Number}, reference::NTuple{2,<:Number})
    relative_interval(interval..., reference...)
end

"""
    join_intervals!(ints::AbstractVector{<:NTuple{2,<:Number}}, min_gap=0) -> AbstractVector
    join_intervals!(f, ints::AbstractVector{<:NTuple{2,<:Number}}, min_gap=0) -> AbstractVector

Merge successive intervals in the sorted vector `ints` when the gap between them is at most
`min_gap`. Mutates and resizes `ints` in-place. An optional function `f` is applied to each
merged interval before storing.

Assumes `ints` is sorted by start time. See [`join_intervals`](@ref) for a non-mutating version.

# Examples
```jldoctest
julia> ints = [(1, 3), (4, 6), (10, 12)];

julia> join_intervals!(ints, 1)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 6)
 (10, 12)
```
"""
function join_intervals! end

function join_intervals!(
    f::Function,
    ints::AbstractVector{<:NTuple{2,<:Number}},
    min_gap::Number = 0,
)
    nint = length(ints)
    if nint == 0
        clipsize!(ints, 0)
        return ints
    end
    outno = 0
    joined_start = ints[1][1]
    prev_end = ints[1][2]
    for intno = 2:nint
        int = ints[intno]
        if int[1] - prev_end > min_gap
            # End last stretch
            outno += 1
            ints[outno] = f((joined_start, prev_end))
            joined_start = int[1]
        end
        prev_end = max(prev_end, int[2])
    end
    outno += 1
    ints[outno] = f((joined_start, prev_end))
    clipsize!(ints, outno)
    ints
end
join_intervals!(ints::AbstractVector, args...) = join_intervals!(identity, ints, args...)

"""
    join_intervals(ints::AbstractVector{<:NTuple{2,<:Number}}, min_gap=0) -> Vector
    join_intervals(f, ints::AbstractVector{<:NTuple{2,<:Number}}, min_gap=0) -> Vector

Non-mutating version of [`join_intervals!`](@ref). Returns a new vector with successive
intervals merged when the gap between them is at most `min_gap`.

# Examples
```jldoctest
julia> join_intervals([(1, 3), (4, 6), (10, 12)], 1)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 6)
 (10, 12)
```
"""
join_intervals(f, ints::AbstractVector, args...) = join_intervals!(f, copy(ints), args...)

join_intervals(ints::AbstractVector, args...) = join_intervals(identity, ints, args...)

"""
    interval_complements(start, stop, intervals, contraction=0) -> Vector{NTuple{2}}

Find the gaps between `intervals` within the bounding interval `[start, stop]`. Each
complement boundary is contracted inward by `contraction` (useful for excluding edge
artifacts).

Assumes `intervals` is sorted and non-overlapping.

# Examples
```jldoctest
julia> interval_complements(0, 10, [(2, 4), (6, 8)])
3-element Vector{Tuple{Int64, Int64}}:
 (0, 2)
 (4, 6)
 (8, 10)
```
"""
function interval_complements(
    start::T,
    stop::T,
    intervals::AbstractVector{<:NTuple{2,T}},
    contraction::Number = 0,
) where {T}
    nint = length(intervals)
    if nint == 0
        if stop - start > 2 * contraction
            return NTuple{2,T}[(start + contraction, stop - contraction)]
        else
            return Vector{NTuple{2,T}}()
        end
    end
    complement = Vector{NTuple{2,T}}(undef, nint + 1)
    gapno = 0
    if intervals[1][1] - start > contraction
        gapno += 1
        complement[gapno] = (start + contraction, intervals[1][1] - contraction)
    end
    for i = 1:(nint-1)
        if intervals[i+1][1] - intervals[i][2] > 2 * contraction
            gapno += 1
            complement[gapno] =
                (intervals[i][2] + contraction, intervals[i+1][1] - contraction)
        end
    end
    if stop - intervals[end][2] > contraction
        gapno += 1
        complement[gapno] = (intervals[end][2] + contraction, stop - contraction)
    end
    clipsize!(complement, gapno)
    complement
end

function interval_complements(
    start,
    stop,
    intervals::AbstractVector{<:NTuple{2,T}},
    args...,
) where {T}
    interval_complements(convert(T, start), convert(T, stop), intervals, args...)
end

"""
    mask_events(event_times::AbstractVector{<:Number}, start, stop) -> SubArray

Return a view of `event_times` containing only elements within `[start, stop]`.
Uses binary search via [`interval_indices`](@ref). Assumes `event_times` is sorted.
"""
function mask_events(event_times::AbstractVector{<:Number}, start, stop)
    i_b, i_e = interval_indices(event_times, start, stop)
    view(event_times, i_b:i_e)
end

"""
    interval_indices(basis, start::Number, stop::Number) -> (i_b, i_e)
    interval_indices(basis, bounds::NTuple{2,<:Number}) -> (i_b, i_e)

Find the first and last indices in the sorted collection `basis` whose values fall within
`[start, stop]`. Uses binary search (`searchsortedfirst` / `searchsortedlast`).
"""
function interval_indices(
    basis::Union{<:AbstractVector,AbstractRange},
    start::Number,
    stop::Number,
)
    i_b = searchsortedfirst(basis, start)
    i_e = searchsortedlast(basis, stop)
    i_b, i_e
end
interval_indices(basis::AbstractVector, bnds::NTuple{2,<:Number}) =
    interval_indices(basis, bnds[1], bnds[2])

"""
    throttle(xs::AbstractVector{T}, min_gap::Number) -> Vector{NTuple{2,T}} where T<:Number

Group sorted points into intervals by merging consecutive points that are within `min_gap`
of each other. Each output interval spans from the first to last point in its group.

Assumes `xs` is sorted.

# Examples
```jldoctest
julia> throttle([1, 2, 3, 10, 11, 20], 2)
3-element Vector{Tuple{Int64, Int64}}:
 (1, 3)
 (10, 11)
 (20, 20)
```
"""
function throttle(xs::AbstractVector{T}, min_gap::Number) where {T<:Number}
    nx = length(xs)
    out = Vector{NTuple{2,T}}(undef, nx)
    nx == 0 && return out
    @inbounds joined_start = xs[1]
    last_x = joined_start
    nout = 0
    @inbounds for i = 2:nx
        x = xs[i]
        if x - last_x > min_gap
            nout += 1
            out[nout] = (joined_start, last_x)
            joined_start = x
        end
        last_x = x
    end
    nout += 1
    @inbounds out[nout] = (joined_start, last_x)
    clipsize!(out, nout)
    out
end

"""
    intervals_diff(ints_a, ints_b) -> Vector{NTuple{2}}

Compute the set difference of two sorted interval collections: the portions of `ints_a`
that are not covered by any interval in `ints_b`. Both inputs must be sorted and
non-overlapping.

# Examples
```jldoctest
julia> intervals_diff([(1, 10)], [(3, 5), (7, 8)])
3-element Vector{Tuple{Int64, Int64}}:
 (1, 3)
 (5, 7)
 (8, 10)
```
"""
function intervals_diff(
    ints_a::AbstractVector{<:NTuple{2,S}},
    ints_b::AbstractVector{<:NTuple{2,T}},
) where {S<:Number,T<:Number}
    na = length(ints_a)
    nb = length(ints_b)
    nb == 0 && return copy(ints_a)
    out_type = promote_type(S, T)
    ints_out = Vector{NTuple{2,out_type}}(undef, na + nb + 1)
    out_no = 0
    b_no = 1
    for a_no = 1:na
        # Advance b_no until overlap with the current a interval is possible
        while b_no <= nb && ints_b[b_no][2] <= ints_a[a_no][1]
            b_no += 1
        end
        # Break if no intervals in b could overlap with a
        if b_no > nb
            n_a_rest = na - a_no + 1
            ints_out[(out_no+1):(out_no+n_a_rest)] .=
                convert.(NTuple{2,out_type}, view(ints_a, a_no:na))
            out_no += n_a_rest
            break
        end
        last_b = b_no - 1 # Allow for no overlap by using b_no - 1
        # Find which b intervals overlap with this a interval
        while last_b < nb && check_overlap(ints_a[a_no], ints_b[last_b+1])
            last_b += 1
        end
        # Find complement between this a interval and overlapping b intervals
        complements = interval_complements(
            ints_a[a_no][1],
            ints_a[a_no][2],
            view(ints_b, b_no:last_b),
        )
        nc = length(complements)
        ints_out[(out_no+1):(out_no+nc)] = complements
        out_no += nc
        # Skip over used intervals in b
        b_no = ifelse(last_b > b_no, last_b, b_no)
    end
    clipsize!(ints_out, out_no)
    ints_out
end

"""
    expand_intervals!(ints_in::AbstractVector{<:NTuple{2,<:Number}}, expand::Number) -> AbstractVector

Expand each interval by `expand / 2` on each side, then merge any resulting overlaps.
Mutates and resizes `ints_in` in-place.

See [`expand_intervals`](@ref) for a non-mutating version.

# Examples
```jldoctest
julia> ints = [(2, 4), (8, 10)];

julia> expand_intervals!(ints, 2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 5)
 (7, 11)
```
"""
function expand_intervals!(ints_in::AbstractVector{<:NTuple{2,<:Number}}, expand::Number)
    half_exp = expand / 2
    f = ((b, e),) -> (b - half_exp, e + half_exp)
    join_intervals!(f, ints_in, expand)
    ints_in
end
"""
    expand_intervals(ints_in, expand::Number) -> Vector

Non-mutating version of [`expand_intervals!`](@ref).
"""
expand_intervals(ints_in, expand) = expand_intervals!(copy(ints_in), expand)

"""
    parse_ranges_str(s::AbstractString) -> Vector{Int}

Parse a comma-separated string of integers and integer ranges into a sorted vector of
unique integers. Ranges are specified with hyphens.

# Examples
```jldoctest
julia> parse_ranges_str("1-3, 7, 10-12")
7-element Vector{Int64}:
  1
  2
  3
  7
 10
 11
 12
```
"""
function parse_ranges_str(s::AbstractString)
    bnds = map(split(s, ',', keepempty = false)) do c
        parse.(Int, strip.(split(c, '-', keepempty = false)))
    end
    isempty(bnds) && return Int[]
    ranges = Set{Int}()
    for i in eachindex(bnds)
        nbnd = length(bnds[i])
        if nbnd == 1
            @inbounds push!(ranges, bnds[i][1])
        elseif nbnd == 2
            @inbounds b, e = extrema(bnds[i])
            for v = b:e
                push!(ranges, v)
            end
        else
            error("Too many hyphens")
        end
    end
    sort!(collect(ranges))
end

"""
    measure_to_bounds(start::Number, duration::Number) -> NTuple{2}
    measure_to_bounds(t::NTuple{2}) -> NTuple{2}
    measure_to_bounds(ts::AbstractArray{<:NTuple{2}}) -> AbstractArray
    measure_to_bounds(starts::AbstractArray, durations::AbstractArray) -> AbstractArray

Convert a `(start, duration)` representation to a `(start, stop)` representation, where
`stop = start + duration`.
"""
measure_to_bounds(a::Number, b::Number) = (a, a + b)
measure_to_bounds(t::NTuple{2}) = measure_to_bounds(t[1], t[2])
measure_to_bounds(ts::AbstractArray{<:NTuple{2}}) = measure_to_bounds.(ts)
measure_to_bounds(a::AbstractArray, b::AbstractArray) = measure_to_bounds.(a, b)

"""
    clip_interval_duration(reqb, reqe, boundmin, boundmax) -> NTuple{2}
    clip_interval_duration(int::NTuple{2}, bounds::NTuple{2}) -> NTuple{2}

Clip an interval to lie within `[boundmin, boundmax]` while attempting to preserve its
duration. The interval is shifted to fit within the bounds when possible. If the requested
interval is longer than the bounds, it is clamped to exactly `[boundmin, boundmax]`.

# Examples
```jldoctest
julia> clip_interval_duration(8, 12, 0, 10)
(6, 10)

julia> clip_interval_duration(-2, 5, 0, 10)
(0, 7)
```
"""
function clip_interval_duration(
    reqb::T,
    reqe::T,
    boundmin::T,
    boundmax::T,
) where {T<:Number}
    adj_b = max(zero(T), boundmin - reqb)
    adj_e = -max(zero(T), reqe - boundmax)
    adj = adj_b + adj_e
    clipped_b = reqb + adj
    clipped_e = reqe + adj
    req_int_smaller = reqe - reqb < boundmax - boundmin
    actual_b = ifelse(req_int_smaller, clipped_b, boundmin)
    actual_e = ifelse(req_int_smaller, clipped_e, boundmax)
    return (actual_b, actual_e)
end

clip_interval_duration(a, b, c, d) = clip_interval_duration(promote(a, b, c, d)...)
clip_interval_duration(reqb::Number, reqe, boundmax) =
    clip_interval_duration(reqb, reqe, 0, boundmax)
clip_interval_duration(int::NTuple{2}, intb::NTuple{2}) =
    clip_interval_duration(int..., intb...)
clip_interval_duration(int::NTuple{2}, boundmin, boundmax) =
    clip_interval_duration(int[1], int[2], boundmin, boundmax)

"""
    maximum_interval_overlap(xs::AbstractVector{NTuple{2,T}}, y::NTuple{2,T}) -> (index, overlap)

Find the interval in `xs` that has the greatest overlap with the target interval `y`.
Returns the index of the best-matching interval and the overlap measure. Returns `(0, typemin(T))`
if `xs` is empty.
"""
function maximum_interval_overlap(xs::AbstractVector{NTuple{2,T}}, y::NTuple{2,T}) where {T}
    best_ndx = 0
    best_overlap = typemin(T)
    for i in eachindex(xs)
        overlap = interval_intersect_measure(xs[i], y)
        better = overlap > best_overlap
        best_overlap = ifelse(better, overlap, best_overlap)
        best_ndx = ifelse(better, i, best_ndx)
    end
    best_ndx, best_overlap
end

end
