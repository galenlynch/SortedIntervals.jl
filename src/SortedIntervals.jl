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

clipsize!(a::AbstractVector, n::Integer) = sizehint!(resize!(a, n), n)

# =============================================================================
# Intervals
# =============================================================================

check_overlap(start1, stop1, start2, stop2) = (start1 <= stop2) & (start2 <= stop1)

function check_overlap(tupa::NTuple{2,<:Number}, tupb::NTuple{2,<:Number})
    check_overlap(tupa[1], tupa[2], tupb[1], tupb[2])
end

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

# Assumes sorted
function find_overlaps(a::AbstractVector{<:Tuple{<:Any,<:Any}})
    na = length(a)
    overlap_idx = Vector{Vector{Int}}(undef, na)
    @inbounds @simd for i = 1:na
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

function find_all_overlapping(fa, fb, intsa, intsb)
    intervals_are_ordered(fa, intsa) || error("intsa not valid")
    intervals_are_ordered(fb, intsb) || error("intsb not valid")
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

intervals_are_ordered(f, ints) = _intervals_are_ordered(f, well_ordered_crit, ints)
intervals_are_ordered(ints) = intervals_are_ordered(identity, ints)

intervals_are_partially_ordered(f, ints) =
    _intervals_are_ordered(f, partially_ordered_crit, ints)
intervals_are_partially_ordered(ints) = intervals_are_partially_ordered(identity, ints)

"""
Assumes each list is sorted and non-overlapping
"""
function interval_intersections(intsa, intsb)
    intervals_are_ordered(intsa) || error("intsa not valid")
    intervals_are_ordered(intsb) || error("intsb not valid")
    na = length(intsa)
    nb = length(intsb)
    outs = similar(intsa, max(na, nb))
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

function interval_intersections_overlapping(intsa, intsb)
    intervals_are_partially_ordered(intsa) || error("intsa not valid")
    intervals_are_partially_ordered(intsb) || error("intsb not valid")
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

@inline measure(a::NTuple{2,<:Number}) = a[2] - a[1]
@inline measure(::Nothing) = 0

@inline midpoint(a::NTuple{2,<:Number}) = (a[1] + a[2]) / 2

function reduce_extrema(s1::T, s2::T, t1::T, t2::T) where {T<:Number}
    return (min(s1, t1), max(s2, t2))
end
function reduce_extrema(s::NTuple{2,T}, t::NTuple{2,T}) where {T<:Number}
    return reduce_extrema(s..., t...)
end

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
    join_intervals!(ints::Vector{NTuple{2, <:Number}}, min_gap)

Join a list of sorted intervals, `ints`, if the gap between successive intervals
is less than `min_gap`. Mutates input in-place

Assumes ints are sorted by their first index.
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
        prev_end = int[2]
    end
    outno += 1
    ints[outno] = f((joined_start, prev_end))
    clipsize!(ints, outno)
    ints
end
join_intervals!(ints::AbstractVector, args...) = join_intervals!(identity, ints, args...)

"""
    join_intervals(f, ints::Vector{NTuple{2, <:Number}}, min_gap)

Like [`join_intervals!`](@ref), but does not mutate input.
"""
join_intervals(f, ints::AbstractVector, args...) = join_intervals!(f, copy(ints), args...)

join_intervals(ints::AbstractVector, args...) = join_intervals(identity, ints, args...)

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

function mask_events(event_times::AbstractVector{<:Number}, start, stop)
    i_b, i_e = interval_indices(event_times, start, stop)
    view(event_times, i_b:i_e)
end

"""
    interval_indices(
        basis::Union{<:AbstractVector, AbstractRange}, start::Number, stop::Number
    ) -> i_b, i_e

Find the indices in `basis` that correspond to the interval specified by `start`
 and `stop`.
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
    throttle(xs::AbstractVector{T}, min_gap::Number) where T<:Number ->
    Vector{NTuple{2, T}}

Join points in `x` into ranges, if the difference between neighboring elements
is less than `min_gap`. Assumes `x` is sorted.
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
    intervals_diff(ints_a, ints_b) -> ints_out

Find the 'set diff' of the intervals in `ints_a` and `ints_b`. Assumes both
inputs are sorted.
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

function expand_intervals!(ints_in::AbstractVector{<:NTuple{2,<:Number}}, expand::Number)
    half_exp = expand / 2
    f = ((b, e),) -> (b - half_exp, e + half_exp)
    join_intervals!(f, ints_in, expand)
    ints_in
end
expand_intervals(ints_in, expand) = expand_intervals!(copy(ints_in), expand)

"""
    parse_ranges_str(s::AbstractString)

Parse strings like "3", "1-5", or "1, 2-10, 12-30" into a vector of the integers
in specified ranges.
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

measure_to_bounds(a::Number, b::Number) = (a, a + b)
measure_to_bounds(t::NTuple{2}) = measure_to_bounds(t[1], t[2])
measure_to_bounds(ts::AbstractArray{<:NTuple{2}}) = measure_to_bounds.(ts)
measure_to_bounds(a::AbstractArray, b::AbstractArray) = measure_to_bounds.(a, b)

"Clip an interval while trying to maintain its duration"
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
    clip_interval_duration(int[1], int[2], boundmax)

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
