using SortedIntervals
using Test

@testset "SortedIntervals.jl" begin
    @testset "clipsize!" begin
        v = collect(1:10)
        ret = clipsize!(v, 3)
        @test ret === v
        @test v == [1, 2, 3]
        @test length(v) == 3
    end

    @testset "check_overlap" begin
        @test check_overlap(1, 3, 2, 3)
        @test !check_overlap(1, 3, 4, 5)
        @test check_overlap((1, 3), (2, 4))
        @test !check_overlap((1, 2), (3, 4))
        @test check_overlap([(1, 3), (2, 4)])
        @test !check_overlap([(1, 2), (3, 4)])
        # touching endpoints (closed intervals)
        @test check_overlap(1, 3, 3, 5)
        @test check_overlap((1, 3), (3, 5))
        # empty vector
        @test !check_overlap(NTuple{2,Int}[])
    end

    @testset "is_subinterval" begin
        @test is_subinterval(2, 3, 1, 5)
        @test is_subinterval((2, 3), (1, 5))
        @test !is_subinterval(0, 3, 1, 5)
        @test !is_subinterval((0, 3), (1, 5))
        # equal intervals
        @test is_subinterval((1, 5), (1, 5))
    end

    @testset "find_overlaps" begin
        ints = [(1.0, 3.0), (2.0, 4.0), (5.0, 6.0)]
        overlaps = find_overlaps(ints)
        @test 2 in overlaps[1]
        @test 1 in overlaps[2]
        @test isempty(overlaps[3])
        # single element
        @test find_overlaps([(1, 2)]) == [Int[]]
        # empty
        @test find_overlaps(NTuple{2,Int}[]) == Vector{Int}[]
    end

    @testset "find_all_overlapping" begin
        @test find_all_overlapping([(1, 3), (5, 7), (10, 12)], [(2, 6)]) == BitVector([1, 1, 0])
        # no overlaps
        @test find_all_overlapping([(1, 3), (5, 7)], [(8, 10)]) == BitVector([0, 0])
        # empty b
        @test find_all_overlapping([(1, 3), (5, 7)], NTuple{2,Int}[]) == BitVector([0, 0])
        # accessor function form
        data_a = [(t=(1, 3),), (t=(5, 7),)]
        data_b = [(t=(2, 6),)]
        @test find_all_overlapping(x -> x.t, x -> x.t, data_a, data_b) == BitVector([1, 1])
    end

    @testset "interval_intersect" begin
        @test interval_intersect(1, 3, 4, 5) === nothing
        @test interval_intersect(1, 4, 3, 5) == (3, 4)
        @test interval_intersect((1, 4), (3, 5)) == (3, 4)
        @test interval_intersect((1, 2), (3, 4)) === nothing
        # type promotion
        @test interval_intersect(1, 4, 3.0, 5.0) == (3.0, 4.0)
    end

    @testset "interval_intersect_measure" begin
        @test interval_intersect_measure(1, 4, 3, 5) == 1
        @test interval_intersect_measure((0, 10), (5, 15)) == 5
        @test interval_intersect_measure(1, 2, 3, 4) == 0
        # tuple form, no overlap
        @test interval_intersect_measure((1, 2), (5, 6)) == 0
    end

    @testset "interval_intersections" begin
        a = [(0.0, 5.0), (10.0, 15.0)]
        b = [(3.0, 12.0)]
        result = interval_intersections(a, b)
        @test result == [(3.0, 5.0), (10.0, 12.0)]
        # empty inputs
        @test interval_intersections(NTuple{2,Float64}[], NTuple{2,Float64}[]) == NTuple{2,Float64}[]
        # no overlap
        @test interval_intersections([(1.0, 2.0)], [(3.0, 4.0)]) == NTuple{2,Float64}[]
        # more intersections than max(na, nb)
        @test interval_intersections([(0, 3), (4, 7), (8, 11)], [(1, 5), (6, 9)]) ==
              [(1, 3), (4, 5), (6, 7), (8, 9)]
    end

    @testset "interval_intersections_overlapping" begin
        a = [(1.0, 5.0), (3.0, 7.0)]
        b = [(2.0, 10.0)]
        result = interval_intersections_overlapping(a, b)
        @test result == [(2.0, 7.0)]
        # non-overlapping inputs should match interval_intersections
        a2 = [(1.0, 3.0), (5.0, 7.0)]
        b2 = [(2.0, 6.0)]
        @test interval_intersections_overlapping(a2, b2) == interval_intersections(a2, b2)
    end

    @testset "overlap_interval_union" begin
        @test overlap_interval_union((1, 5), (3, 8)) == (1, 8)
        @test overlap_interval_union(1, 5, 3, 8) == (1, 8)
        # contained interval
        @test overlap_interval_union((1, 10), (3, 5)) == (1, 10)
    end

    @testset "intervals_diff" begin
        ints_a = [(1, 2), (3, 4)]
        @test intervals_diff(ints_a, [(1, 2)]) == [(3, 4)]
        @test intervals_diff(ints_a, [(3, 5)]) == [(1, 2)]
        @test intervals_diff(ints_a, NTuple{2,Int}[]) == ints_a
        @test intervals_diff(ints_a, [(0, 1)]) == ints_a
        @test intervals_diff(ints_a, [(4, 5)]) == ints_a
        @test intervals_diff(ints_a, [(1.2, 1.7)]) == [(1.0, 1.2), (1.7, 2.0), (3.0, 4.0)]
        @test intervals_diff(ints_a, [(1.2, 1.7), (3.2, 3.7)]) ==
              [(1.0, 1.2), (1.7, 2.0), (3.0, 3.2), (3.7, 4.0)]
        # empty a
        @test intervals_diff(NTuple{2,Int}[], [(1, 2)]) == NTuple{2,Int}[]
    end

    @testset "interval_complements" begin
        ints = [(2.0, 3.0), (5.0, 6.0)]
        @test interval_complements(0.0, 10.0, ints) == [(0.0, 2.0), (3.0, 5.0), (6.0, 10.0)]
        @test interval_complements(0.0, 10.0, NTuple{2,Float64}[]) == [(0.0, 10.0)]
        # with contraction
        @test interval_complements(0.0, 10.0, [(4.0, 6.0)], 1.0) == [(1.0, 3.0), (7.0, 9.0)]
    end

    @testset "join_intervals" begin
        ints = [(1.0, 2.0), (2.5, 3.0), (5.0, 6.0)]
        @test join_intervals(ints, 0.0) == ints
        @test join_intervals(ints, 1.0) == [(1.0, 3.0), (5.0, 6.0)]
        @test join_intervals(ints, 3.0) == [(1.0, 6.0)]
        # empty
        @test join_intervals(NTuple{2,Float64}[], 1.0) == NTuple{2,Float64}[]
        # single element
        @test join_intervals([(1.0, 2.0)], 0.0) == [(1.0, 2.0)]
    end

    @testset "join_intervals!" begin
        ints = [(1, 3), (4, 6), (10, 12)]
        ret = join_intervals!(ints, 1)
        @test ret === ints
        @test ints == [(1, 6), (10, 12)]
        # overlapping intervals: contained interval must not shrink the result
        @test join_intervals([(1, 10), (3, 5), (4, 6)], 0) == [(1, 10)]
    end

    @testset "expand_intervals" begin
        ints = [(1.0, 2.0), (5.0, 6.0)]
        expanded = expand_intervals(ints, 2.0)
        @test expanded == [(0.0, 3.0), (4.0, 7.0)]
        # original not mutated
        @test ints == [(1.0, 2.0), (5.0, 6.0)]
    end

    @testset "expand_intervals!" begin
        ints = [(2, 4), (8, 10)]
        ret = expand_intervals!(ints, 2)
        @test ret === ints
        @test ints == [(1, 5), (7, 11)]
    end

    @testset "throttle" begin
        points = [1.0, 1.1, 1.2, 5.0, 5.1]
        result = throttle(points, 0.5)
        @test result == [(1.0, 1.2), (5.0, 5.1)]
        # single point
        @test throttle([5.0], 1.0) == [(5.0, 5.0)]
        # empty
        @test throttle(Float64[], 1.0) == NTuple{2,Float64}[]
    end

    @testset "measure" begin
        @test measure((1.0, 5.0)) == 4.0
        @test measure((3, 7)) == 4
        @test measure(nothing) == 0
    end

    @testset "midpoint" begin
        @test midpoint((1.0, 5.0)) == 3.0
        @test midpoint((2, 6)) == 4.0
    end

    @testset "reduce_extrema" begin
        @test reduce_extrema(1, 3, 2, 5) == (1, 5)
        @test reduce_extrema((1, 3), (2, 5)) == (1, 5)
    end

    @testset "extrema_red" begin
        @test extrema_red([1.0, 2.0, 3.0]) == (1.0, 3.0)
        @test extrema_red([(1.0, 3.0), (2.0, 5.0)]) == (1.0, 5.0)
        # 2×N matrix method
        @test extrema_red([1.0 3.0; 5.0 8.0]) == (1.0, 8.0)
        # error on empty
        @test_throws ArgumentError extrema_red(NTuple{2,Float64}[])
    end

    @testset "clip_int" begin
        @test clip_int(0.0, 10.0, 2.0, 8.0) == (2.0, 8.0)
        @test clip_int((0.0, 10.0), (2.0, 8.0)) == (2.0, 8.0)
        @test clip_int(3.0, 7.0, 2.0, 8.0) == (3.0, 7.0)
    end

    @testset "intervals_are_ordered" begin
        @test intervals_are_ordered([(1.0, 2.0), (3.0, 4.0)])
        @test !intervals_are_ordered([(1.0, 5.0), (3.0, 4.0)])
        @test !intervals_are_ordered([(3.0, 4.0), (1.0, 2.0)])
        @test intervals_are_ordered(NTuple{2,Float64}[])
        # malformed interval (start > stop)
        @test !intervals_are_ordered([(5, 1)])
        # accessor function form
        @test intervals_are_ordered(x -> x.t, [(t=(1, 3),), (t=(4, 6),)])
        @test !intervals_are_ordered(x -> x.t, [(t=(1, 5),), (t=(3, 7),)])
    end

    @testset "intervals_are_partially_ordered" begin
        @test intervals_are_partially_ordered([(1.0, 5.0), (3.0, 7.0)])
        @test !intervals_are_partially_ordered([(3.0, 4.0), (1.0, 2.0)])
        # malformed interval
        @test !intervals_are_partially_ordered([(5, 1)])
        # accessor function form
        @test intervals_are_partially_ordered(x -> x.t, [(t=(1, 5),), (t=(3, 7),)])
    end

    @testset "interval_indices" begin
        basis = [1.0, 2.0, 3.0, 4.0, 5.0]
        ib, ie = interval_indices(basis, 2.0, 4.0)
        @test ib == 2
        @test ie == 4
        # tuple overload
        @test interval_indices(basis, (2.0, 4.0)) == (2, 4)
    end

    @testset "mask_events" begin
        events = [1.0, 2.0, 3.0, 4.0, 5.0]
        masked = mask_events(events, 2.0, 4.0)
        @test masked == [2.0, 3.0, 4.0]
        @test masked isa SubArray
    end

    @testset "maximum_interval_overlap" begin
        intervals = [(0.0, 5.0), (10.0, 20.0), (15.0, 25.0)]
        target = (12.0, 18.0)
        idx, overlap = maximum_interval_overlap(intervals, target)
        @test idx == 2
        @test overlap == 6.0
    end

    @testset "parse_ranges_str" begin
        @test parse_ranges_str("1-5") == [1, 2, 3, 4, 5]
        @test parse_ranges_str("1, 3, 5") == [1, 3, 5]
        @test parse_ranges_str("1-3, 7-9") == [1, 2, 3, 7, 8, 9]
        @test parse_ranges_str("") == Int[]
        # error on too many hyphens
        @test_throws ErrorException parse_ranges_str("1-2-3")
    end

    @testset "measure_to_bounds" begin
        @test measure_to_bounds(5, 3) == (5, 8)
        @test measure_to_bounds((5, 3)) == (5, 8)
        @test measure_to_bounds([(1, 2), (3, 4)]) == [(1, 3), (3, 7)]
        @test measure_to_bounds([1, 3], [2, 4]) == [(1, 3), (3, 7)]
    end

    @testset "clip_interval_duration" begin
        # shift right to fit
        @test clip_interval_duration(-2, 5, 0, 10) == (0, 7)
        # shift left to fit
        @test clip_interval_duration(8, 12, 0, 10) == (6, 10)
        # already fits
        @test clip_interval_duration(3, 7, 0, 10) == (3, 7)
        # interval larger than bounds
        @test clip_interval_duration(-5, 20, 0, 10) == (0, 10)
        # tuple form
        @test clip_interval_duration((8, 12), (0, 10)) == (6, 10)
        # 3-arg tuple form with nonzero boundmin
        @test clip_interval_duration((1, 4), 5, 20) == (5, 8)
    end
end
