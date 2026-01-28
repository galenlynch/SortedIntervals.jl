using SortedIntervals
using Test

@testset "SortedIntervals.jl" begin
    @testset "check_overlap" begin
        @test check_overlap(1, 3, 2, 3)
        @test !check_overlap(1, 3, 4, 5)
        @test check_overlap((1, 3), (2, 4))
        @test !check_overlap((1, 2), (3, 4))
        @test check_overlap([(1, 3), (2, 4)])  # vector with overlapping intervals
        @test !check_overlap([(1, 2), (3, 4)])  # vector with non-overlapping intervals
    end

    @testset "is_subinterval" begin
        @test is_subinterval(2, 3, 1, 5)
        @test is_subinterval((2, 3), (1, 5))
        @test !is_subinterval(0, 3, 1, 5)
        @test !is_subinterval((0, 3), (1, 5))
    end

    @testset "interval_intersect" begin
        @test interval_intersect(1, 3, 4, 5) === nothing
        @test interval_intersect(1, 4, 3, 5) == (3, 4)
        @test interval_intersect((1, 4), (3, 5)) == (3, 4)
        @test interval_intersect((1, 2), (3, 4)) === nothing
    end

    @testset "interval_intersect_measure" begin
        @test interval_intersect_measure(1, 4, 3, 5) == 1
        @test interval_intersect_measure((0, 10), (5, 15)) == 5
        @test interval_intersect_measure(1, 2, 3, 4) == 0
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
    end

    @testset "interval_intersections" begin
        a = [(0.0, 5.0), (10.0, 15.0)]
        b = [(3.0, 12.0)]
        result = interval_intersections(a, b)
        @test result == [(3.0, 5.0), (10.0, 12.0)]
    end

    @testset "interval_complements" begin
        ints = [(2.0, 3.0), (5.0, 6.0)]
        @test interval_complements(0.0, 10.0, ints) == [(0.0, 2.0), (3.0, 5.0), (6.0, 10.0)]
        @test interval_complements(0.0, 10.0, NTuple{2,Float64}[]) == [(0.0, 10.0)]
    end

    @testset "join_intervals" begin
        ints = [(1.0, 2.0), (2.5, 3.0), (5.0, 6.0)]
        @test join_intervals(ints, 0.0) == ints
        @test join_intervals(ints, 1.0) == [(1.0, 3.0), (5.0, 6.0)]
        @test join_intervals(ints, 3.0) == [(1.0, 6.0)]
    end

    @testset "expand_intervals" begin
        ints = [(1.0, 2.0), (5.0, 6.0)]
        expanded = expand_intervals(ints, 2.0)
        @test expanded == [(0.0, 3.0), (4.0, 7.0)]
    end

    @testset "throttle" begin
        points = [1.0, 1.1, 1.2, 5.0, 5.1]
        result = throttle(points, 0.5)
        @test result == [(1.0, 1.2), (5.0, 5.1)]
    end

    @testset "measure_midpoint" begin
        @test measure((1.0, 5.0)) == 4.0
        @test midpoint((1.0, 5.0)) == 3.0
        @test measure(nothing) == 0
    end

    @testset "clip_int" begin
        @test clip_int(0.0, 10.0, 2.0, 8.0) == (2.0, 8.0)
        @test clip_int((0.0, 10.0), (2.0, 8.0)) == (2.0, 8.0)
        @test clip_int(3.0, 7.0, 2.0, 8.0) == (3.0, 7.0)
    end

    @testset "intervals_are_ordered" begin
        @test intervals_are_ordered([(1.0, 2.0), (3.0, 4.0)])
        @test !intervals_are_ordered([(1.0, 5.0), (3.0, 4.0)])  # overlap
        @test !intervals_are_ordered([(3.0, 4.0), (1.0, 2.0)])  # not sorted
        @test intervals_are_ordered(NTuple{2,Float64}[])  # empty
    end

    @testset "intervals_are_partially_ordered" begin
        @test intervals_are_partially_ordered([(1.0, 5.0), (3.0, 7.0)])  # overlaps OK
        @test !intervals_are_partially_ordered([(3.0, 4.0), (1.0, 2.0)])  # not sorted
    end

    @testset "interval_indices" begin
        basis = [1.0, 2.0, 3.0, 4.0, 5.0]
        ib, ie = interval_indices(basis, 2.0, 4.0)
        @test ib == 2
        @test ie == 4
    end

    @testset "mask_events" begin
        events = [1.0, 2.0, 3.0, 4.0, 5.0]
        masked = mask_events(events, 2.0, 4.0)
        @test masked == [2.0, 3.0, 4.0]
    end

    @testset "parse_ranges_str" begin
        @test parse_ranges_str("1-5") == [1, 2, 3, 4, 5]
        @test parse_ranges_str("1, 3, 5") == [1, 3, 5]
        @test parse_ranges_str("1-3, 7-9") == [1, 2, 3, 7, 8, 9]
        @test parse_ranges_str("") == Int[]
    end

    @testset "reduce_extrema" begin
        @test reduce_extrema(1, 3, 2, 5) == (1, 5)
        @test reduce_extrema((1, 3), (2, 5)) == (1, 5)
    end

    @testset "extrema_red" begin
        @test extrema_red([1.0, 2.0, 3.0]) == (1.0, 3.0)
        @test extrema_red([(1.0, 3.0), (2.0, 5.0)]) == (1.0, 5.0)
    end

    @testset "find_overlaps" begin
        ints = [(1.0, 3.0), (2.0, 4.0), (5.0, 6.0)]
        overlaps = find_overlaps(ints)
        @test 2 in overlaps[1]  # interval 1 overlaps with interval 2
        @test 1 in overlaps[2]  # interval 2 overlaps with interval 1
        @test isempty(overlaps[3])  # interval 3 doesn't overlap with others
    end

    @testset "maximum_interval_overlap" begin
        intervals = [(0.0, 5.0), (10.0, 20.0), (15.0, 25.0)]
        target = (12.0, 18.0)
        idx, overlap = maximum_interval_overlap(intervals, target)
        @test idx == 2  # (10, 20) has most overlap with (12, 18)
        @test overlap == 6.0
    end
end
