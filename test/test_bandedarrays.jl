using Base.Test

using Rifraf


@testset "BandedArrays" begin
    @testset "inband" begin
        m = BandedArray(Int, (13, 11), 5)
        @test inband(m, 1, 1)
        @test inband(m, 8, 1)
        @test !inband(m, 9, 1)
        @test inband(m, 1, 6)
        @test !inband(m, 1, 7)
    end

    @testset "data_row" begin
        m = BandedArray(Int, (3, 3), 1)
        @test data_row(m, 1, 1) == 2
        @test data_row(m, 2, 1) == 3
        @test data_row(m, 1, 2) == 1
        @test data_row(m, 2, 2) == 2
        @test data_row(m, 3, 2) == 3
        @test data_row(m, 2, 3) == 1
        @test data_row(m, 3, 3) == 2

        m = BandedArray(Int, (3, 5), 1)
        @test data_row(m, 1, 1) == 4
        @test data_row(m, 3, 5) == 2

        m = BandedArray(Int, (5, 3), 1)
        @test data_row(m, 1, 1) == 2
        @test data_row(m, 5, 3) == 4
    end

    @testset "row_range" begin
        m = BandedArray(Int, (3, 5), 1)
        @test row_range(m, 1) == (1, 2)
        @test row_range(m, 2) == (1, 3)

        m = BandedArray(Int, (5, 3), 1)
        @test row_range(m, 1) == (1, 4)
        @test row_range(m, 2) == (1, 5)
    end

    @testset "data_row_range" begin
        m = BandedArray(Int, (3, 5), 1)
        @test data_row_range(m, 1) == (4, 5)
        @test data_row_range(m, 2) == (3, 5)

        m = BandedArray(Int, (5, 3), 1)
        @test data_row_range(m, 1) == (2, 5)
        @test data_row_range(m, 2) == (1, 5)
    end

    @testset "sparsecol" begin
        m = BandedArray(Int, (5, 3), 1)
        m[1, 1] = 1
        @test sparsecol(m, 1) == [1, 0, 0, 0]
    end

    @testset "flip" begin
        m = BandedArray(Int, (5, 3), 1, row_padding=5, col_padding=5)
        m[1, 1] = 1
        flip!(m)
        @test m[5, 3] == 1
    end

    @testset "sym_band" begin
        m = BandedArray(Int, (3, 3), 1)
        m.data[:] = 1
        expected = ones(Int, (3, 3))
        expected[3, 1] = 0
        expected[1, 3] = 0
        @test full(m) == expected
    end

    @testset "test_wide" begin
        m = BandedArray(Int, (3, 4), 1)
        m.data[:] = 1
        expected = ones(Int, (3, 4))
        expected[1, end] = 0
        expected[end, 1] = 0
        @test full(m) == expected
    end

    @testset "wide col" begin
        m = BandedArray(Int, (3, 5), 1)
        m.data[:] = 1
        first = ones(Int, 2)
        middle = ones(Int, 3)
        last = first
        @test sparsecol(m, 1) == first
        @test sparsecol(m, 2) == middle
        @test sparsecol(m, 3) == middle
        @test sparsecol(m, 4) == middle
        @test sparsecol(m, 5) == last
    end

    @testset "tall" begin
        m = BandedArray(Int, (4, 3), 1)
        m.data[:] = 1
        expected = ones(Int, (4, 3))
        expected[1, end] = 0
        expected[end, 1] = 0
        @test full(m) == expected
    end

    @testset "tall band" begin
        m = BandedArray(Int, (5, 3), 1)
        m.data[:] = 1
        expected = ones(Int, (5, 3))
        expected[5, 1] = 0
        expected[1, 3] = 0
        @test full(m) == expected
    end

    @testset "individual setting" begin
        m = BandedArray(Int, (3, 3), 1)
        m[1, 2] = 3
        m[2, 1] = 5
        expected = zeros(Int, (3, 3))
        expected[1, 2] = 3
        expected[2, 1] = 5
        @test full(m) == expected
    end

    @testset "individual setting of entire band" begin
        m = BandedArray(Int, (3, 3), 1)
        m[1, 1] = 1
        m[2, 1] = 1
        m[1, 2] = 2
        m[2, 2] = 2
        m[3, 2] = 2
        m[2, 3] = 3
        m[3, 3] = 3
        expected = zeros(Int, (3, 3))
        expected[1:2, 1] = 1
        expected[1:3, 2] = 2
        expected[2:3, 3] = 3
        @test full(m) == expected
    end

    @testset "individual setting of entire band, wide array" begin
        m = BandedArray(Int, (3, 5), 1)
        m[1, 1] = 1
        m[1, 2] = 1
        m[1, 3] = 1
        m[3, 5] = 2
        expected = zeros(Int, (3, 5))
        expected[1, 1:3] = 1
        expected[3, 5] = 2
        @test full(m) == expected
    end
    @testset "individual setting of entire band, tall array" begin
        m = BandedArray(Int, (5, 3), 1)
        m[1, 1] = 1
        m[2, 1] = 1
        m[3, 1] = 1
        m[4, 1] = 1
        m[5, 3] = 2
        expected = zeros(Int, (5, 3))
        expected[1:4, 1] = 1
        expected[5, 3] = 2
        @test full(m) == expected
    end

    @testset "resize up within padding" begin
        m = BandedArray(Int, (3, 3), 1, row_padding=10, col_padding=10)
        old = m.data
        resize!(m, (5, 5))
        @test m.data == old
    end

    @testset "resize down" begin
        m = BandedArray(Int, (5, 5), 1)
        old = m.data
        resize!(m, (3, 3))
        @test m.data == old
    end

    @testset "resize with reallocate" begin
        m = BandedArray(Int, (3, 3), 1, row_padding=5, col_padding=5)
        @test Rifraf.row_range(m, 1) == (1, 2)
        old = m.data
        resize!(m, (5, 10))
        @test m.data != old
        @test Rifraf.row_range(m, 1) == (1, 2)
        @test Rifraf.row_range(m, 3) == (1, 4)
        @test Rifraf.row_range(m, 5) == (1, 5)
        @test Rifraf.row_range(m, 10) == (4, 5)
    end

    @testset "resize columns" begin
        m = BandedArray(Int, (3, 3), 1, row_padding=5, col_padding=10)
        @test Rifraf.row_range(m, 1) == (1, 2)
        old = m.data
        resize!(m, (5, 10))
        @test m.data == old
    end

    @testset "test get and set" begin
        default = Score(-Inf)
        m = BandedArray(Score, (13, 11), 5, default=default)
        m[1, 1] = Score(1.0)
        @test m[1, 1] == Score(1.0)
        @test m[end, 1] == default
    end

    @testset "equal_ranges" begin
        @test Rifraf.equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
        @test Rifraf.equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
        @test Rifraf.equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
    end

end
