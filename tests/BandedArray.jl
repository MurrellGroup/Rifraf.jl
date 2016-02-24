include("../BandedArray.jl")

using BandedArrayModule

using Base.Test

# test inband
begin
    m = BandedArray(Int32, (13, 11), 5)
    @test inband(m, 1, 1)
    @test inband(m, 8, 1)
    @test ! inband(m, 9, 1)
    @test inband(m, 1, 6)
    @test ! inband(m, 1, 7)
end

# test data_row
begin
    m = BandedArray(Int32, (3, 3), 1)
    @test data_row(m, 1, 1) == 2
    @test data_row(m, 2, 1) == 3
    @test data_row(m, 1, 2) == 1
    @test data_row(m, 2, 2) == 2
    @test data_row(m, 3, 2) == 3
    @test data_row(m, 2, 3) == 1
    @test data_row(m, 3, 3) == 2

    m = BandedArray(Int32, (3, 5), 0)
    @test data_row(m, 1, 1) == 3
    @test data_row(m, 3, 5) == 1

    m = BandedArray(Int32, (3, 5), 1)
    @test data_row(m, 1, 1) == 4
    @test data_row(m, 3, 5) == 2

    m = BandedArray(Int32, (5, 3), 1)
    @test data_row(m, 1, 1) == 2
    @test data_row(m, 5, 3) == 4
end

# test row_range
begin
    m = BandedArray(Int32, (3, 5), 0)
    @test row_range(m, 1) == (1, 1)
    @test row_range(m, 2) == (1, 2)
    @test row_range(m, 3) == (1, 3)
    @test row_range(m, 4) == (2, 3)
    @test row_range(m, 5) == (3, 3)

    m = BandedArray(Int32, (3, 5), 1)
    @test row_range(m, 1) == (1, 2)
    @test row_range(m, 2) == (1, 3)

    m = BandedArray(Int32, (5, 3), 1)
    @test row_range(m, 1) == (1, 4)
    @test row_range(m, 2) == (1, 5)
end

# test data_row_range
begin
    m = BandedArray(Int32, (3, 5), 0)
    @test data_row_range(m, 1) == (3, 3)
    @test data_row_range(m, 2) == (2, 3)
    @test data_row_range(m, 3) == (1, 3)
    @test data_row_range(m, 4) == (1, 2)
    @test data_row_range(m, 5) == (1, 1)

    m = BandedArray(Int32, (3, 5), 1)
    @test data_row_range(m, 1) == (4, 5)
    @test data_row_range(m, 2) == (3, 5)

    m = BandedArray(Int32, (5, 3), 1)
    @test data_row_range(m, 1) == (2, 5)
    @test data_row_range(m, 2) == (1, 5)
end

# test sparsecol
begin
    m = BandedArray(Int32, (5, 3), 1)
    m[1, 1] = 1
    @test sparsecol(m, 1) == [1, 0, 0, 0]
end

# test flip
begin
    m = BandedArray(Int32, (5, 3), 1)
    m[1, 1] = 1
    m = flip(m)
    @test m[5, 3] == 1
end

#test_sym
begin
    m = BandedArray(Int32, (3, 3), 0)
    m.data[:] = 1
    @test full(m) == eye(3)
end

#test_sym_band
begin
    m = BandedArray(Int32, (3, 3), 1)
    m.data[:] = 1
    expected = ones(Int32, (3, 3))
    expected[3, 1] = 0
    expected[1, 3] = 0
    @test full(m) == expected
end

# #test_wide
begin
    m = BandedArray(Int32, (3, 4), 0)
    m.data[:] = 1
    expected = zeros(Int32, (3, 4))
    expected[1, 1:2] = 1
    expected[2, 2:3] = 1
    expected[3, 3:4] = 1
    @test full(m) == expected
end

#test_wide_band
begin
    m = BandedArray(Int32, (3, 5), 0)
    m.data[:] = 1
    expected = zeros(Int32, (3, 5))
    expected[1, 1:3] = 1
    expected[2, 2:4] = 1
    expected[3, 3:5] = 1
    @test full(m) == expected
end

#test_wide_col
begin
    m = BandedArray(Int32, (3, 5), 1)
    m.data[:] = 1
    first = ones(Int32, 2)
    middle = ones(Int32, 3)
    last = first
    @test sparsecol(m, 1) == first
    @test sparsecol(m, 2) == middle
    @test sparsecol(m, 3) == middle
    @test sparsecol(m, 4) == middle
    @test sparsecol(m, 5) == last
end

#test_tall
begin
    m = BandedArray(Int32, (4, 3), 0)
    m.data[:] = 1
    expected = zeros(Int32, (4, 3))
    expected[1, 1:1] = 1
    expected[2, 1:2] = 1
    expected[3, 2:3] = 1
    expected[4, 3:3] = 1
    @test full(m) == expected
end

#test_tall_band
begin
    m = BandedArray(Int32, (5, 3), 1)
    m.data[:] = 1
    expected = ones(Int32, (5, 3))
    expected[5, 1] = 0
    expected[1, 3] = 0
    @test full(m) == expected
end

#test individual setting
begin
    m = BandedArray(Int32, (3, 3), 1)
    m[1, 2] = 3
    m[2, 1] = 5
    expected = zeros(Int32, (3, 3))
    expected[1, 2] = 3
    expected[2, 1] = 5
    @test full(m) == expected
end

#test individual setting of entire band
begin
    m = BandedArray(Int32, (3, 3), 1)
    m[1, 1] = 1
    m[2, 1] = 1
    m[1, 2] = 2
    m[2, 2] = 2
    m[3, 2] = 2
    m[2, 3] = 3
    m[3, 3] = 3
    expected = zeros(Int32, (3, 3))
    expected[1:2, 1] = 1
    expected[1:3, 2] = 2
    expected[2:3, 3] = 3
    @test full(m) == expected
end

#test individual setting of entire band, wide array
begin
    m = BandedArray(Int32, (3, 5), 1)
    m[1, 1] = 1
    m[1, 2] = 1
    m[1, 3] = 1
    m[3, 5] = 2
    expected = zeros(Int32, (3, 5))
    expected[1, 1:3] = 1
    expected[3, 5] = 2
    @test full(m) == expected
end

#test individual setting of entire band, tall array
begin
    m = BandedArray(Int32, (5, 3), 1)
    m[1, 1] = 1
    m[2, 1] = 1
    m[3, 1] = 1
    m[4, 1] = 1
    m[5, 3] = 2
    expected = zeros(Int32, (5, 3))
    expected[1:4, 1] = 1
    expected[5, 3] = 2
    @test full(m) == expected
end

# TODO: test bounds checks
