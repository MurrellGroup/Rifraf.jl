include("../BandedArray.jl")

using BandedArrayModule

using Base.Test

# test inband
begin
    m = BandedArray{Int32}((13, 11), 5)
    @test inband(m, 1, 1)
    @test inband(m, 8, 1)
    @test ! inband(m, 9, 1)
    @test inband(m, 1, 6)
    @test ! inband(m, 1, 7)
end

#test_sym
begin
    m = BandedArray{Int32}((3, 3), 0)
    m.data[:] = 1
    @test full(m) == eye(3)
end

#test_sym_band
begin
    m = BandedArray{Int32}((3, 3), 1)
    m.data[:] = 1
    expected = ones(Int32, (3, 3))
    expected[3, 1] = 0
    expected[1, 3] = 0
    @test full(m) == expected
end

# #test_wide
begin
    m = BandedArray{Int32}((3, 4), 0)
    m.data[:] = 1
    expected = zeros(Int32, (3, 4))
    expected[1, 1:2] = 1
    expected[2, 2:3] = 1
    expected[3, 3:4] = 1
    @test full(m) == expected
end

#test_wide_band
begin
    m = BandedArray{Int32}((3, 5), 0)
    m.data[:] = 1
    expected = zeros(Int32, (3, 5))
    expected[1, 1:3] = 1
    expected[2, 2:4] = 1
    expected[3, 3:5] = 1
    @test full(m) == expected
end

#test_wide_col
begin
    m = BandedArray{Int32}((3, 5), 1)
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
    m = BandedArray{Int32}((4, 3), 0)
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
    m = BandedArray{Int32}((5, 3), 1)
    m.data[:] = 1
    expected = ones(Int32, (5, 3))
    expected[5, 1] = 0
    expected[1, 3] = 0
    @test full(m) == expected
end

#test individual setting
begin
    m = BandedArray{Int32}((3, 3), 1)
    m[1, 2] = 3
    m[2, 1] = 5
    expected = zeros(Int32, (3, 3))
    expected[1, 2] = 3
    expected[2, 1] = 5
    @test full(m) == expected
end
