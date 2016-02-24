module BandedArrayModule

export BandedArray, sparsecol, data_row, row_range, data_row_range, inband, full, flip!

immutable BandedArray{T} <: AbstractArray{T,1}
    data::Array{T}
    shape::NTuple{2,Int}
    bandwidth::Int
    diff::Int
    h_offset::Int
    v_offset::Int

    function BandedArray(shape, bandwidth)
        nrows, ncols = shape
        diff = abs(nrows - ncols)
        data_nrows = 2 * bandwidth + diff + 1
        dshape = (data_nrows, ncols)
        h_offset = max(ncols - nrows, 0)
        v_offset = max(nrows - ncols, 0)

        return new(zeros(T, dshape), shape, bandwidth, diff, h_offset, v_offset)
    end
end

Base.size(A::BandedArray) = A.shape

Base.getindex{T}(A::BandedArray{T}, i::Int, j::Int) = get(A.data, (data_row(A, i, j), j), zero(T))

Base.setindex!{T}(A::BandedArray{T}, v, i::Int, j::Int) = (A.data[data_row(A, i, j), j] = v)

function sparsecol(A::BandedArray, j)
    start, stop = data_row_range(A, j)
    return A.data[start:stop, j]
end

# FIXME: these are wrong. need to be tested.
data_row(A::BandedArray, i, j) = (i - j) + A.h_offset + A.bandwidth + 1

function row_range(A::BandedArray, j)
    start = max(1, j - A.h_offset - A.bandwidth)
    stop = min(j + A.v_offset + A.bandwidth, size(A)[1])
    return start, stop
end

function data_row_range(A::BandedArray, j)
    start = max(1, A.h_offset + A.bandwidth - j)
    a, b = row_range(A, j)
    stop = start + (b - a)
    return start, stop
end

inband(A::BandedArray, i, j) = 1 <= data_row(A, i, j) <= size(A.data)[1]

function Base.full{T}(A::BandedArray{T})
    result = zeros(T, size(A))
    for j = 1:size(A)[2]
        start, stop = row_range(A, j)
        dstart, dstop = data_row_range(A, j)
        result[start:stop, j] = A.data[dstart:dstop, j]
    end
    return result
end

function flip{T}(A::BandedArray{T})
    newdata = flipdim(flipdim(A.data, 1), 2)
    return BandedArray{T}(newdata, A.shape, A.bandwidth, A.diff, A.h_offset, A.v_offset)
end

end
