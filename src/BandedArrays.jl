module BandedArrays

export BandedArray, sparsecol, data_row, row_range, data_row_range, inband, full, flip

immutable BandedArray{T} <: AbstractArray{T,2}
    data::Array{T,2}
    shape::Tuple{Int,Int}
    bandwidth::Int
    h_offset::Int
    v_offset::Int

    function BandedArray(data, shape, bandwidth)
        # TODO: check bandwidth
        drows = datarows(shape, bandwidth)
        nrows, ncols = shape
        h_offset = max(ncols - nrows, 0)
        v_offset = max(nrows - ncols, 0)

        if size(data) != (drows, ncols)
            error("data is wrong shape")
        end

        return new(data, shape, bandwidth, h_offset, v_offset)
    end

end

function BandedArray(T::Type, shape::Tuple{Int,Int}, bandwidth::Int)
    drows = datarows(shape, bandwidth)
    dshape = (drows, shape[2])
    data = zeros(T, dshape)
    return BandedArray{T}(data, shape, bandwidth)
end

Base.size(A::BandedArray) = A.shape

function datarows(shape::Tuple{Int,Int}, bandwidth::Int)
    diff = abs(shape[1] - shape[2])
    return 2 * bandwidth + diff + 1
end

Base.getindex{T}(A::BandedArray{T}, i::Int, j::Int) = A.data[data_row(A, i, j), j]

Base.setindex!{T}(A::BandedArray{T}, v, i::Int, j::Int) = (A.data[data_row(A, i, j), j] = v)

function sparsecol(A::BandedArray, j::Int)
    start, stop = data_row_range(A, j)
    return sub(A.data, start:stop, j)
end

# TODO: bounds checks
data_row(A::BandedArray, i::Int, j::Int) = (i - j) + A.h_offset + A.bandwidth + 1

function row_range(A::BandedArray, j::Int)
    start = max(1, j - A.h_offset - A.bandwidth)
    stop = min(j + A.v_offset + A.bandwidth, size(A)[1])
    return start, stop
end

function data_row_range(A::BandedArray, j::Int)
    a, b = row_range(A, j)
    return data_row(A, a, j), data_row(A, b, j)
end

function inband(A::BandedArray, i::Int, j::Int)
    nrows, ncols = size(A.data)
    if i < 1 || j < 1 || i > nrows || j > ncols
        return false
    end
    return 1 <= data_row(A, i, j) <= size(A.data)[1]::Int
end

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
    newdata = similar(A.data)
    nrows, ncols = size(A.data)
    for j in 1:ncols
        for i in 1:nrows
            newdata[i, j] = A.data[nrows - i + 1, ncols - j + 1]
        end
    end
    return BandedArray{T}(newdata, size(A), A.bandwidth)
end

end
