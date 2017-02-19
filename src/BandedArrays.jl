"""Sparse array with a band of nonzeroes."""
immutable BandedArray{T} <: AbstractArray{T,2}
    data::Array{T,2}
    shape::Tuple{Int,Int}
    bandwidth::Int
    # offsets of band
    h_offset::Int
    v_offset::Int
    # limits for determining if [i, j] is in band
    lower::Int
    upper::Int

    function BandedArray(data::Array{T, 2}, shape::Tuple{Int, Int},
                         bandwidth::Int)
        if bandwidth < 1
            error("bandwidth must be positive")
        end
        drows = ndatarows(shape, bandwidth)
        nrows, ncols = shape
        if size(data) != (drows, ncols)
            error("data is wrong shape")
        end
        h_offset = max(ncols - nrows, 0)
        v_offset = max(nrows - ncols, 0)
        if ncols > nrows
            lower = nrows - ncols - bandwidth
            upper = bandwidth
        else
            lower = -bandwidth
            upper = nrows - ncols + bandwidth
        end
        return new(data, shape, bandwidth,
                   h_offset, v_offset,
                   lower, upper)
    end
end

function BandedArray(T::Type, shape::Tuple{Int,Int}, bandwidth::Int;
                     initialize::Bool=true)
    ndrows = ndatarows(shape, bandwidth)
    if initialize
        data = zeros(T, ndrows, shape[2])
    else
        data = Array(T, ndrows, shape[2])
    end
    return BandedArray{T}(data, shape, bandwidth)
end

function ndatarows(shape::Tuple{Int,Int}, bandwidth::Int)
    length_diff = abs(shape[1] - shape[2])
    return 2 * bandwidth + length_diff + 1
end

Base.size(A::BandedArray) = A.shape

"""The row in `data` which contains element [i, j]"""
function data_row(A::BandedArray, i::Int, j::Int)
    return (i - j) + A.h_offset + A.bandwidth + 1
end

"""Macro to help index into `data` matching A[i, j]"""
macro banddata(A, i, j)
    :($A.data[data_row($A, $i, $j), $j])
end

Base.getindex{T}(A::BandedArray{T}, i::Int, j::Int) = @banddata(A, i, j,)
Base.setindex!{T}(A::BandedArray{T}, v, i::Int, j::Int) = (@banddata(A, i, j) = v)

"""The rows in A's column `j` that are dense"""
function row_range(A::BandedArray, j::Int)
    start = max(1, j - A.h_offset - A.bandwidth)
    stop = min(j + A.v_offset + A.bandwidth, size(A)[1])
    return start, stop
end

"""The rows in data's column `j` that can be used"""
function data_row_range(A::BandedArray, j::Int)
    a, b = row_range(A, j)
    return data_row(A, a, j), data_row(A, b, j)
end

"""Get a view of the in-band elements of column `j`."""
function sparsecol(A::BandedArray, j::Int)
    start, stop = data_row_range(A, j)
    return view(A.data, start:stop, j)
end

"""Is [i, j] in the banded region?"""
function inband(A::BandedArray, i::Int, j::Int)
    nrows, ncols = size(A)
    if i < 1 || j < 1 || i > nrows || j > ncols
        return false
    end
    k = A.bandwidth
    return A.lower <= i- j <= A.upper
end

"""Return a dense representation of A"""
function Base.full{T}(A::BandedArray{T})
    result = zeros(T, size(A))
    for j = 1:size(A)[2]
        start, stop = row_range(A, j)
        dstart, dstop = data_row_range(A, j)
        result[start:stop, j] = A.data[dstart:dstop, j]
    end
    return result
end

"""Reverse the rows and the columns of `A`.

If `A` has shape `m` by `n`, then element [i, j] goes to
position [m - i, n - j].

"""
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
