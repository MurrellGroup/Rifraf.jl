# TODO: finish implementing AbstractArray interface.
# TODO: better bounds checking for getting and setting

"""Sparse array with a band of nonzeroes."""
type BandedArray{T} <: AbstractArray{T,2}
    data::Array{T,2}
    nrows::Int
    ncols::Int
    bandwidth::Int
    # offsets of band
    h_offset::Int
    v_offset::Int
    # limits for determining if [i, j] is in band
    lower::Int
    upper::Int
    # used for initial creation and resizing
    padding::Int  # set by user
    initialize::Bool

    default::T

    function BandedArray(data::Array{T, 2}, shape::Tuple{Int, Int},
                         bandwidth::Int, padding::Int,
                         initialize::Bool, default::T)
        if bandwidth < 1
            error("bandwidth must be positive")
        end
        nrows, ncols = shape
        drows = ndatarows(nrows, ncols, bandwidth)
        if size(data) != (drows + padding, ncols + padding)
            error("data is wrong shape")
        end
        h_offset = max(ncols - nrows, 0)
        v_offset = max(nrows - ncols, 0)
        lower, upper = bandlimits(nrows, ncols, bandwidth)
        return new(data, nrows, ncols, bandwidth,
                   h_offset, v_offset,
                   lower, upper, padding, initialize, default)
    end
end

function bandlimits(nrows, ncols, bandwidth)
    if ncols > nrows
        lower = nrows - ncols - bandwidth
        upper = bandwidth
    else
        lower = -bandwidth
        upper = nrows - ncols + bandwidth
    end
    return lower, upper
end

function BandedArray(T::Type, shape::Tuple{Int, Int}, bandwidth::Int;
                     padding::Int=0, initialize::Bool=true, default=zero(T))
    nrows, ncols = shape
    data = allocate_data(T, nrows, ncols, bandwidth, padding, initialize)
    return BandedArray{T}(data, shape, bandwidth, padding, initialize, default)
end

function allocate_data(T::Type,
                       nrows::Int, ncols::Int, bandwidth::Int,
                       padding::Int, initialize::Bool)
    ndrows = ndatarows(nrows, ncols, bandwidth)
    if initialize
        data = zeros(T, ndrows + padding, ncols + padding)
    else
        data = Array(T, ndrows + padding, ncols + padding)
    end
    return data
end

function reallocate!{T}(A::BandedArray{T})
    # TODO: copy if requested
    A.data = allocate_data(T, A.nrows, A.ncols, A.bandwidth, A.padding, A.initialize)
end

function Base.resize!{T}(A::BandedArray{T}, shape::Tuple{Int, Int})
    nrows, ncols = shape
    A.nrows = nrows
    A.ncols = ncols
    A.h_offset = max(ncols - nrows, 0)
    A.v_offset = max(nrows - ncols, 0)
    A.lower, A.upper = bandlimits(nrows, ncols, A.bandwidth)

    drows, dcols = size(A.data)
    new_ndrows = ndatarows(nrows, ncols, A.bandwidth)
    if new_ndrows > drows || ncols > dcols
        reallocate!(A)
    end
end

function newbandwidth!{T}(A::BandedArray{T}, bandwidth::Int)
    A.bandwidth = bandwidth
    reallocate!(A)
end

"""Number of used data rows"""
function ndatarows(nrows::Int, ncols::Int, bandwidth::Int)
    length_diff = abs(nrows - ncols)
    return 2 * bandwidth + length_diff + 1
end

Base.size(A::BandedArray) = (A.nrows, A.ncols)

"""The row in `data` which contains element [i, j]"""
function data_row(A::BandedArray, i::Int, j::Int)
    if !inband(A, i, j)
        error("[$i, $j] is not in band")
    end
    return (i - j) + A.h_offset + A.bandwidth + 1
end

"""Macro to help index into `data` matching A[i, j]"""
macro banddata(A, i, j)
    :($A.data[data_row($A, $i, $j), $j])
end

function Base.getindex{T}(A::BandedArray{T}, i::Int, j::Int)
    if inband(A, i, j)
        return @banddata(A, i, j)
    else
        return A.default
    end
end

function Base.setindex!{T}(A::BandedArray{T}, v, i::Int, j::Int)
    if inband(A, i, j)
        @banddata(A, i, j) = v
    else
        error("Cannot set out-of-band element [$i, $j].")
    end
end

"""The rows in A's column `j` that are dense"""
function row_range(A::BandedArray, j::Int)
    start = max(1, j - A.h_offset - A.bandwidth)
    stop = min(j + A.v_offset + A.bandwidth, A.nrows)
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
    if i < 1 || j < 1 || i > A.nrows || j > A.ncols
        return false
    end
    return A.lower <= i- j <= A.upper
end

"""Return a dense representation of A"""
function Base.full{T}(A::BandedArray{T})
    result = zeros(T, size(A))
    for j = 1:A.ncols
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
function flip!{T}(A::BandedArray{T})
    nrows = ndatarows(A.nrows, A.ncols, A.bandwidth)
    a, b = divrem(A.ncols, 2)
    for j in 1:a
        for i in 1:nrows
            first = A.data[i, j]
            second = A.data[nrows - i + 1, A.ncols - j + 1]
            A.data[i, j] = second
            A.data[nrows - i + 1, A.ncols - j + 1] = first
        end
    end
    # handle the middle column
    if b == 1
        c = div(nrows, 2)
        j = a + 1
        for i in 1:c
            first = A.data[i, j]
            second = A.data[nrows - i + 1, A.ncols - j + 1]
            A.data[i, j] = second
            A.data[nrows - i + 1, A.ncols - j + 1] = first
        end
    end
end
