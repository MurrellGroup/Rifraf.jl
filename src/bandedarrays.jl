# TODO: finish implementing AbstractArray interface.
# TODO: better bounds checking for getting and setting

"""Sparse array with a band of nonzeroes."""
mutable struct BandedArray{T} <: AbstractArray{T,2}
    data::Array{T,2}
    nrows::Int
    ncols::Int
    bandwidth::Int
    # offsets of band
    h_offset::Int  # how many more columns than rows
    v_offset::Int  # rows than columns
    # limits for determining if [i, j] is in band
    lower::Int
    upper::Int
    # used for initial creation and resizing
    row_padding::Int  # set by user
    col_padding::Int  # set by user
    default::T
    initialize::Bool

    function BandedArray{T}(data::Array{T,2}, shape::Tuple{Int,Int},
                            bandwidth::Int,
                            row_padding::Int, col_padding::Int,
                            default::T, initialize::Bool) where T
        if bandwidth < 1
            error("bandwidth must be positive")
        end
        nrows, ncols = shape
        drows = ndatarows(nrows, ncols, bandwidth)
        if size(data) != (drows + row_padding, ncols + col_padding)
            error("data is wrong shape")
        end
        h_offset = max(ncols - nrows, 0)
        v_offset = max(nrows - ncols, 0)
        lower, upper = bandlimits(nrows, ncols, bandwidth)
        new(data, nrows, ncols, bandwidth,
            h_offset, v_offset,
            lower, upper, row_padding, col_padding,
            default, initialize)
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

function BandedArray(T::Type, shape::Tuple{Int,Int}, bandwidth::Int;
                     row_padding::Int=0, col_padding::Int=0,
                     default=zero(T), initialize::Bool=true)
    nrows, ncols = shape
    data = allocate_data(T, nrows, ncols, bandwidth, row_padding, col_padding, initialize)
    return BandedArray{T}(data, shape, bandwidth, row_padding, col_padding, default, initialize)
end

function allocate_data(T::Type, nrows::Int, ncols::Int, bandwidth::Int,
                       row_padding::Int, col_padding::Int, initialize::Bool)
    ndrows = ndatarows(nrows, ncols, bandwidth)
    if initialize
        data = zeros(T, ndrows + row_padding, ncols + col_padding)
    else
        data = Array{T}(ndrows + row_padding, ncols + col_padding)
    end
    return data
end

function reallocate!(A::BandedArray{T}) where {T}
    # TODO: copy if requested
    A.data = allocate_data(T, A.nrows, A.ncols, A.bandwidth,
                           A.row_padding, A.col_padding, A.initialize)
end

function Base.resize!(A::BandedArray{T}, shape::Tuple{Int,Int}) where {T}
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

function newbandwidth!(A::BandedArray{T}, bandwidth::Int) where {T}
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

function Base.getindex(A::BandedArray{T}, i::Int, j::Int) where {T}
    if inband(A, i, j)
        return A.data[data_row(A, i, j), j]
    else
        return A.default
    end
end

function Base.setindex!(A::BandedArray{T}, v, i::Int, j::Int) where {T}
    if inband(A, i, j)
        A.data[data_row(A, i, j), j] = v
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
    return A.lower <= i - j <= A.upper
end

"""Return a dense representation of A"""
function Base.full(A::BandedArray{T}) where {T}
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
function flip!(A::BandedArray{T}) where {T}
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

"""
    equal_ranges(a_range, b_range)

Takes the true ranges of two sub-columns of a banded array and returns
the start:stop range for their overlapping elements.

# Examples
```julia-repl
julia> equal_ranges((3, 5), (4, 6))
((2, 3), (1, 2))

julia> equal_ranges((1, 5), (1, 2))
((1, 2), (1, 2))

julia> equal_ranges((1, 5), (4, 5))
((4, 5), (1, 2))

```

"""
function equal_ranges(a_range::Tuple{Int,Int},
                      b_range::Tuple{Int,Int})
    a_start, a_stop = a_range
    b_start, b_stop = b_range
    alen = a_stop - a_start + 1
    blen = b_stop - b_start + 1
    amin = max(b_start - a_start + 1, 1)
    amax = alen - max(a_stop - b_stop, 0)
    bmin = max(a_start - b_start + 1, 1)
    bmax = blen - max(b_stop - a_stop, 0)
    return (amin, amax), (bmin, bmax)
end
