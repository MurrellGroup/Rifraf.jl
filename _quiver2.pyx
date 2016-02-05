import cython
cimport cython
import numpy as np
cimport numpy as np

class OutsideBandException(Exception):
    pass


cdef class BandedMatrix(object):
    """A sparse banded matrix.

    Currently only supports limited slicing and no arithmetic
    operations.

    """

    cdef int nrows
    cdef int ncols
    cdef int bandwidth
    cdef public int offset
    cdef public float[:, :] data

    def __init__(self, shape, bandwidth, offset=0):
        nrows, ncols = shape
        if not nrows > 0 and ncols > 0:
            raise Exception()
        self.nrows = nrows
        self.ncols = ncols
        self.bandwidth = bandwidth
        self.data = np.zeros((2 * bandwidth + 1, ncols), dtype=np.float32)
        self.offset = offset

    property shape:
        def __get__(self):
            return (self.nrows, self.ncols)

    cdef _convert_row(self, int i, int j, int o=0):
        if i < max(j - self.bandwidth - self.offset, 0):
            raise OutsideBandException()
        if i > min(j + self.bandwidth - self.offset + o, self.shape[0] - (1 - o)):
            raise OutsideBandException()
        return i - (j - self.bandwidth - self.offset)

    cpdef range(self, j):
        start = max(j - self.bandwidth - self.offset, 0)
        stop = min(j - self.offset + self.bandwidth + 1, self.shape[0])
        assert start < stop
        return start, stop

    cpdef data_range(self, j):
        start, stop = self.range(j)
        dstart = self._convert_row(start, j)
        dstop = self._convert_row(stop, j, o=1)
        assert dstart < dstop
        return dstart, dstop

    def __getitem__(self, key):
        # TODO: only supports limited slicing
        cdef int j
        cdef int row
        i, j = key
        if not isinstance(j, int):
            raise Exception()
        if j < 0:
            j = self.shape[1] + j
        if isinstance(i, int):
            if i < 0:
                i = self.shape[0] + i
            try:
                row = self._convert_row(i, j)
                return self.data[row, j]
            except OutsideBandException:
                return 0
        elif i == slice(None, None, None):
            start, stop = self.data_range(j)
            # TODO: return BandedMatrix
            return self.data[start:stop, j]
        else:
            raise Exception()

    def __setitem__(self, key, val):
        # TODO: implement multidim slicing and negative indices
        cdef int i, j, row
        i, j = key
        row = self._convert_row(i, j)
        self.data[row, j] = val

    def todense(self):
        result = np.zeros(self.shape)
        cdef int start, stop, dstart, dstop, j
        for j in range(self.shape[1]):
            start, stop = self.range(j)
            dstart, dstop = self.data_range(j)
            result[start:stop, j] = self.data[dstart:dstop, j]
        return result


cpdef update(BandedMatrix arr, int i, int j, int actual_j, cython.str s_base, cython.str t_base,
             float log_p, float log_ins, float log_del, int bandwidth):
    # TODO: log(1-p) may not be negligible
    cdef float sub = (0 if s_base == t_base else log_p)
    if i == actual_j - bandwidth:
        return max([arr[i, j - 1] + log_del,  # deletion
                    arr[i - 1, j - 1] + sub])
    elif i == actual_j + bandwidth:
        return max([arr[i - 1, j] + log_ins,  # insertion
                    arr[i - 1, j - 1] + sub])
    elif actual_j - bandwidth < i < actual_j + bandwidth:
        return max([arr[i - 1, j] + log_ins,  # insertion
                    arr[i, j - 1] + log_del,  # deletion
                    arr[i - 1, j - 1] + sub])
    else:
        raise Exception('update called outside bandwidth')


cpdef forward(cython.str s, float[:] log_p, cython.str t, float log_ins, float log_del, int bandwidth):
    result = BandedMatrix((len(s) + 1, len(t) + 1), bandwidth)
    cdef int i, j
    for i in range(min(bandwidth + 1, len(s) + 1)):
        result[i, 0] = log_ins * i
    for j in range(min(bandwidth + 1, len(t) + 1)):
        result[0, j] = log_del * j
    for i in range(1, len(s) + 1):
        for j in range(max(1, i - bandwidth), min(len(t) + 1, i + bandwidth + 1)):
            result[i, j] = update(result, i, j, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del, bandwidth)
    return result
