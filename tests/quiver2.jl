include("../BandedArray.jl")
include("../quiver2.jl")

using BandedArrayModule
using Quiver2

using Base.Test

# perfect_forward
begin
    log_del = -10.0
    log_ins = -5.0
    bandwidth = 1
    template = "AA"
    seq = "AA"
    log_p = fill(-3.0, length(seq))
    A = Quiver2.forward(seq, log_p, template, log_ins, log_del, bandwidth)
    expected = reshape([[0.0, -10.0, 0.0];
                        [-5.0, 0.0, -10.0];
                        [0.0,-5.0, 0.0]],
                       (3, 3))
    @test full(A) == expected
end

# # backward
# begin
#     log_del = -10
#     log_ins = -5
#     template = 'AA'
#     seq = 'AT'
#     log_p = np.repeat(-3, len(seq))
#     B = backward(seq, log_p, template, log_ins, log_del, bandwidth=1).full()
#     expected = np.array([[-3, -5, 0],
#                          [-13, -3, -5],
#                          [0, -10, 0]])
#     assert_array_equal(B, expected)
# end

# # imperfect_forward
# begin
#     log_del = -10
#     log_ins = -5
#     template = 'AA'
#     seq = 'AT'
#     log_p = np.repeat(-3, len(seq))
#     A = forward(seq, log_p, template, log_ins, log_del, bandwidth=1).full()
#     expected = np.array([[  0, -10, 0],
#                          [ -5,  0,  -10],
#                          [0, -5,  -3]])
#     assert_array_equal(A, expected)
# end

