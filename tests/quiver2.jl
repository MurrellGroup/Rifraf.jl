include("../BandedArray.jl")
include("../quiver2.jl")

using BandedArrayModule
using Quiver2

using Base.Test

# perfect_forward
begin
    log_ins = -5.0
    log_del = -10.0
    bandwidth = 1
    template = "AA"
    seq = "AA"
    log_p = fill(-3.0, length(seq))
    A = Quiver2.forward(seq, log_p, template, log_ins, log_del, bandwidth)
    # transpose because of column-major order
    expected = transpose(reshape([[0.0, -10.0, 0.0];
                                  [-5.0, 0.0, -10.0];
                                  [0.0,-5.0, 0.0]],
                                 (3, 3)))
    @test full(A) == expected
end

# backward
begin
    log_del = -10.0
    log_ins = -5.0
    bandwidth = 1
    template = "AA"
    seq = "AT"
    log_p = fill(-3.0, length(seq))
    B = Quiver2.backward(seq, log_p, template, log_ins, log_del, bandwidth)
    expected = transpose(reshape([[-3, -5, 0];
                                  [-13, -3, -5];
                                  [0, -10, 0]],
                                 (3, 3)))
    @test full(B) == expected
end

# imperfect_forward
begin
    log_del = -10.0
    log_ins = -5.0
    bandwidth = 1
    template = "AA"
    seq = "AT"
    log_p = fill(-3.0, length(seq))
    A = Quiver2.forward(seq, log_p, template, log_ins, log_del, bandwidth)
    expected = transpose(reshape([[  0, -10, 0];
                                  [ -5,  0,  -10];
                                  [0, -5,  -3]],
                                 (3, 3)))
    @test full(A) == expected
end

# test updated_col substitution
begin
    log_del = -10.0
    log_ins = -5.0
    bandwidth = 1
    template = "TAAG"
    seq = "ATAG"
    pos = 2
    log_p = fill(-3.0, length(seq))
    A = Quiver2.forward(seq, log_p, template, log_ins, log_del, bandwidth)
    result = Quiver2.updated_col(pos, pos + 1, 'T', template, seq, log_p, A, log_ins, log_del)
    Ae = Quiver2.forward(seq, log_p, "TTAG", log_ins, log_del, bandwidth)
    expected = sparsecol(Ae, pos + 1)
    @test result == expected
end

function rbase()
    bases = ['A', 'C', 'G', 'T']
    return bases[rand(1:4)]
end

function mutate_base(base::Char)
    result = rbase()
    while result == base
        result = rbase()
    end
    return result
end

function random_seq(n)
    return join([rbase() for i in 1:n])
end

function coinflip(p)
    return rand() < p
end

function sample_from_template(template, point_rate, insertion_rate, deletion_rate)
    result = []
    for base in template
        while coinflip(insertion_rate)
            push!(result, rbase())
        end
        if coinflip(deletion_rate)
            continue
        end
        if coinflip(point_rate)
            push!(result, mutate_base(base))
        else
            push!(result, base)
        end
    end
    while coinflip(insertion_rate)
        push!(result, rbase())
    end
    return string(result)
end

function phred(p)
    return -10 * log10(p)
end

# test_random_mutations
begin
    point_rate = 0.1
    insertion_rate = 0.01
    deletion_rate = 0.01
    log_ins = log10(insertion_rate)
    log_del = log10(deletion_rate)
    for i = 1:1000
        template_len = rand(10:20)
        template_seq = random_seq(template_len)
        template = template_seq
        seq = sample_from_template(template_seq, point_rate, insertion_rate, deletion_rate)
        bandwidth = max(2 * abs(length(template) - length(seq)), 5)
        m = ["substitution", "insertion", "deletion"][rand(1:3)]
        maxpos = (m == "insertion" ? length(template) + 1: length(template))
        mutation = Quiver2.Mutation(m, rand(1:maxpos), rbase())
        new_template = Quiver2.update_template(template, mutation)
        phreds = phred(fill(point_rate + insertion_rate + deletion_rate, length(seq)))
        log_p = -phreds / 10
        A = Quiver2.forward(seq, log_p, template, log_ins, log_del, bandwidth)
        B = Quiver2.backward(seq, log_p, template, log_ins, log_del, bandwidth)
        M = Quiver2.forward(seq, log_p, new_template, log_ins, log_del, bandwidth)
        if mutation.m == "substitution"
            # this is the only case where the updated column exactly matches the full result
            col = Quiver2.updated_col(pos, pos + 1, mutation.base, template, seq, log_p, A, log_ins, log_del)
            exp_col = sparsecol(M, pos + 1)
            @test col == exp_col
        end
        score = Quiver2.score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth)
        @test_approx_eq score M[end, end]
        @test_approx_eq score Quiver2.backward(seq, log_p, new_template, log_ins, log_del, bandwidth)[1, 1]
    end
end
