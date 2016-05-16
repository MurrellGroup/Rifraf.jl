using Quiver2.BandedArrays
using Quiver2.Sample
using Quiver2.Model
using Quiver2.Mutations

using Base.Test

# tests
function test_perfect_forward()
    log_ins = -5.0
    log_del = -10.0
    bandwidth = 1
    template = "AA"
    seq = "AA"
    log_p = fill(-3.0, length(seq))
    A = Model.forward_seq(template, seq, log_p, log_ins, log_del, bandwidth)
    # transpose because of column-major order
    expected = transpose(reshape([[0.0, -10.0, 0.0];
                                  [-5.0, 0.0, -10.0];
                                  [0.0,-5.0, 0.0]],
                                 (3, 3)))
    @test full(A) == expected
end

function test_perfect_backward()
    log_del = -10.0
    log_ins = -5.0
    bandwidth = 1
    template = "AA"
    seq = "AT"
    log_p = fill(-3.0, length(seq))
    B = Model.backward_seq(template, seq, log_p, log_ins, log_del, bandwidth)
    expected = transpose(reshape([[-3, -5, 0];
                                  [-13, -3, -5];
                                  [0, -10, 0]],
                                 (3, 3)))
    @test full(B) == expected
end

function test_imperfect_forward()
    log_del = -10.0
    log_ins = -5.0
    bandwidth = 1
    template = "AA"
    seq = "AT"
    log_p = fill(-3.0, length(seq))
    A = Model.forward_seq(template, seq, log_p, log_ins, log_del, bandwidth)
    expected = transpose(reshape([[  0, -10, 0];
                                  [ -5,  0,  -10];
                                  [0, -5,  -3]],
                                 (3, 3)))
    @test full(A) == expected
end

function test_equal_ranges()
    @test Model.equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
    @test Model.equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
    @test Model.equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
end

function test_random_mutations()
    error_rate = 0.1
    sub_ratio = 2 / 10
    ins_ratio = 4 / 10
    del_ratio = 4 / 10
    log_ins = log10(error_rate * ins_ratio)
    log_del = log10(error_rate * del_ratio)
    for i = 1:100
        template_len = rand(10:20)
        template_seq = random_seq(template_len)
        template = convert(AbstractString, template_seq)
        bioseq, log_p = sample_from_template(template_seq, error_rate,
                                      sub_ratio, ins_ratio, del_ratio)
        seq = convert(AbstractString, bioseq)
        bandwidth = max(2 * abs(length(template) - length(seq)), 5)
        T = [Mutations.Substitution, Mutations.Insertion, Mutations.Deletion][rand(1:3)]
        maxpos = (T == Mutations.Insertion ? length(template) + 1: length(template))
        if T == Mutations.Deletion
            mutation = T(rand(1:maxpos))
        else
            mutation = T(rand(1:maxpos), rbase())
        end
        new_template = Model.update_template(template, mutation)
        A = Model.forward_seq(template, seq, log_p, log_ins, log_del, bandwidth)
        B = Model.backward_seq(template, seq, log_p, log_ins, log_del, bandwidth)
        M = Model.forward_seq(new_template, seq, log_p, log_ins, log_del, bandwidth)
        score = Model.score_mutation(mutation, template, seq, log_p,
                                     A, B, log_ins, log_del, bandwidth, false, 0.0, 0.0)
        @test_approx_eq score M[end, end]
    end
end

function test_apply_mutations()
    template = "ACG"
    mutations = [Mutations.Insertion(1, 'T'),
                 Mutations.Deletion(3),
                 Mutations.Substitution(2, 'T')]
    expected = "TAT"
    result = Model.apply_mutations(template, mutations)
    @test result == expected
end

function test_quiver2()
    # TODO: can't guarantee this test actually passes, since it is random
    mean_error_rate = 0.1
    max_error_rate = 0.2
    sub_ratio = 1.0 / 7.0
    ins_ratio = 3.0 / 7.0
    del_ratio = 3.0 / 7.0
    for i in 1:100
        reference, template, template_log_p, reads, log_ps = sample(30, 30, max_error_rate,
                                                                    sub_ratio, ins_ratio, del_ratio)
        initial_template = reads[1]
        result, info = Model.quiver2(reference, initial_template, reads,
                                     log_ps,
                                     log10(ins_ratio * mean_error_rate),
                                     log10(del_ratio * mean_error_rate),
                                     use_ref=true,
                                     bandwidth=3, min_dist=9, batch=20,
                                     verbose=false)
        @test result == template
    end
end

srand(1)

test_perfect_forward()
test_perfect_backward()
test_imperfect_forward()
test_equal_ranges()
test_random_mutations()
test_apply_mutations()
test_quiver2()
