using Bio.Seq

using Quiver2.BandedArrays
using Quiver2.Sample
using Quiver2.Model
using Quiver2.Mutations

using Base.Test

function test_perfect_forward()
    log_ins = -5.0
    log_del = -10.0
    bandwidth = 1
    template = "AA"
    seq = "AA"
    log_p = fill(-3.0, length(seq))
    A = Model.forward(template, seq, log_p, log_ins, log_del, bandwidth)
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
    B = Model.backward(template, seq, log_p, log_ins, log_del, bandwidth)
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
    A = Model.forward(template, seq, log_p, log_ins, log_del, bandwidth)
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

function test_random_mutation(mutation, template_len)
    error_rate = 0.1
    sub_ratio = 2 / 10
    ins_ratio = 4 / 10
    del_ratio = 4 / 10
    error_std = 0.01
    log_ins = log10(error_rate * ins_ratio)
    log_del = log10(error_rate * del_ratio)
    template_seq = random_seq(template_len)
    template = convert(AbstractString, template_seq)
    bioseq, log_p = sample_from_template(template_seq, error_rate,
                                         sub_ratio, ins_ratio, del_ratio, error_std)
    seq = convert(AbstractString, bioseq)
    bandwidth = max(3 * abs(length(template) - length(seq)), 5)
    new_template = Mutations.update_template(template, mutation)
    M = Model.forward(new_template, seq, log_p, log_ins, log_del, bandwidth)

    i = mutation.pos + 2
    A = Model.forward(template, seq, log_p, log_ins, log_del, bandwidth)
    B = Model.backward(template, seq, log_p, log_ins, log_del, bandwidth)
    score = Model.score_mutation(mutation, template, seq, log_p,
                                 A, B, log_ins, log_del, bandwidth, false, 0.0, 0.0)
    @test_approx_eq score M[end, end]
    # TODO: test that inband values are equal in A and Acols.
end

function test_random_substitutions()
    for i = 1:30
        template_len = rand(10:20)
        pos = rand(1:template_len)
        mutation = Mutations.Substitution(pos, rbase())
        test_random_mutation(mutation, template_len)
    end
end

function test_random_insertions()
    for i = 1:30
        template_len = rand(10:20)
        pos = rand(0:template_len)
        mutation = Mutations.Insertion(pos, rbase())
        test_random_mutation(mutation, template_len)
    end
end

function test_random_codon_insertions()
    for i = 1:30
        template_len = rand(10:20)
        pos = rand(0:template_len)
        mutation = Mutations.CodonInsertion(pos, random_codon())
        test_random_mutation(mutation, template_len)
    end
end

function test_random_deletions()
    for i = 1:30
        template_len = rand(10:20)
        pos = rand(1:template_len)
        mutation = Mutations.Deletion(pos)
        test_random_mutation(mutation, template_len)
    end
end


function test_random_codon_deletions()
    for i = 1:30
        template_len = rand(10:20)
        pos = rand(1:(template_len - 2))
        mutation = Mutations.CodonDeletion(pos)
        test_random_mutation(mutation, template_len)
    end
end

function test_is_inframe()
    reference = "AAAGGGTTT"
    ref_log_p = ones(length(reference))
    ref_log_p = -2.0 * ones(length(reference))
    log_ins = -5.0
    log_del = -5.0
    codon_ins = -2.0
    codon_del = -2.0
    bandwidth = 6

    template = "AAACCCGGGTTT"
    @test Quiver2.Model.is_inframe(false, template, reference,
                                   ref_log_p, log_ins, log_del,
                                   bandwidth,
                                   true,
                                   codon_ins, codon_del)

    template = "AAACCCGGGTTTT"
    @test !Quiver2.Model.is_inframe(false, template, reference,
                                    ref_log_p, log_ins, log_del,
                                    bandwidth,
                                    true,
                                    codon_ins, codon_del)

    template = "AAA"
    @test Quiver2.Model.is_inframe(true, template, reference,
                                   ref_log_p, log_ins, log_del,
                                   bandwidth,
                                   true,
                                   codon_ins, codon_del)
end

function test_quiver2()
    # TODO: can't guarantee this test actually passes, since it is random
    n_seqs=10
    ref_len=90
    template_error_rate = 0.1
    t_sub_part = 8.0
    t_ins_part = 1.0
    t_del_part = 1.0

    mean_error_rate = 0.01
    max_error_rate = 0.01
    sub_ratio = 1.0 / 7.0
    ins_ratio = 3.0 / 7.0
    del_ratio = 3.0 / 7.0
    error_std = 0.01

    n = 100
    n_wrong = 0
    n_wrong_length = 0
    n_out_frame = 0

    for i in 1:n
        use_ref = rand([true, false])
        allow_stutters = rand([true, false])
        reference, template, template_log_p, reads, log_ps, error_rates =
            sample(n_seqs, ref_len,
                   template_error_rate,
                   t_sub_part, t_ins_part, t_del_part,
                   max_error_rate,
                   sub_ratio, ins_ratio, del_ratio,
                   error_std=error_std,
                   error_rate_alpha=3.0, error_rate_beta=1.0)
        if !use_ref
            reference = DNASequence("")
        end
        initial_template = reads[1]
        result, info = Model.quiver2(initial_template, reads,
                                     log_ps,
                                     log10(ins_ratio * mean_error_rate),
                                     log10(del_ratio * mean_error_rate);
                                     reference=reference,
                                     allow_stutters=allow_stutters,
                                     bandwidth=3, min_dist=9, batch=5,
                                     max_iters=100)
        if length(result) % 3 != 0
            n_out_frame += 1
        end
        if length(result) != length(template)
            n_wrong_length += 1
        end
        if result != template
            n_wrong += 1
        end
    end
    if n_wrong > 0
        println("wrong length : $(n_wrong_length) / $(n)")
        println("out of frame : $(n_out_frame) / $(n)")
        println("wrong        : $(n_wrong) / $(n)")
        @test false
    end
end

srand(12345)

test_perfect_forward()
test_perfect_backward()
test_imperfect_forward()
test_equal_ranges()
test_random_substitutions()
test_random_insertions()
test_random_codon_insertions()
test_random_deletions()
test_random_codon_deletions()
test_is_inframe()
test_quiver2()
