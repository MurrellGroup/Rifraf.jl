using Bio.Seq

using Quiver2.BandedArrays
using Quiver2.Sample
using Quiver2.Model
using Quiver2.Mutations
using Quiver2.Util

using Base.Test


const errors = Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0)


function test_perfect_forward()
    bandwidth = 1
    template = "AA"
    seq = "AA"
    lp = -3.0
    match = Util.inv_log10(lp)
    log_p = fill(lp, length(seq))
    A = Model.forward(template, seq, log_p, errors, bandwidth)
    # transpose because of column-major order
    expected = transpose(reshape([[0.0, lp + errors.deletion, 0.0];
                                  [lp + errors.insertion, match, match + lp + errors.deletion];
                                  [0.0, match + lp + errors.insertion, 2 * match]],
                                 (3, 3)))
    @test full(A) == expected
end

function test_imperfect_backward()
    bandwidth = 1
    template = "AA"
    seq = "AT"
    lp = -3.0
    match = Util.inv_log10(lp)
    log_p = fill(lp, length(seq))
    B = Model.backward(template, seq, log_p, errors, bandwidth)
    expected = transpose(reshape([[lp + errors.mismatch + match, lp + errors.insertion + match, 0.0];
                                  [2*lp + errors.deletion + errors.mismatch, lp + errors.mismatch, lp + errors.insertion];
                                  [0.0, lp + errors.deletion, 0.0]],
                                 (3, 3)))
    @test full(B) == expected
end

function test_imperfect_forward()
    bandwidth = 1
    template = "AA"
    seq = "AT"
    lp = -3.0
    match = Util.inv_log10(lp)
    log_p = fill(lp, length(seq))
    A = Model.forward(template, seq, log_p, errors, bandwidth)
    B = Model.backward(template, seq, log_p, errors, bandwidth)
    expected = transpose(reshape([[0.0, lp + errors.deletion, 0.0];
                                  [lp + errors.insertion, match, match + lp + errors.deletion];
                                  [0.0, match + lp + errors.insertion, match + lp + errors.mismatch]],
                                 (3, 3)))

    @test full(A) == expected
    @test A[end, end] == B[1, 1]
end

function test_equal_ranges()
    @test Model.equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
    @test Model.equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
    @test Model.equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
end


function test_insertion_agreement()
    template = "AA"
    seq = "ATA"
    bandwidth = 10
    log_p = [-5.0, -1.0, -6.0]
    A = Model.forward(template, seq, log_p, errors, bandwidth)
    B = Model.backward(template, seq, log_p, errors, bandwidth)
    score = (Util.inv_log10(log_p[1]) +
             log_p[2] + errors.insertion +
             Util.inv_log10(log_p[3]))
    @test_approx_eq A[end, end] score
    @test_approx_eq B[1, 1] score
end

function test_deletion_agreement()
    template = "GATAG"
    seq = "GAAG"
    bandwidth = 10
    log_p = [-5.0, -2.0, -1.0, -6.0]
    A = Model.forward(template, seq, log_p, errors, bandwidth)
    B = Model.backward(template, seq, log_p, errors, bandwidth)
    score = (Util.inv_log10(log_p[1]) +
             Util.inv_log10(log_p[2]) +
             maximum(log_p[2:3]) + errors.deletion +
             Util.inv_log10(log_p[3]) +
             Util.inv_log10(log_p[4]))
    @test_approx_eq A[end, end] score
    @test_approx_eq B[1, 1] score
end

function test_deletion_agreement2()
    template = "ATA"
    seq = "AA"
    bandwidth = 10
    log_p = [-2.0, -3.0]
    A = Model.forward(template, seq, log_p, errors, bandwidth)
    B = Model.backward(template, seq, log_p, errors, bandwidth)
    score = (Util.inv_log10(log_p[1]) +
             maximum(log_p[1:2]) + errors.deletion +
             Util.inv_log10(log_p[2]))
    @test_approx_eq A[end, end] score
    @test_approx_eq B[1, 1] score
end

function test_random_mutation(mutation, template_len)
    template_seq = random_seq(template_len)
    template = convert(AbstractString, template_seq)

    template_error_p = Float64[0.1 for i=1:length(template)]
    codon_moves = rand([true, false])
    if codon_moves
        error_model = ErrorModel(2.0, 0.1, 0.1, 3.0, 3.0)
    else
        error_model = ErrorModel(2.0, 4.0, 4.0, 0.0, 0.0)
    end
    log_actual_error_std = 0.5
    log_reported_error_std = 0.5

    bioseq, actual, phreds = sample_from_template(template_seq,
                                                  template_error_p,
                                                  errors,
                                                  log_actual_error_std,
                                                  log_reported_error_std)
    log_p = Float64[Float64(q) / (-10.0) for q in phreds]
    seq = convert(AbstractString, bioseq)
    bandwidth = max(3 * abs(length(template) - length(seq)), 5)
    new_template = Mutations.update_template(template, mutation)
    Anew = Model.forward(new_template, seq, log_p, error_model, bandwidth)
    Bnew = Model.backward(new_template, seq, log_p, error_model, bandwidth)
    @test_approx_eq Anew[end, end] Bnew[1, 1]

    A = Model.forward(template, seq, log_p, errors, bandwidth)
    B = Model.backward(template, seq, log_p, errors, bandwidth)
    score = Model.seq_score_mutation(mutation, A, B, template, seq, log_p,
                                     error_model, codon_moves)
    @test_approx_eq score Anew[end, end]
    # TODO: test that inband values are equal in A and Acols.
end

function test_random_substitutions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(1:template_len)
        mutation = Mutations.Substitution(pos, rbase())
        test_random_mutation(mutation, template_len)
    end
end

function test_random_insertions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(0:template_len)
        mutation = Mutations.Insertion(pos, rbase())
        test_random_mutation(mutation, template_len)
    end
end

function test_random_codon_insertions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(0:template_len)
        mutation = Mutations.CodonInsertion(pos, random_codon())
        test_random_mutation(mutation, template_len)
    end
end

function test_random_deletions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(1:template_len)
        mutation = Mutations.Deletion(pos)
        test_random_mutation(mutation, template_len)
    end
end


function test_random_codon_deletions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(1:(template_len - 2))
        mutation = Mutations.CodonDeletion(pos)
        test_random_mutation(mutation, template_len)
    end
end


function test_no_single_indels()
    reference = "AAAGGGTTT"
    ref_log_p = fill(log10(0.01), length(reference))
    bandwidth = 6

    template = "AAACCCGGGTTT"
    @test Quiver2.Model.no_single_indels(template, reference,
                                         ref_log_p, errors,
                                         bandwidth)

    template = "AAACCCGGGTTTT"
    @test !Quiver2.Model.no_single_indels(template, reference,
                                          ref_log_p, errors,
                                          bandwidth)

    template = "AAA"
    @test Quiver2.Model.no_single_indels(template, reference,
                                         ref_log_p, errors,
                                         bandwidth)
end

function test_quiver2()
    # TODO: can't guarantee this test actually passes, since it is random
    n_seqs=10
    ref_len=90
    template_error_rate = 0.03
    template_errors = Model.ErrorModel(8.0, 0.0, 0.0, 1.0, 1.0)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)

    template_error_mean = 0.005
    template_error_std = 0.001
    log_seq_actual_std = 0.2
    log_seq_reported_std = 0.2
    seq_errors = Model.ErrorModel(1.0 / 7.0, 3.0 / 7.0, 3.0 / 7.0, 0.0, 0.0)

    n = 100
    n_wrong = 0
    n_wrong_length = 0
    n_out_frame = 0

    for i in 1:n
        use_ref = rand([true, false])
        (reference, template, template_error, reads,
         actual, phreds) = sample(n_seqs, ref_len,
                                  template_error_rate,
                                  template_errors,
                                  template_error_mean,
                                  template_error_std,
                                  log_seq_actual_std,
                                  log_seq_reported_std,
                                  seq_errors)
        if !use_ref
            reference = DNASequence("")
        end
        initial_template = reads[1]

        result, q1, q2, info = Model.quiver2(initial_template, reads,
                                             phreds, seq_errors,
                                             reference=reference,
                                             ref_log_p=log10(template_error_rate),
                                             ref_errors=ref_errors,
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


function test_base_probs()
    template = "CGTAC"
    seqs = ["CGAC",
            "CGAC",
            "CGAC"]
    lps = Vector{Float64}[[-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0]]
    bandwidth = 3
    state = Quiver2.Model.initial_state(template, seqs, lps, errors, bandwidth)
    base, ins = Quiver2.Model.estimate_probs(state, seqs, lps, errors,
                                             "", Float64[], errors)
    @test base[1, 2] > 0.9
    del_probs = base[:, end]
    @test del_probs[1] < 1e-9
    @test del_probs[3] > 0.9
end


function test_ins_probs()
    template = "CGAT"
    seqs = ["CGTAT",
            "CGTAT",
            "CGTAT"]
    bandwidth = 3
    lps = Vector{Float64}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    state = Quiver2.Model.initial_state(template, seqs, lps, errors, bandwidth)
    base, ins = Quiver2.Model.estimate_probs(state, seqs, lps, errors,
                                             "", Float64[], errors)
    @test maximum(ins[1, :]) < 1e-9
    @test ins[3, 4] > 0.9
end

function test_indel_probs()
    template = "CGAT"
    seqs = ["CGTAT",
            "CGTAT",
            "CGTAT"]
    bandwidth = 3
    lps = Vector{Float64}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    state = Quiver2.Model.initial_state(template, seqs, lps, errors, bandwidth)
    base, ins = Quiver2.Model.estimate_probs(state, seqs, lps, errors,
                                             "", Float64[], errors)
    probs = Quiver2.Model.estimate_indel_probs(base, ins)
    @test probs[1] == probs[4]
    @test probs[1] < 0.5
    @test probs[2] == probs[3]
    @test probs[2] > 0.5
end

function test_align()
    template = "ATA"
    seq = "AA"
    bandwidth = 10
    log_p = [-2.0, -3.0]
    t, s = Model.align(template, seq, log_p, errors, bandwidth)
    @test t == "ATA"
    @test s == "A-A"
end

srand(1234)

test_perfect_forward()
test_imperfect_backward()
test_imperfect_forward()
test_equal_ranges()
test_insertion_agreement()
test_deletion_agreement()
test_deletion_agreement2()
test_random_substitutions()
test_random_insertions()
test_random_codon_insertions()
test_random_deletions()
test_random_codon_deletions()
test_no_single_indels()
test_quiver2()
test_base_probs()
test_ins_probs()
test_indel_probs()
test_align()
