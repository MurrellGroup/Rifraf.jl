using Bio.Seq

using Quiver2.BandedArrays
using Quiver2.Sample
using Quiver2.Model
using Quiver2.Proposals
using Quiver2.Util

using Base.Test


const errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
const scores = Model.Scores(errors)


function inv_log10(logp::Float64)
    log10(1.0 - exp10(logp))
end


function test_cols(A, B, codon_moves)
    expected = A[end, end]
    @test_approx_eq A[end, end] B[1, 1]
    ncols = size(A)[2]
    # if codon_moves is true, we cannot expect every column to contain
    # the correct score
    # TODO: every three columns should
    if !codon_moves
        for j in 1:ncols
            Acol = sparsecol(A, j)
            Bcol = sparsecol(B, j)
            score = maximum(Acol + Bcol)
            @test_approx_eq expected score
        end
    end
    @test_approx_eq A[end, end] B[1, 1]
end

function test_perfect_forward()
    bandwidth = 1
    template = "AA"
    seq = "AA"
    lp = -3.0
    match = inv_log10(lp)
    log_p = fill(lp, length(seq))
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    A = Model.forward(template, pseq, scores)
    # transpose because of column-major order
    expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                  [lp + scores.insertion, match, match + lp + scores.deletion];
                                  [0.0, match + lp + scores.insertion, 2 * match]],
                                 (3, 3)))
    @test full(A) == expected
end

function test_imperfect_backward()
    bandwidth = 1
    template = "AA"
    seq = "AT"
    lp = -3.0
    match = inv_log10(lp)
    log_p = fill(lp, length(seq))
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    B = Model.backward(template, pseq, scores)
    expected = transpose(reshape([[lp + scores.mismatch + match, lp + scores.insertion + match, 0.0];
                                  [2*lp + scores.deletion + scores.mismatch, lp + scores.mismatch, lp + scores.insertion];
                                  [0.0, lp + scores.deletion, 0.0]],
                                 (3, 3)))
    @test full(B) == expected
end

function test_imperfect_forward()
    bandwidth = 1
    template = "AA"
    seq = "AT"
    lp = -3.0
    match = inv_log10(lp)
    log_p = fill(lp, length(seq))
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    A = Model.forward(template, pseq, scores)
    B = Model.backward(template, pseq, scores)
    test_cols(A, B, false)
    expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                  [lp + scores.insertion, match, match + lp + scores.deletion];
                                  [0.0, match + lp + scores.insertion, match + lp + scores.mismatch]],
                                 (3, 3)))

    @test full(A) == expected
end

function test_equal_ranges()
    @test Model.equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
    @test Model.equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
    @test Model.equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
end

function test_forward_backward_agreement()
    template = "TG"
    seq = "GTCG"
    log_p = [-1.2, -0.8, -0.7, -1.0]
    bandwidth = 5
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    local_scores = Model.Scores(Model.ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
    A = Model.forward(template, pseq, local_scores)
    B = Model.backward(template, pseq, local_scores)
    test_cols(A, B, true)
end

function test_forward_backward_agreement_2()
    template = "GCACGGTC"
    seq = "GACAC"
    log_p = [-1.1, -1.1, -0.4, -1.0, -0.7]
    bandwidth = 5
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    local_scores = Model.Scores(Model.ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
    A = Model.forward(template, pseq, local_scores)
    B = Model.backward(template, pseq, local_scores)
    test_cols(A, B, true)
end

function test_insertion_agreement()
    template = "AA"
    seq = "ATA"
    bandwidth = 10
    log_p = [-5.0, -1.0, -6.0]
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    A = Model.forward(template, pseq, scores)
    B = Model.backward(template, pseq, scores)
    score = (inv_log10(log_p[1]) +
             log_p[2] + scores.insertion +
             inv_log10(log_p[3]))
    @test_approx_eq A[end, end] score
    test_cols(A, B, false)
end

function test_deletion_agreement()
    template = "GATAG"
    seq = "GAAG"
    bandwidth = 10
    log_p = [-5.0, -2.0, -1.0, -6.0]
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    A = Model.forward(template, pseq, scores)
    B = Model.backward(template, pseq, scores)
    score = (inv_log10(log_p[1]) +
             inv_log10(log_p[2]) +
             maximum(log_p[2:3]) + scores.deletion +
             inv_log10(log_p[3]) +
             inv_log10(log_p[4]))
    @test_approx_eq A[end, end] score
    test_cols(A, B, false)
end

function test_deletion_agreement2()
    template = "ATA"
    seq = "AA"
    bandwidth = 10
    log_p = [-2.0, -3.0]
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    A = Model.forward(template, pseq, scores)
    B = Model.backward(template, pseq, scores)
    score = (inv_log10(log_p[1]) +
             maximum(log_p[1:2]) + scores.deletion +
             inv_log10(log_p[2]))
    @test_approx_eq A[end, end] score
    test_cols(A, B, false)
end

function test_random_proposal(proposal, template_len)
    template_seq = random_seq(template_len)
    template = convert(String, template_seq)

    template_error_p = Float64[0.1 for i=1:length(template)]
    codon_moves = rand([true, false])
    if codon_moves
        local_errors = Model.ErrorModel(2.0, 0.1, 0.1, 3.0, 3.0)
    else
        local_errors = Model.ErrorModel(2.0, 4.0, 4.0, 0.0, 0.0)
    end
    local_scores = Model.Scores(local_errors)
    phred_scale = 3.0
    log_actual_error_std = 0.5
    log_reported_error_std = 0.5

    (bioseq, actual, phreds, sbools,
     tbools) = sample_from_template(template_seq,
                                    template_error_p,
                                    errors,
                                    phred_scale,
                                    log_actual_error_std,
                                    log_reported_error_std)
    log_p = Float64[Float64(q) / (-10.0) for q in phreds]
    seq = convert(String, bioseq)
    bandwidth = max(5 * abs(length(template) - length(seq)), 30)
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)

    new_template = Proposals.update_template(template, proposal)
    Anew = Model.forward(new_template, pseq, local_scores)
    Bnew = Model.backward(new_template, pseq, local_scores)
    test_cols(Anew, Bnew, codon_moves)

    A = Model.forward(template, pseq, local_scores)
    B = Model.backward(template, pseq, local_scores)
    test_cols(A, B, codon_moves)
    newcols = zeros(Float64, size(A)[1], 6)
    score = Model.seq_score_proposal(proposal, A, B, template, pseq,
                                     local_scores, newcols)
    @test_approx_eq_eps score Anew[end, end] 0.1
end

function test_random_substitutions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(1:template_len)
        proposal = Proposals.Substitution(pos, rbase())
        test_random_proposal(proposal, template_len)
    end
end

function test_random_insertions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(0:template_len)
        proposal = Proposals.Insertion(pos, rbase())
        test_random_proposal(proposal, template_len)
    end
end

function test_random_deletions()
    for i = 1:1000
        template_len = rand(30:50)
        pos = rand(1:template_len)
        proposal = Proposals.Deletion(pos)
        test_random_proposal(proposal, template_len)
    end
end


function test_no_single_indels()
    reference = "AAAGGGTTT"
    ref_log_p = fill(log10(0.01), length(reference))
    bandwidth = 6
    rseq = Quiver2.Model.PString(reference, ref_log_p, bandwidth)
    local_errors = Model.normalize(Model.ErrorModel(2.0, 0.5, 0.5, 1.0, 1.0))
    local_scores = Model.Scores(local_errors)


    template = "AAACCCGGGTTT"
    @test !Quiver2.Model.has_single_indels(template, rseq, local_scores)

    template = "AAACCCGGGTTTT"
    @test Quiver2.Model.has_single_indels(template, rseq,
                                          local_scores)

    template = "AAA"
    @test !Quiver2.Model.has_single_indels(template, rseq,
                                           local_scores)
end

function _test_surgery_proposals(template, seqs, lps, expected)
    bandwidth = 5
    rseq = Quiver2.Model.PString("", Float64[], bandwidth)
    mult = 2
    errors = Model.ErrorModel(1.0, 4.0, 5.0, 0.0, 0.0)
    scores = Model.Scores(errors)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    do_subs = true
    do_indels = true

    if length(lps) != length(seqs)
        lps =  Vector{Float64}[fill(-9.0, length(s)) for s in seqs]
    end

    pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                  for (s, p) in zip(seqs, lps)]
    state = Quiver2.Model.initial_state(template, pseqs)
    Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    proposals, deltas = Quiver2.Model.surgery_proposals(state, pseqs, scores, do_subs, do_indels)
    @test length(symdiff(Set(proposals), Set(expected))) == 0
end


function test_surgery_proposals()
    template = "ACGAG"
    seqs = ["CGTAC",
            "CGAC",
            "CGTAG"]
    expected = [Proposals.Deletion(1),
                Proposals.Insertion(3, 'T'),
                Proposals.Substitution(5, 'C')]
    _test_surgery_proposals(template, seqs, [], expected)

    template = "AA"
    seqs = ["AAG",
            "AA",
            "AAG"]
    expected = [Proposals.Insertion(2, 'G')]
    _test_surgery_proposals(template, seqs, [], expected)

    template = "AA"
    seqs = ["GAA",
            "AA",
            "GAA"]
    expected = [Proposals.Insertion(0, 'G')]
    _test_surgery_proposals(template, seqs, [], expected)

    # test that highly confident base overwhelms others
    template = "AA"
    seqs = ["GAA",
            "AA",
            "AA"]
    lps = Vector{Float64}[[-10.0, -5.0, -5.0],
                          [-3.0, -5.0],
                          [-3.0, -50]]
    expected = [Proposals.Insertion(0, 'G')]
    _test_surgery_proposals(template, seqs, lps, expected)

end

function test_quiver2()
    # can't guarantee this test actually passes, since it is random
    n_seqs = 3
    len = 30
    ref_error_rate = 0.3
    ref_sample_errors = Model.ErrorModel(8.0, 0.0, 0.0, 1.0, 1.0)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    template_error_mean = 0.0001
    template_alpha = 1.0
    phred_scale = 3.0
    log_seq_actual_std = 3.0
    log_seq_reported_std = 0.3

    seq_errors = Model.ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0)
    seq_scores = Model.Scores(seq_errors)

    n = 100
    n_wrong = 0
    n_wrong_length = 0
    n_out_frame = 0

    for i in 1:n
        use_ref = rand([true, false])
        do_alignment_proposals = rand([true, false])
        do_surgery_proposals = rand([true, false])
        trust_proposals = rand([true, false])
        fix_indels_stat = rand([true, false])
        indel_correction_only = rand([true, false])

        (reference, template, template_error, reads, actual, phreds, sbools,
         tbools) = sample(n_seqs, len,
                          ref_error_rate,
                          ref_sample_errors,
                          template_error_mean,
                          template_alpha,
                          phred_scale,
                          log_seq_actual_std,
                          log_seq_reported_std,
                          seq_errors)
        if !use_ref
            reference = DNASequence("")
        end

        (result, info) = Model.quiver2(reads,
                                       phreds, seq_scores;
                                       reference=reference,
                                       ref_scores=ref_scores,
                                       do_alignment_proposals=do_alignment_proposals,
                                       do_surgery_proposals=do_surgery_proposals,
                                       trust_proposals=trust_proposals,
                                       fix_indels_stat=fix_indels_stat,
                                       indel_correction_only=indel_correction_only,
                                       bandwidth=10, min_dist=9, batch=5,
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
    bandwidth = 5
    mult = 2
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                  for (s, p) in zip(seqs, lps)]
    rseq = Quiver2.Model.PString("", Float64[], bandwidth)
    state = Quiver2.Model.initial_state(template, pseqs)
    Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    probs = Quiver2.Model.estimate_probs(state, pseqs, scores,
                                         rseq, scores, false)
    @test probs.sub[1, 2] > 0.9
    @test probs.del[1] < 1e-9
    @test probs.del[3] > 0.9
end


function test_ins_probs()
    template = "CGAT"
    seqs = ["CGTAT",
            "CGTAT",
            "CGTAT"]
    bandwidth = 5
    mult = 2
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    lps = Vector{Float64}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                  for (s, p) in zip(seqs, lps)]
    rseq = Quiver2.Model.PString("", Float64[], bandwidth)
    state = Quiver2.Model.initial_state(template, pseqs)
    Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    probs = Quiver2.Model.estimate_probs(state, pseqs, scores,
                                         rseq, scores, false)
    @test maximum(probs.ins[1, :]) < 1e-9
    @test probs.ins[3, 4] > 0.9
end

function test_align()
    template = "ATAA"
    seq = "AAA"
    bandwidth = 10
    log_p = [-2.0, -3.0, -3.0]
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    moves = Model.align_moves(template, pseq, scores)
    t, s = Model.moves_to_alignment_strings(moves, template, seq)
    @test t == "ATAA"
    @test s == "A-AA"
end

function test_align_2()
    template = "AACCTT"
    seq = "AAACCCTT"
    bandwidth = 10
    log_p = fill(log10(0.1), length(seq))
    pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
    moves = Model.align_moves(template, pseq, scores)
    t, s = Model.moves_to_alignment_strings(moves, template, seq)
    @test t[end-1:end] == "TT"
end


srand(1234)

test_perfect_forward()
test_imperfect_backward()
test_imperfect_forward()
test_equal_ranges()
test_forward_backward_agreement()
test_forward_backward_agreement_2()
test_insertion_agreement()
test_deletion_agreement()
test_deletion_agreement2()
test_random_substitutions()
test_random_insertions()
test_random_deletions()
test_no_single_indels()
test_surgery_proposals()
test_quiver2()
test_base_probs()
test_ins_probs()
test_align()
test_align_2()
