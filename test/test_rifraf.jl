using Bio.Seq
using Base.Test

using Rifraf

import Rifraf.sample_from_template,
       Rifraf.random_seq,
       Rifraf.CandProposal,
       Rifraf.Proposal,
       Rifraf.Substitution,
       Rifraf.Insertion,
       Rifraf.Deletion,
       Rifraf.apply_proposal,
       Rifraf.rbase,
       Rifraf.forward,
       Rifraf.backward,
       Rifraf.phred_to_log_p,
       Rifraf.equal_ranges,
       Rifraf.initial_state,
       Rifraf.estimate_probs,
       Rifraf.has_single_indels,
       Rifraf.score_proposal,
       Rifraf.single_indel_proposals,
       Rifraf.recompute!,
       Rifraf.get_candidate_proposals,
       Rifraf.alignment_proposals,
       Rifraf.surgery_proposals,
       Rifraf.align_moves,
       Rifraf.moves_to_aligned_seqs,
       Rifraf.alignment_error_probs,
       Rifraf.Prob,
       Rifraf.LogProb,
       Rifraf.LogProb

srand(1234)


function inv_log10(logvals)
    log10(1.0 - exp10(logvals))
end

function check_all_cols(A, B, codon_moves)
    expected = A[end, end]
    if !(A[end, end] ≈ B[1, 1])
        return false
    end
    ncols = size(A)[2]
    # if codon_moves is true, we cannot expect every column to contain
    # the correct score
    # TODO: every three columns should
    if !codon_moves
        for j in 1:ncols
            Acol = sparsecol(A, j)
            Bcol = sparsecol(B, j)
            score = maximum(Acol + Bcol)
            if !(expected ≈ score)
                return false
            end
        end
    end
    return A[end, end] ≈ B[1, 1]
end


@testset "equal_ranges" begin
    @test equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
    @test equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
    @test equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
end


@testset "test forward and backward" begin
    const errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Scores(errors)

    @testset "perfect_forward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AA")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth)
        A = forward(template, pseq, scores)
        # transpose because of column-major order
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, 2 * match]],
                                     (3, 3)))
        @test full(A) == expected
    end

    @testset "imperfect_backward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AT")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth)
        B = backward(template, pseq, scores)
        expected = transpose(reshape([[lp + scores.mismatch + match,
                                       lp + scores.insertion + match, 0.0];
                                      [2*lp + scores.deletion + scores.mismatch,
                                       lp + scores.mismatch,
                                       lp + scores.insertion];
                                      [0.0, lp + scores.deletion, 0.0]],
                                     (3, 3)))
        @test full(B) == expected
    end

    @testset "imperfect_forward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AT")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth)
        A = forward(template, pseq, scores)
        B = backward(template, pseq, scores)
        check_all_cols(A, B, false)
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, match + lp + scores.mismatch]],
                                     (3, 3)))

        @test full(A) == expected
    end

    @testset "forward/backward agreement" begin
        @testset "forward/backward agreement 1" begin
            template = DNASeq("TG")
            seq = DNASeq("GTCG")
            log_p = [-1.2, -0.8, -0.7, -1.0]
            bandwidth = 5
            pseq = RifrafSequence(seq, log_p, bandwidth)
            local_scores = Scores(ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            A = forward(template, pseq, local_scores)
            B = backward(template, pseq, local_scores)
            check_all_cols(A, B, true)
        end

        @testset "forward/backward agreement" begin
            template = DNASeq("GCACGGTC")
            seq = DNASeq("GACAC")
            log_p = [-1.1, -1.1, -0.4, -1.0, -0.7]
            bandwidth = 5
            pseq = RifrafSequence(seq, log_p, bandwidth)
            local_scores = Scores(ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            A = forward(template, pseq, local_scores)
            B = backward(template, pseq, local_scores)
            check_all_cols(A, B, true)
        end
    end

    @testset "insertion_agreement" begin
        template = DNASeq("AA")
        seq = DNASeq("ATA")
        bandwidth = 10
        log_p = [-5.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth)
        A = forward(template, pseq, scores)
        B = backward(template, pseq, scores)
        score = (inv_log10(log_p[1]) +
                 log_p[2] + scores.insertion +
                 inv_log10(log_p[3]))
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)
    end

    @testset "deletion agreement" begin
        template = DNASeq("GATAG")
        seq = DNASeq("GAAG")
        bandwidth = 10
        log_p = [-5.0, -2.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth)
        A = forward(template, pseq, scores)
        B = backward(template, pseq, scores)
        score = (inv_log10(log_p[1]) +
                 inv_log10(log_p[2]) +
                 maximum(log_p[2:3]) + scores.deletion +
                 inv_log10(log_p[3]) +
                 inv_log10(log_p[4]))
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)
    end

    @testset "deletion agreement 2" begin
        template = DNASeq("ATA")
        seq = DNASeq("AA")
        bandwidth = 10
        log_p = [-2.0, -3.0]
        pseq = RifrafSequence(seq, log_p, bandwidth)
        A = forward(template, pseq, scores)
        B = backward(template, pseq, scores)
        score = (inv_log10(log_p[1]) +
                 maximum(log_p[1:2]) + scores.deletion +
                 inv_log10(log_p[2]))
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)
    end
end


@testset "scoring proposals" begin
    const errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Scores(errors)

    function test_random_proposal(proposal, template_len)
        template = random_seq(template_len)

        template_error_p = Prob[0.1 for i=1:length(template)]
        codon_moves = rand([true, false])
        if codon_moves
            local_errors = ErrorModel(2.0, 0.1, 0.1, 3.0, 3.0)
        else
            local_errors = ErrorModel(2.0, 4.0, 4.0, 0.0, 0.0)
        end
        local_scores = Scores(local_errors)
        phred_scale = 3.0
        log_actual_error_std = 0.5
        log_reported_error_std = 0.5

        (seq, actual, phreds, sbools,
         tbools) = sample_from_template(template,
                                        template_error_p,
                                        errors,
                                        phred_scale,
                                        log_actual_error_std,
                                        log_reported_error_std)
        log_p = LogProb[float(q) / (-10.0) for q in phreds]
        bandwidth = max(5 * abs(length(template) - length(seq)), 30)
        pseq = RifrafSequence(seq, log_p, bandwidth)

        new_template = apply_proposal(template, proposal)
        Anew = forward(new_template, pseq, local_scores)
        Bnew = backward(new_template, pseq, local_scores)
        check_all_cols(Anew, Bnew, codon_moves)

        A = forward(template, pseq, local_scores)
        B = backward(template, pseq, local_scores)
        check_all_cols(A, B, codon_moves)
        newcols = zeros(size(A)[1], 6)
        score = score_proposal(proposal, A, B, template, pseq,
                               local_scores, newcols)
        @test_approx_eq_eps score Anew[end, end] 0.1
    end

    @testset "random_substitutions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(1:template_len)
            proposal = Substitution(pos, rbase())
            test_random_proposal(proposal, template_len)
        end
    end

    @testset "random_insertions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(0:template_len)
            proposal = Insertion(pos, rbase())
            test_random_proposal(proposal, template_len)
        end
    end

    @testset "random_deletions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(1:template_len)
            proposal = Deletion(pos)
            test_random_proposal(proposal, template_len)
        end
    end
end

@testset "no single indels" begin
    reference = DNASeq("AAAGGGTTT")
    ref_log_p = fill(log10(0.01), length(reference))
    bandwidth = 6
    rseq = RifrafSequence(reference, ref_log_p, bandwidth)
    local_errors = Rifraf.normalize(ErrorModel(2.0, 0.5, 0.5, 1.0, 1.0))
    local_scores = Scores(local_errors)

    template = DNASeq("AAACCCGGGTTT")
    @test !has_single_indels(template, rseq, local_scores)

    template = DNASeq("AAACCCGGGTTTT")
    @test has_single_indels(template, rseq, local_scores)

    template = DNASeq("AAA")
    @test !has_single_indels(template, rseq, local_scores)
end


@testset "single indel proposals" begin
    ref_errors = ErrorModel(10.0, 1e-10, 1e-10, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    ref = DNASeq("CGGCGATTT")
    consensus_errors = LogProb[-8.04822,-5.10032,-5.09486,-1.0,-2.68901,-6.52537,-5.20094]
    consensus = RifrafSequence(DNASeq("CTGCCGA"), consensus_errors, 10)

    proposals = single_indel_proposals(ref, consensus, ref_scores)
    expected = [Deletion(4)]
    @test expected == proposals
end


@testset "fast proposals" begin
    function _test_fast_proposals(consensus, seqs, lps, expected;
                                  do_alignment::Bool=true,
                                  do_surgery::Bool=true,
                                  do_subs::Bool=true,
                                  do_indels::Bool=true)
        bandwidth = 6
        padding = 3
        rseq = RifrafSequence(DNASeq(), LogProb[], bandwidth)
        mult = 2
        errors = ErrorModel(1.0, 5.0, 5.0, 0.0, 0.0)
        scores = Scores(errors)
        ref_errors = ErrorModel(10.0, 1e-10, 1e-10, 1.0, 1.0)
        ref_scores = Scores(ref_errors)

        if length(lps) == 0
            lps =  Vector{LogProb}[fill(-9.0, length(s)) for s in seqs]
        end

        pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                               for (s, p) in zip(seqs, lps)]
        state = initial_state(consensus, pseqs, DNASeq(""), bandwidth, padding)
        recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)

        if do_alignment
            proposals = alignment_proposals(state.Amoves, state.consensus,
                                            pseqs, do_subs, do_indels)
            @test length(symdiff(Set(proposals), Set(expected))) == 0
        end
        if do_surgery
            proposals, deltas = surgery_proposals(state, pseqs, scores, do_subs, do_indels)
            @test length(symdiff(Set(proposals), Set(expected))) == 0
        end
    end

    @testset "fast proposals 1" begin
        consensus = DNASeq("ACGAG")
        seqs = [DNASeq("CGTAC"),
                DNASeq("CGAC"),
                DNASeq("CGTAG")]
        expected = [Deletion(1),
                    Insertion(3, DNA_T),
                    Substitution(5, DNA_C)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals 2" begin
        consensus = DNASeq("AA")
        seqs = [DNASeq("AAG"),
                DNASeq("AA"),
                DNASeq("AAG")]
        expected = [Insertion(2, DNA_G)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals 3" begin
        consensus = DNASeq("AA")
        seqs = [DNASeq("GAA"),
                DNASeq("AA"),
                DNASeq("GAA")]
        expected = [Insertion(0, DNA_G)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals 4" begin
        consensus = DNASeq("AA")
        seqs = [DNASeq("AGA"),
                DNASeq("AA"),
                DNASeq("AGA")]
        expected = [Insertion(1, DNA_G)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals - highly confident base" begin
        consensus = DNASeq("AA")
        seqs = [DNASeq("GAA"),
                DNASeq("AA"),
                DNASeq("AA")]
        lps = Vector{LogProb}[[-10.0, -5.0, -5.0],
                                   [-3.0, -5.0],
                                   [-3.0, -50]]
        expected = [Insertion(0, DNA_G)]
        _test_fast_proposals(consensus, seqs, lps, expected)
    end
    @testset "push deletions" begin
        # test that deletions get pushed
        consensus = DNASeq("CCGTAAAC")
        seqs = [DNASeq("CGTAAAC"), DNASeq("CCGTAAAC"), DNASeq("CGTAAAC")]
        lps = Vector{LogProb}[[-1.0,-0.6,-1.1,-0.5,-0.5,-1.2,-2.0],
                                   [-1.0,-1.0,-1.6,-2.2,-0.8,-0.5,-1.2,-1.8],
                                   [-1.0,-1.8,-1.0,-0.5,-0.6,-1.0,-1.4]]
        expected = [Deletion(2)]
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true)
    end
    @testset "push deletions - simple" begin
        # test that deletions get pushed to end
        consensus = DNASeq("CCC")
        seqs = [DNASeq("CC")]
        expected = [Deletion(3)]
        _test_fast_proposals(consensus, seqs, [], expected,
                             do_alignment=false, do_surgery=true)
    end

    @testset "push insertions" begin
        # test that insertions get pushed to end
        consensus = DNASeq("TT")
        seqs = [DNASeq("TTT"),
                DNASeq("CTT")]
        expected = [Insertion(0, DNA_C),
                    Insertion(2, DNA_T)]
        _test_fast_proposals(consensus, seqs, [], expected,
                             do_alignment=false, do_surgery=true)
    end

    @testset "push subs" begin
        # test that subs get pushed backwards
        consensus = DNASeq("ATG")
        seqs = [DNASeq("ATG"), DNASeq("CATG"), DNASeq("CTG")]
        lps = [[-0.8, -1.7, -0.8],
               [-0.7, -0.5, -1.2, -1.0],
               [-0.9, -1.9, -1.0]]
        expected = [Substitution(1, DNA_C)]
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true,
                             do_subs=true, do_indels=false)
    end

    @testset "test fast proposals converged" begin
        consensus = DNASeq("GTTCGGCTC")
        seqs = [DNASeq("GTTCGGCTTC"),
                DNASeq("GTTCGGCTC"),
                DNASeq("GTTCCTG")]
        phreds = Vector{Phred}[[28,16,13,21,15,13,13,12,20,16],
                               [21,16,9,17,6,15,6,16,12],
                               [26,14,5,24,8,12,7]]
        lps = map(phred_to_log_p, phreds)
        expected = []
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true)

    end

@testset "test fast proposals converged 2" begin
    consensus = DNASeq("CAGTGCCGG")
    seqs = [DNASeq("CATGCCGG"),
            DNASeq("CATGCCCTGG"),
            DNASeq("CAGGGCCGG")]
    phreds = Vector{Phred}[[13,7,12,13,7,11,6,14],
                           [16,14,14,20,5,5,15,12,10,20],
                           [23,9,7,6,9,10,10,10,23]]
    lps = map(phred_to_log_p, phreds)
    expected = []
    # no indels, because this happened during refinement
    _test_fast_proposals(consensus, seqs, lps, expected,
                         do_alignment=false, do_surgery=true,
                         do_subs=true, do_indels=false)
end
@testset "test fast proposals insertion/mismatch swap" begin
    # inserting 'T' after position 5 should swap with G/T
    # mismatch in second sequence
    consensus = DNASeq("GGAAGTCC")
    seqs = [DNASeq("GGAAGTCC"),
            DNASeq("GGAATTCC"),
            DNASeq("GGAAGTCTACC")]
    phreds = Vector{Phred}[[18,9,12,14,11,11,15,14],
                           [25,10,6,8,12,11,19,13],
                           [24,12,9,15,8,8,8,8,8,23,19]]
    lps = map(phred_to_log_p, phreds)
    expected = [Substitution(5, DNA_T),
                Insertion(6, DNA_T)]
    _test_fast_proposals(consensus, seqs, lps, expected,
                         do_alignment=false, do_surgery=true,
                         do_subs=true, do_indels=true)
end
end

@testset "candidate scores" begin
    bandwidth = 9
    padding = 3

    function _test_candidate_scores(consensus, pseqs, scores, expected)
        ref_scores = Scores(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
        state = initial_state(consensus, pseqs, DNASeq(""), bandwidth, padding)
        rseq = RifrafSequence(DNASeq(), LogProb[], 1)
        redo_as = true
        redo_bs = true
        use_ref = false
        mult = 2
        recompute!(state, pseqs, scores, rseq, ref_scores, mult,
                   redo_as, redo_bs, 2, use_ref)
        do_alignment_proposals = true
        do_surgery_proposals = true
        trust_proposals = true
        indels_only = false
        cands = get_candidate_proposals(state, pseqs, scores,
                                        rseq, ref_scores,
                                        do_alignment_proposals,
                                        do_surgery_proposals,
                                        trust_proposals,
                                        indels_only)
        @test length(cands) == length(expected)
        if length(cands) == length(expected)
            cand_scores = sort([c.score for c in cands])
            exp_scores = sort([c.score for c in expected])
            @test all([a ≈ b for (a, b) in zip(cand_scores, exp_scores)])
        end
    end

    @testset "substitutions" begin
        consensus = DNASeq("TTT")
        seqs = [DNASeq("TAT")]
        lps =  Vector{LogProb}[fill(-1.0, length(s)) for s in seqs]

        pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                               for (s, p) in zip(seqs, lps)]
        scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))

        expected = [CandProposal(Substitution(2, DNA_A),
                                 sum(pseqs[1].match_log_p)),
                    CandProposal(Insertion(1, DNA_A),
                                 sum(pseqs[1].match_log_p) + pseqs[1].error_log_p[1] + scores.deletion),
                    CandProposal(Insertion(2, DNA_A),
                                 sum(pseqs[1].match_log_p) + pseqs[1].error_log_p[2] + scores.deletion),
                    CandProposal(Deletion(3),
                                 sum([pseqs[1].match_log_p[1],
                                      pseqs[1].error_log_p[2] + scores.insertion,
                                      pseqs[1].match_log_p[3]]))]
        _test_candidate_scores(consensus, pseqs, scores, expected)
    end

    @testset "deletion" begin
        consensus = DNASeq("TTT")
        seqs = [DNASeq("TT")]
        lps =  Vector{LogProb}[fill(-1.0, length(s)) for s in seqs]

        pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                               for (s, p) in zip(seqs, lps)]
        scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))

        expected = [CandProposal(Deletion(3),
                                 sum(pseqs[1].match_log_p))]
        _test_candidate_scores(consensus, pseqs, scores, expected)
    end

    @testset "insertion" begin
        # test insertion
        consensus = DNASeq("TT")
        seqs = [DNASeq("TAT")]
        lps =  Vector{LogProb}[fill(-1.0, length(s)) for s in seqs]

        pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                               for (s, p) in zip(seqs, lps)]
        scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))

        expected = [CandProposal(Insertion(2, DNA_A),
                                 sum(pseqs[1].match_log_p))]
        _test_candidate_scores(consensus, pseqs, scores, expected)
    end
end

@testset "full model" begin
    # can't guarantee this test actually passes, since it is random
    n_seqs = 6
    len = 30
    ref_error_rate = 0.1
    ref_sample_errors = ErrorModel(8.0, 0.0, 0.0, 1.0, 1.0)
    ref_errors = ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    template_error_mean = 0.01
    template_alpha = 1.0
    phred_scale = 3.0
    log_seq_actual_std = 3.0
    log_seq_reported_std = 0.3

    seq_errors = ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0)
    seq_scores = Scores(seq_errors)

    @testset "full model $i" for i in 1:100
        use_ref = rand([true, false])
        do_alignment_proposals = rand([true, false])
        do_surgery_proposals = rand([true, false])
        trust_proposals = rand([true, false])
        seed_indels = rand([true, false])
        fix_indels_stat = rand([true, false])
        indel_correction_only = rand([true, false])
        batch = rand([3, 6])

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
            reference = DNASeq()
        end

        (result, info) = rifraf(reads,
                                phreds, seq_scores;
                                reference=reference,
                                ref_scores=ref_scores,
                                do_alignment_proposals=do_alignment_proposals,
                                do_surgery_proposals=do_surgery_proposals,
                                trust_proposals=trust_proposals,
                                seed_indels=seed_indels,
                                fix_indels_stat=fix_indels_stat,
                                indel_correction_only=indel_correction_only,
                                bandwidth=10, min_dist=9, batch=batch,
                                max_iters=100)
        @test result == template
    end
end


@testset "base_probs" begin
    consensus = DNASeq("CGTAC")
    seqs = [DNASeq("CGAC"),
            DNASeq("CGAC"),
            DNASeq("CGAC")]
    lps = Vector{LogProb}[[-9.0, -9.0, -9.0, -9.0],
                               [-9.0, -9.0, -9.0, -9.0],
                               [-9.0, -9.0, -9.0, -9.0]]
    bandwidth = 6
    padding = 3
    mult = 2
    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)
    ref_errors = ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                           for (s, p) in zip(seqs, lps)]
    rseq = RifrafSequence(DNASeq(), LogProb[], bandwidth)
    state = initial_state(consensus, pseqs, rseq.seq, bandwidth, padding)
    recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    probs = estimate_probs(state, pseqs, scores, rseq, scores, false)
    @test probs.sub[1, 2] > 0.9
    @test probs.del[1] < 1e-9
    @test probs.del[3] > 0.9
end


@testset "ins_probs" begin
    consensus = DNASeq("CGAT")
    seqs = [DNASeq("CGTAT"),
            DNASeq("CGTAT"),
            DNASeq("CGTAT")]
    bandwidth = 6
    padding = 3
    mult = 2
    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)
    ref_errors = ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    lps = Vector{LogProb}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                           for (s, p) in zip(seqs, lps)]
    rseq = RifrafSequence(DNASeq(), LogProb[], bandwidth)
    state = initial_state(consensus, pseqs, rseq.seq, bandwidth, padding)
    recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    probs = estimate_probs(state, pseqs, scores, rseq, scores, false)
    @test maximum(probs.ins[1, :]) < 1e-9
    @test probs.ins[3, 4] > 0.9
end

@testset "alignment_probs" begin
    consensus = DNASeq("ACGT")
    seqs = [DNASeq("ACGT"),
            DNASeq("CGT"),
            DNASeq("CCGT")]
    bandwidth = 6
    padding = 3
    mult = 2

    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)
    ref_errors = ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    ps = [[0.1, 0.1, 0.1, 0.1],
          [0.2, 0.1, 0.1],
          [0.2, 0.1, 0.1, 0.1]]
    lps = map(log10, ps)
    pseqs = RifrafSequence[RifrafSequence(s, p, bandwidth)
                           for (s, p) in zip(seqs, lps)]
    rseq = RifrafSequence(DNASeq(), LogProb[], bandwidth)
    state = initial_state(consensus, pseqs, rseq.seq, bandwidth, padding)
    recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)

    result = alignment_error_probs(length(consensus),
                                   pseqs, state.Amoves)
    indices = sortperm(result)
    expected = [4, 3, 2, 1]
    @test indices == expected
end

@testset "align" begin
    const errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Scores(errors)

    @testset "align 1" begin
        template = DNASeq("ATAA")
        seq = DNASeq("AAA")
        bandwidth = 10
        log_p = [-2.0, -3.0, -3.0]
        pseq = RifrafSequence(seq, log_p, bandwidth)
        moves = align_moves(template, pseq, scores)
        t, s = moves_to_aligned_seqs(moves, template, seq)
        @test t == DNASequence("ATAA")
        @test s == DNASequence("A-AA")
    end

    @testset "align 2" begin
        template = DNASeq("AACCTT")
        seq = DNASeq("AAACCCTT")
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth)
        moves = align_moves(template, pseq, scores)
        t, s = moves_to_aligned_seqs(moves, template, seq)
        @test t[end-1:end] == DNASequence("TT")
    end
end
