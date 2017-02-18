using Bio.Seq

using Quiver2.BandedArrays
using Quiver2.Sample
using Quiver2.Model
using Quiver2.Proposals
using Quiver2.Util

using Base.Test


srand(1234)


function inv_log10(logp::Float64)
    log10(1.0 - exp10(logp))
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
    @test Model.equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
    @test Model.equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
    @test Model.equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
end


@testset "test forward and backward" begin
    const errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Model.Scores(errors)

    @testset "perfect_forward" begin
        bandwidth = 1
        template = dna"AA"
        seq = dna"AA"
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

    @testset "imperfect_backward" begin
        bandwidth = 1
        template = dna"AA"
        seq = dna"AT"
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
        B = Model.backward(template, pseq, scores)
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
        template = dna"AA"
        seq = dna"AT"
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
        A = Model.forward(template, pseq, scores)
        B = Model.backward(template, pseq, scores)
        check_all_cols(A, B, false)
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, match + lp + scores.mismatch]],
                                     (3, 3)))

        @test full(A) == expected
    end

    @testset "forward/backward agreement" begin
        @testset "forward/backward agreement 1" begin
            template = dna"TG"
            seq = dna"GTCG"
            log_p = [-1.2, -0.8, -0.7, -1.0]
            bandwidth = 5
            pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
            local_scores = Model.Scores(Model.ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            A = Model.forward(template, pseq, local_scores)
            B = Model.backward(template, pseq, local_scores)
            check_all_cols(A, B, true)
        end

        @testset "forward/backward agreement" begin
            template = dna"GCACGGTC"
            seq = dna"GACAC"
            log_p = [-1.1, -1.1, -0.4, -1.0, -0.7]
            bandwidth = 5
            pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
            local_scores = Model.Scores(Model.ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            A = Model.forward(template, pseq, local_scores)
            B = Model.backward(template, pseq, local_scores)
            check_all_cols(A, B, true)
        end
    end

    @testset "insertion_agreement" begin
        template = dna"AA"
        seq = dna"ATA"
        bandwidth = 10
        log_p = [-5.0, -1.0, -6.0]
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
        A = Model.forward(template, pseq, scores)
        B = Model.backward(template, pseq, scores)
        score = (inv_log10(log_p[1]) +
                 log_p[2] + scores.insertion +
                 inv_log10(log_p[3]))
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)
    end

    @testset "deletion agreement" begin
        template = dna"GATAG"
        seq = dna"GAAG"
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
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)
    end

    @testset "deletion agreement 2" begin
        template = dna"ATA"
        seq = dna"AA"
        bandwidth = 10
        log_p = [-2.0, -3.0]
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
        A = Model.forward(template, pseq, scores)
        B = Model.backward(template, pseq, scores)
        score = (inv_log10(log_p[1]) +
                 maximum(log_p[1:2]) + scores.deletion +
                 inv_log10(log_p[2]))
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)
    end
end

@testset "scoring proposals" begin
    const errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Model.Scores(errors)

    function test_random_proposal(proposal, template_len)
        template = random_seq(template_len)

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

        (seq, actual, phreds, sbools,
         tbools) = sample_from_template(template,
                                        template_error_p,
                                        errors,
                                        phred_scale,
                                        log_actual_error_std,
                                        log_reported_error_std)
        log_p = Float64[Float64(q) / (-10.0) for q in phreds]
        bandwidth = max(5 * abs(length(template) - length(seq)), 30)
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)

        new_template = Proposals.apply_proposal(template, proposal)
        Anew = Model.forward(new_template, pseq, local_scores)
        Bnew = Model.backward(new_template, pseq, local_scores)
        check_all_cols(Anew, Bnew, codon_moves)

        A = Model.forward(template, pseq, local_scores)
        B = Model.backward(template, pseq, local_scores)
        check_all_cols(A, B, codon_moves)
        newcols = zeros(Float64, size(A)[1], 6)
        score = Model.seq_score_proposal(proposal, A, B, template, pseq,
                                         local_scores, newcols)
        @test_approx_eq_eps score Anew[end, end] 0.1
    end

    @testset "random_substitutions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(1:template_len)
            proposal = Proposals.Substitution(pos, rbase())
            test_random_proposal(proposal, template_len)
        end
    end

    @testset "random_insertions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(0:template_len)
            proposal = Proposals.Insertion(pos, rbase())
            test_random_proposal(proposal, template_len)
        end
    end

    @testset "random_deletions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(1:template_len)
            proposal = Proposals.Deletion(pos)
            test_random_proposal(proposal, template_len)
        end
    end
end

@testset "no single indels" begin
    reference = dna"AAAGGGTTT"
    ref_log_p = fill(log10(0.01), length(reference))
    bandwidth = 6
    rseq = Quiver2.Model.PString(reference, ref_log_p, bandwidth)
    local_errors = Model.normalize(Model.ErrorModel(2.0, 0.5, 0.5, 1.0, 1.0))
    local_scores = Model.Scores(local_errors)

    template = dna"AAACCCGGGTTT"
    @test !Quiver2.Model.has_single_indels(template, rseq, local_scores)

    template = dna"AAACCCGGGTTTT"
    @test Quiver2.Model.has_single_indels(template, rseq,
                                          local_scores)

    template = dna"AAA"
    @test !Quiver2.Model.has_single_indels(template, rseq,
                                           local_scores)
end


@testset "single_indel_proposals" begin
    ref_errors = Quiver2.Model.ErrorModel(10.0, 1e-10, 1e-10, 1.0, 1.0)
    ref_scores = Quiver2.Model.Scores(ref_errors)

    ref = dna"CGGCGATTT"
    consensus_errors = Float64[-8.04822,-5.10032,-5.09486,-1.0,-2.68901,-6.52537,-5.20094]
    consensus = Quiver2.Model.PString(dna"CTGCCGA", consensus_errors, 10)

    proposals = Quiver2.Model.single_indel_proposals(ref, consensus, ref_scores)
    expected = [Quiver2.Proposals.Deletion(4)]
    @test expected == proposals
end


@testset "fast proposals" begin
    function _test_fast_proposals(consensus, seqs, lps, expected;
                                  do_alignment::Bool=true,
                                  do_surgery::Bool=true,
                                  do_subs::Bool=true,
                                  do_indels::Bool=true)
        bandwidth = 5
        rseq = Quiver2.Model.PString(DNASequence(), Float64[], bandwidth)
        mult = 2
        errors = Model.ErrorModel(1.0, 5.0, 5.0, 0.0, 0.0)
        scores = Model.Scores(errors)
        ref_errors = Model.ErrorModel(10.0, 1e-10, 1e-10, 1.0, 1.0)
        ref_scores = Model.Scores(ref_errors)

        if length(lps) == 0
            lps =  Vector{Float64}[fill(-9.0, length(s)) for s in seqs]
        end

        pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                      for (s, p) in zip(seqs, lps)]
        state = Quiver2.Model.initial_state(consensus, pseqs)
        Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)

        if do_alignment
            proposals = Quiver2.Model.alignment_proposals(state, pseqs, do_subs, do_indels)
            @test length(symdiff(Set(proposals), Set(expected))) == 0
        end
        if do_surgery
            proposals, deltas = Quiver2.Model.surgery_proposals(state, pseqs, scores, do_subs, do_indels)
            @test length(symdiff(Set(proposals), Set(expected))) == 0
        end
    end

    @testset "fast proposals 1" begin
        consensus = dna"ACGAG"
        seqs = [dna"CGTAC",
                dna"CGAC",
                dna"CGTAG"]
        expected = [Proposals.Deletion(1),
                    Proposals.Insertion(3, DNA_T),
                    Proposals.Substitution(5, DNA_C)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals 2" begin
        consensus = dna"AA"
        seqs = [dna"AAG",
                dna"AA",
                dna"AAG"]
        expected = [Proposals.Insertion(2, DNA_G)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals 3" begin
        consensus = dna"AA"
        seqs = [dna"GAA",
                dna"AA",
                dna"GAA"]
        expected = [Proposals.Insertion(0, DNA_G)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals 4" begin
        consensus = dna"AA"
        seqs = [dna"AGA",
                dna"AA",
                dna"AGA"]
        expected = [Proposals.Insertion(1, DNA_G)]
        _test_fast_proposals(consensus, seqs, [], expected)
    end

    @testset "fast proposals - highly confident base" begin
        consensus = dna"AA"
        seqs = [dna"GAA",
                dna"AA",
                dna"AA"]
        lps = Vector{Float64}[[-10.0, -5.0, -5.0],
                              [-3.0, -5.0],
                              [-3.0, -50]]
        expected = [Proposals.Insertion(0, DNA_G)]
        _test_fast_proposals(consensus, seqs, lps, expected)
    end
    @testset "push deletions" begin
        # test that deletions get pushed
        consensus = dna"CCGTAAAC"
        seqs = [dna"CGTAAAC", dna"CCGTAAAC", dna"CGTAAAC"]
        lps = Vector{Float64}[[-1.0,-0.6,-1.1,-0.5,-0.5,-1.2,-2.0],
                              [-1.0,-1.0,-1.6,-2.2,-0.8,-0.5,-1.2,-1.8],
                              [-1.0,-1.8,-1.0,-0.5,-0.6,-1.0,-1.4]]
        expected = [Proposals.Deletion(2)]
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true)
    end
    @testset "push deletions - simple" begin
        # test that deletions get pushed to end
        consensus = dna"CCC"
        seqs = [dna"CC"]
        expected = [Proposals.Deletion(3)]
        _test_fast_proposals(consensus, seqs, [], expected,
                             do_alignment=false, do_surgery=true)
    end

    @testset "push insertions" begin
        # test that insertions get pushed to end
        consensus = dna"TT"
        seqs = [dna"TTT",
                dna"CTT"]
        expected = [Proposals.Insertion(0, DNA_C),
                    Proposals.Insertion(2, DNA_T)]
        _test_fast_proposals(consensus, seqs, [], expected,
                             do_alignment=false, do_surgery=true)
    end

    @testset "push subs" begin
        # test that subs get pushed backwards
        consensus = dna"ATG"
        seqs = [dna"ATG", dna"CATG", dna"CTG"]
        lps = [[-0.8, -1.7, -0.8],
               [-0.7, -0.5, -1.2, -1.0],
               [-0.9, -1.9, -1.0]]
        expected = [Proposals.Substitution(1, DNA_C)]
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true,
                             do_subs=true, do_indels=false)
    end

    @testset "test fast proposals converged" begin
        consensus = dna"GTTCGGCTC"
        seqs = [dna"GTTCGGCTTC",
                dna"GTTCGGCTC",
                dna"GTTCCTG"]
        phreds = Vector{Int8}[[28,16,13,21,15,13,13,12,20,16],
                              [21,16,9,17,6,15,6,16,12],
                              [26,14,5,24,8,12,7]]
        lps = map(Quiver2.Util.phred_to_log_p, phreds)
        expected = []
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true)

    end

    @testset "test fast proposals converged 2" begin
        consensus = dna"CAGTGCCGG"
        seqs = [dna"CATGCCGG",
                dna"CATGCCCTGG",
                dna"CAGGGCCGG"]
        phreds = Vector{Int8}[[13,7,12,13,7,11,6,14],
                              [16,14,14,20,5,5,15,12,10,20],
                              [23,9,7,6,9,10,10,10,23]]
        lps = map(Quiver2.Util.phred_to_log_p, phreds)
        expected = []
        # no indels, because this happened during refinement
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true,
                             do_subs=true, do_indels=false)
    end
    @testset "test fast proposals insertion/mismatch swap" begin
        # inserting 'T' after position 5 should swap with G/T
        # mismatch in second sequence
        consensus = dna"GGAAGTCC"
        seqs = [dna"GGAAGTCC",
                dna"GGAATTCC",
                dna"GGAAGTCTACC"]
        phreds = Vector{Int8}[[18,9,12,14,11,11,15,14],
                              [25,10,6,8,12,11,19,13],
                              [24,12,9,15,8,8,8,8,8,23,19]]
        lps = map(Quiver2.Util.phred_to_log_p, phreds)
        expected = [Proposals.Substitution(5, DNA_T),
                    Proposals.Insertion(6, DNA_T)]
        _test_fast_proposals(consensus, seqs, lps, expected,
                             do_alignment=false, do_surgery=true,
                             do_subs=true, do_indels=true)
    end
end

@testset "candidate scores" begin
    bandwidth = 10

    function _test_candidate_scores(template, pseqs, scores, expected)
        ref_scores = Model.Scores(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
        state = Quiver2.Model.initial_state(template, pseqs)
        rseq = Quiver2.Model.PString(DNASequence(), Float64[], 1)
        redo_as = true
        redo_bs = true
        use_ref = false
        mult = 2
        Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult,
                                 redo_as, redo_bs, 2, use_ref)
        do_alignment_proposals = true
        do_surgery_proposals = true
        trust_proposals = true
        indels_only = false
        cands = Quiver2.Model.get_candidate_proposals(state, pseqs, scores,
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
        template = dna"TTT"
        seqs = [dna"TAT"]
        lps =  Vector{Float64}[fill(-1.0, length(s)) for s in seqs]

        pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                      for (s, p) in zip(seqs, lps)]
        scores = Model.Scores(Model.ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))

        expected = [CandProposal(Proposals.Substitution(2, DNA_A),
                                 sum(pseqs[1].match_log_p)),
                    CandProposal(Proposals.Insertion(1, DNA_A),
                                 sum(pseqs[1].match_log_p) + pseqs[1].error_log_p[1] + scores.deletion),
                    CandProposal(Proposals.Insertion(2, DNA_A),
                                 sum(pseqs[1].match_log_p) + pseqs[1].error_log_p[2] + scores.deletion),
                    CandProposal(Proposals.Deletion(3),
                                 sum([pseqs[1].match_log_p[1],
                                      pseqs[1].error_log_p[2] + scores.insertion,
                                      pseqs[1].match_log_p[3]]))]
        _test_candidate_scores(template, pseqs, scores, expected)
    end

    @testset "deletion" begin
        template = dna"TTT"
        seqs = [dna"TT"]
        lps =  Vector{Float64}[fill(-1.0, length(s)) for s in seqs]

        pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                      for (s, p) in zip(seqs, lps)]
        scores = Model.Scores(Model.ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))

        expected = [CandProposal(Proposals.Deletion(3),
                                 sum(pseqs[1].match_log_p))]
        _test_candidate_scores(template, pseqs, scores, expected)
    end

    @testset "insertion" begin
        # test insertion
        template = dna"TT"
        seqs = [dna"TAT"]
        lps =  Vector{Float64}[fill(-1.0, length(s)) for s in seqs]

        pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                      for (s, p) in zip(seqs, lps)]
        scores = Model.Scores(Model.ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))

        expected = [CandProposal(Proposals.Insertion(2, DNA_A),
                                 sum(pseqs[1].match_log_p))]
        _test_candidate_scores(template, pseqs, scores, expected)
    end
end

@testset "full model" begin
    # can't guarantee this test actually passes, since it is random
    n_seqs = 5
    len = 30
    ref_error_rate = 0.1
    ref_sample_errors = Model.ErrorModel(8.0, 0.0, 0.0, 1.0, 1.0)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    template_error_mean = 0.01
    template_alpha = 1.0
    phred_scale = 3.0
    log_seq_actual_std = 3.0
    log_seq_reported_std = 0.3

    seq_errors = Model.ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0)
    seq_scores = Model.Scores(seq_errors)

    @testset "full model $i" for i in 1:100
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
            reference = DNASequence()
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
        @test result == template
    end
end


@testset "base_probs" begin
    template = dna"CGTAC"
    seqs = [dna"CGAC",
            dna"CGAC",
            dna"CGAC"]
    lps = Vector{Float64}[[-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0]]
    bandwidth = 5
    mult = 2
    errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Model.Scores(errors)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                  for (s, p) in zip(seqs, lps)]
    rseq = Quiver2.Model.PString(DNASequence(), Float64[], bandwidth)
    state = Quiver2.Model.initial_state(template, pseqs)
    Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    probs = Quiver2.Model.estimate_probs(state, pseqs, scores,
                                         rseq, scores, false)
    @test probs.sub[1, 2] > 0.9
    @test probs.del[1] < 1e-9
    @test probs.del[3] > 0.9
end


@testset "ins_probs" begin
    template = dna"CGAT"
    seqs = [dna"CGTAT",
            dna"CGTAT",
            dna"CGTAT"]
    bandwidth = 5
    mult = 2
    errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Model.Scores(errors)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    lps = Vector{Float64}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                  for (s, p) in zip(seqs, lps)]
    rseq = Quiver2.Model.PString(DNASequence(), Float64[], bandwidth)
    state = Quiver2.Model.initial_state(template, pseqs)
    Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)
    probs = Quiver2.Model.estimate_probs(state, pseqs, scores,
                                         rseq, scores, false)
    @test maximum(probs.ins[1, :]) < 1e-9
    @test probs.ins[3, 4] > 0.9
end

@testset "alignment_probs" begin
    consensus = dna"ACGT"
    seqs = [dna"ACGT",
            dna"CGT",
            dna"CCGT"]
    bandwidth = 5
    mult = 2

    errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Model.Scores(errors)
    ref_errors = Model.ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Model.Scores(ref_errors)

    ps = [[0.1, 0.1, 0.1, 0.1],
          [0.2, 0.1, 0.1],
          [0.2, 0.1, 0.1, 0.1]]
    lps = map(log10, ps)
    pseqs = Quiver2.Model.PString[Quiver2.Model.PString(s, p, bandwidth)
                                  for (s, p) in zip(seqs, lps)]
    rseq = Quiver2.Model.PString(DNASequence(), Float64[], bandwidth)
    state = Quiver2.Model.initial_state(consensus, pseqs)
    Quiver2.Model.recompute!(state, pseqs, scores, rseq, ref_scores, mult, true, true, 0, false)

    result = Quiver2.Model.alignment_error_probs(length(consensus),
                                                 pseqs, state.Amoves)
    indices = sortperm(result)
    expected = [4, 3, 2, 1]
    @test indices == expected
end

@testset "align" begin
    const errors = Model.normalize(Model.ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Model.Scores(errors)

    @testset "align 1" begin
        template = dna"ATAA"
        seq = dna"AAA"
        bandwidth = 10
        log_p = [-2.0, -3.0, -3.0]
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
        moves = Model.align_moves(template, pseq, scores)
        t, s = Model.moves_to_aligned_seqs(moves, template, seq)
        @test t == dna"ATAA"
        @test s == dna"A-AA"
    end

    @testset "align 2" begin
        template = dna"AACCTT"
        seq = dna"AAACCCTT"
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = Quiver2.Model.PString(seq, log_p, bandwidth)
        moves = Model.align_moves(template, pseq, scores)
        t, s = Model.moves_to_aligned_seqs(moves, template, seq)
        @test t[end-1:end] == dna"TT"
    end
end
