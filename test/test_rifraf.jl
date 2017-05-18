using Bio.Seq
using Iterators
using Base.Test

using Rifraf

include("test_utils.jl")

import Rifraf.sample_from_template,
       Rifraf.random_seq,
       Rifraf.ScoredProposal,
       Rifraf.Proposal,
       Rifraf.Substitution,
       Rifraf.Insertion,
       Rifraf.Deletion,
       Rifraf.apply_proposals,
       Rifraf.rbase,
       Rifraf.forward_moves,
       Rifraf.backward,
       Rifraf.phred_to_log_p,
       Rifraf.equal_ranges,
       Rifraf.initial_state,
       Rifraf.estimate_probs,
       Rifraf.has_single_indels,
       Rifraf.score_proposal,
       Rifraf.single_indel_proposals,
       Rifraf.realign_and_score!,
       Rifraf.get_candidates,
       Rifraf.alignment_proposals,
       Rifraf.surgery_candidates,
       Rifraf.align_moves,
       Rifraf.moves_to_aligned_seqs,
       Rifraf.alignment_error_probs,
       Rifraf.Prob,
       Rifraf.LogProb

srand(1234)


@testset "equal_ranges" begin
    @test equal_ranges((3, 5), (4, 6)) == ((2, 3), (1, 2))
    @test equal_ranges((1, 5), (1, 2)) == ((1, 2), (1, 2))
    @test equal_ranges((1, 5), (4, 5)) == ((4, 5), (1, 2))
end


@testset "scoring proposals" begin
    const errors = ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0)
    const scores = Scores(errors)

    function test_random_proposal(proposal, template_len)
        template = random_seq(template_len)

        template_error_p = Prob[0.1 for _=1:length(template)]
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
        pseq = RifrafSequence(seq, log_p, bandwidth, local_scores)

        new_template = apply_proposals(template, Proposal[proposal])
        new_template_error_p = Prob[0.1 for _=1:length(new_template)]

        Anew = Rifraf.forward(new_template, pseq)
        Bnew = backward(new_template, pseq)
        check_all_cols(Anew, Bnew, codon_moves)

        A = Rifraf.forward(template, pseq)
        B = backward(template, pseq)
        check_all_cols(A, B, codon_moves)

        newcols = zeros(size(A)[1], Rifraf.CODON_LENGTH + 1)
        score = score_proposal(proposal, A, B, template, pseq, newcols)

        @test score ≈ Anew[end, end]
    end

    @testset "random substitutions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(1:template_len)
            proposal = Substitution(pos, rbase())
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random insertions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(0:template_len)
            proposal = Insertion(pos, rbase())
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random deletions" begin
        for i = 1:1000
            template_len = rand(30:50)
            pos = rand(1:template_len)
            proposal = Deletion(pos)
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random deletions at beginning" begin
        for i = 1:10
            template_len = rand(30:50)
            proposal = Deletion(1)
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random deletions at end" begin
        for i = 1:10
            template_len = rand(30:50)
            proposal = Deletion(template_len)
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random substitutions at beginning" begin
        # subs at the beginning
        for i = 1:10
            template_len = rand(30:50)
            proposal = Substitution(1, rbase())
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random substitutions at end" begin
        # subs at the end
        for i = 1:10
            template_len = rand(30:50)
            proposal = Substitution(template_len, rbase())
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random insertions at beginning" begin
        for i = 1:10
            template_len = rand(30:50)
            proposal = Insertion(0, rbase())
            test_random_proposal(proposal, template_len)
        end
    end
    @testset "random insertions at end" begin
        for i = 1:10
            template_len = rand(30:50)
            proposal = Insertion(template_len, rbase())
            test_random_proposal(proposal, template_len)
        end
    end
end

@testset "no single indels" begin
    reference = DNASeq("AAAGGGTTT")
    ref_log_p = fill(log10(0.01), length(reference))
    bandwidth = 6
    local_errors = Rifraf.normalize(ErrorModel(2.0, 0.5, 0.5, 1.0, 1.0))
    local_scores = Scores(local_errors)
    rseq = RifrafSequence(reference, ref_log_p, bandwidth, local_scores)

    template = DNASeq("AAACCCGGGTTT")
    @test !has_single_indels(template, rseq)

    template = DNASeq("AAACCCGGGTTTT")
    @test has_single_indels(template, rseq)

    template = DNASeq("AAA")
    @test !has_single_indels(template, rseq)
end


@testset "single indel proposals" begin
    ref_errors = ErrorModel(10.0, 1e-10, 1e-10, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    ref = DNASeq("CGGCGATTT")
    consensus_errors = LogProb[-8., -8., -8., -1., -8., -10., -10.]
    consensus = RifrafSequence(DNASeq("CTGCCGA"), consensus_errors, 10, ref_scores)

    proposals = single_indel_proposals(ref, consensus)
    expected = [Deletion(4)]
    @test expected == proposals
end

@testset "surgery helpers" begin
    @testset "surrounding insertions" begin
        moves = [Rifraf.TRACE_INSERT,
                 Rifraf.TRACE_MATCH,
                 Rifraf.TRACE_INSERT]
        scores = Scores(ErrorModel(1.0, 5.0, 5.0, 0.0, 0.0))

        # Will want to make the base with higher error prob the insertion
        # so should pick 1st
        seq = RifrafSequence(DNASeq("TAT"), LogProb[-1.0, -2.0, -3.0], 10, scores)
        aln_idx = 2
        seq_idx = 2
        result = Rifraf.best_surrounding_ins_bases(moves, seq, aln_idx, seq_idx)
        @test length(result) == 1
        @test first(keys(result)) == DNA_T
        @test result[DNA_T] == (seq.match_scores[3], seq.ins_scores[3])
    end

    @testset "surrounding deletions" begin
        moves = [Rifraf.TRACE_MATCH,
                 Rifraf.TRACE_DELETE,
                 Rifraf.TRACE_MATCH,
                 Rifraf.TRACE_DELETE,
                 Rifraf.TRACE_MATCH]
        scores = Scores(ErrorModel(1.0, 5.0, 5.0, 0.0, 0.0))

        # gap will want to move towards base with higher log error prob (-1.0)
        consensus = DNASeq("ATGTA")
        seq = RifrafSequence(DNASeq("ATA"), LogProb[-2.0, -10.0, -1.0], 10, scores)
        aln_idx = 3
        cons_idx = 3
        seq_idx = 2
        result = Rifraf.surrounding_del_bases(moves, consensus, seq, aln_idx, cons_idx, seq_idx)
        @test length(result) == 1
        @test first(keys(result)) == DNA_T
        @test result[DNA_T] == (seq.del_scores[2], seq.del_scores[3])
    end
end

@testset "fast proposals" begin
    function _test_fast_proposals(consensus, seqs, lps, expected;
                                  do_alignment::Bool=true,
                                  do_surgery::Bool=true,
                                  do_indels::Bool=true)
        params = RifrafParams(bandwidth=6)

        errors = ErrorModel(1.0, 5.0, 5.0, 0.0, 0.0)
        scores = Scores(errors)
        if length(lps) == 0
            lps =  Vector{LogProb}[fill(-9.0, length(s)) for s in seqs]
        end
        pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                               for (s, p) in zip(seqs, lps)]
        state = initial_state(consensus, pseqs, DNASeq(), params)
        realign_and_score!(state, false)

        if do_alignment
            proposals = alignment_proposals(state.Amoves, state.consensus,
                                            pseqs, do_indels)
            @test length(symdiff(Set(proposals), Set(expected))) == 0
        end
        if do_surgery
            cands = surgery_candidates(state, do_indels)
            proposals = Proposal[c.proposal for c in cands]
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
                             do_indels=false)
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
                         do_indels=false)
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
                         do_indels=true)
end
end

@testset "candidate scores" begin
    params = RifrafParams(bandwidth=9)

    function _test_candidate_scores(consensus, pseqs, expected)
        state = initial_state(consensus, pseqs, DNASeq(), params)
        realign_and_score!(state, false)
        do_alignment_proposals = true
        do_surgery_proposals = false
        indels_only = false
        cands = get_candidates(state,
                               do_alignment_proposals,
                               do_surgery_proposals,
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

        scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))
        pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                               for (s, p) in zip(seqs, lps)]

        expected = [ScoredProposal(Substitution(2, DNA_A),
                                   sum(pseqs[1].match_scores))]
        _test_candidate_scores(consensus, pseqs, expected)
    end

    @testset "deletion" begin
        consensus = DNASeq("TTT")
        seqs = [DNASeq("TT")]
        lps =  Vector{LogProb}[fill(-1.0, length(s)) for s in seqs]

        scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))
        pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                               for (s, p) in zip(seqs, lps)]

        expected = [ScoredProposal(Deletion(3),
                                   sum(pseqs[1].match_scores))]
        _test_candidate_scores(consensus, pseqs, expected)
    end

    @testset "insertion" begin
        # test insertion
        consensus = DNASeq("TT")
        seqs = [DNASeq("TAT")]
        lps =  Vector{LogProb}[fill(-1.0, length(s)) for s in seqs]

        scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))
        pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                               for (s, p) in zip(seqs, lps)]

        expected = [ScoredProposal(Insertion(2, DNA_A),
                                   sum(pseqs[1].match_scores))]
        _test_candidate_scores(consensus, pseqs, expected)
    end
end

@testset "full model" begin
    # can't guarantee this test actually passes, since it is random
    n_seqs = 5
    len = 30

    ref_sample_errors = ErrorModel(8.0, 0.0, 0.0, 1.0, 1.0)
    ref_errors = ErrorModel(8.0, 0.1, 0.1, 1.0, 1.0)
    ref_scores = Scores(ref_errors)

    seq_errors = ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0)
    seq_scores = Scores(seq_errors)

    sample_params = Dict(:ref_error_rate => 0.1,
                         :ref_errors => ref_sample_errors,
                         :error_rate => 0.005,
                         :alpha => 1.0,
                         :phred_scale => 1.5,
                         :actual_std => 3.0,
                         :reported_std => 0.3,
                         :seq_errors => seq_errors)

    @testset "full model" begin
        for (use_ref,
             do_alignment_proposals,
             do_surgery_proposals,
             seed_indels,
             indel_correction_only,
             batch_size) in product([true, false],
                                    [true, false],
                                    [true, false],
                                    [true, false],
                                    [true, false],
                                    [3, 6])

            if do_alignment_proposals && do_surgery_proposals
                continue
            end

            (reference, template, template_error, reads, actual, phreds, sbools,
             tbools) = sample(n_seqs, len; sample_params...)
            if !use_ref
                reference = DNASeq()
            end

            params = RifrafParams(ref_scores=ref_scores,
                                  do_alignment_proposals=do_alignment_proposals,
                                  do_surgery_proposals=do_surgery_proposals,
                                  seed_indels=seed_indels,
                                  indel_correction_only=indel_correction_only,
                                  batch_size=batch_size)

            result = rifraf(reads,
                            phreds, seq_scores;
                            reference=reference,
                            params=params)
            @test result.consensus == template
        end
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

    params = RifrafParams(bandwidth=6)
    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)
    pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                           for (s, p) in zip(seqs, lps)]
    state = initial_state(consensus, pseqs, DNASeq(), params)
    realign_and_score!(state, false)
    probs = estimate_probs(state, false)
    @test probs.sub[1, 2] > 0.9
    @test probs.del[1] < 1e-9
    @test probs.del[3] > 0.9
end


@testset "ins_probs" begin
    consensus = DNASeq("CGAT")
    seqs = [DNASeq("CGTAT"),
            DNASeq("CGTAT"),
            DNASeq("CGTAT")]
    params = RifrafParams(bandwidth=6)
    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)

    lps = Vector{LogProb}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                           for (s, p) in zip(seqs, lps)]
    state = initial_state(consensus, pseqs, DNASeq(), params)
    realign_and_score!(state, false)
    probs = estimate_probs(state, false)
    @test maximum(probs.ins[1, :]) < 1e-9
    @test probs.ins[3, 4] > 0.9
end

@testset "alignment_probs" begin
    consensus = DNASeq("ACGT")
    seqs = [DNASeq("ACGT"),
            DNASeq("CGT"),
            DNASeq("CCGT")]

    params = RifrafParams(bandwidth=6)
    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)

    ps = [[0.1, 0.1, 0.1, 0.1],
          [0.2, 0.1, 0.1],
          [0.2, 0.1, 0.1, 0.1]]
    lps = map(log10, ps)
    pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                           for (s, p) in zip(seqs, lps)]
    state = initial_state(consensus, pseqs, DNASeq(), params)
    realign_and_score!(state, false)

    result = alignment_error_probs(length(consensus),
                                   pseqs, state.Amoves)
    indices = sortperm(result)
    expected = [4, 3, 2, 1]
    @test indices == expected
end
