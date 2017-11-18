using BioSymbols
using BioSequences
using Base.Test

using Rifraf

include("test_utils.jl")

import Rifraf.forward,
       Rifraf.forward_moves,
       Rifraf.backward


# TODO: some codon tests
# TODO: test all update() variants
@testset "test forward and backward" begin
    const scores = Scores(-1.0, -1.0, -1.0, -Inf, -Inf)

    @testset "perfect forward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AA")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A = forward(template, pseq)
        # transpose because of column-major order
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, 2 * match]],
                                     (3, 3)))
        @test full(A) ≈ expected

        A2, _ = forward_moves(template, pseq)
        @test full(A2) ≈ full(A)
    end

    @testset "perfect backward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AA")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A = backward(template, pseq)
        # transpose because of column-major order
        expected = transpose(reshape([[2 * match, match + lp + scores.insertion, 0.0];
                                      [match + lp + scores.deletion, match, lp + scores.insertion];
                                      [0.0, lp + scores.deletion, 0.0]],
                                     (3, 3)))
        @test full(A) ≈ expected
    end

    @testset "imperfect forward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AT")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A1 = forward(template, pseq)
        B = backward(template, pseq)
        check_all_cols(A1, B, false)
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, match + lp + scores.mismatch]],
                                     (3, 3)))

        @test full(A1) ≈ expected atol=0.01

        A2, moves = forward_moves(template, pseq)
        @test full(A1) ≈ full(A2) atol=0.01
    end

    @testset "imperfect backward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AT")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        B = backward(template, pseq)
        expected = transpose(reshape([[lp + scores.mismatch + match, lp + scores.insertion + match, 0.0];
                                      [2*lp + scores.deletion + scores.mismatch, lp + scores.mismatch, lp + scores.insertion];
                                      [0.0, lp + scores.deletion, 0.0]],
                                     (3, 3)))
        @test full(B) ≈ expected atol=0.01
    end

    @testset "forward/backward agreement" begin
        @testset "forward/backward agreement 1" begin
            template = DNASeq("TG")
            seq = DNASeq("GTCG")
            log_p = [-1.2, -0.8, -0.7, -1.0]
            bandwidth = 5
            local_scores = Scores(ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            pseq = RifrafSequence(seq, log_p, bandwidth, local_scores)
            A = forward(template, pseq)
            B = backward(template, pseq)
            check_all_cols(A, B, true)

            A2, _ = forward_moves(template, pseq)
            @test full(A) ≈ full(A2) atol=0.01
        end

        @testset "forward/backward agreement 2" begin
            template = DNASeq("GCACGGTC")
            seq = DNASeq("GACAC")
            log_p = [-1.1, -1.1, -0.4, -1.0, -0.7]
            bandwidth = 5
            local_scores = Scores(ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            pseq = RifrafSequence(seq, log_p, bandwidth, local_scores)
            A = forward(template, pseq)
            B = backward(template, pseq)
            check_all_cols(A, B, true)

            A2, _ = forward_moves(template, pseq)
            @test full(A) ≈ full(A2) atol=0.01
        end
    end

    @testset "insertion_agreement" begin
        template = DNASeq("AA")
        seq = DNASeq("ATA")
        bandwidth = 10
        log_p = [-5.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A = forward(template, pseq)
        B = backward(template, pseq)
        score = (inv_log10(log_p[1]) +
                 log_p[2] + scores.insertion +
                 inv_log10(log_p[3]))
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)

        A2, _ = forward_moves(template, pseq)
        @test full(A) == full(A2)
    end

    @testset "deletion agreement 1" begin
        template = DNASeq("GATAG")
        seq = DNASeq("GAAG")
        bandwidth = 10
        log_p = [-5.0, -2.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A = forward(template, pseq)
        B = backward(template, pseq)
        score = (pseq.match_scores[1] +
                 pseq.match_scores[2] +
                 pseq.del_scores[3] +
                 pseq.match_scores[3] +
                 pseq.match_scores[4])
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)

        A2, _ = forward_moves(template, pseq)
        @test full(A) == full(A2)
    end

    @testset "deletion agreement 2" begin
        template = DNASeq("ATA")
        seq = DNASeq("AA")
        bandwidth = 10
        log_p = [-2.0, -3.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A = forward(template, pseq)
        B = backward(template, pseq)
        score = (pseq.match_scores[1] +
                 pseq.del_scores[2] +
                 pseq.match_scores[2])
        @test A[end, end] ≈ score
        check_all_cols(A, B, false)

        A2, _ = forward_moves(template, pseq)
        @test full(A) == full(A2)
    end
end

@testset "alignment" begin
    const errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    const scores = Scores(errors)

    @testset "align 1" begin
        template = DNASeq("ATAA")
        seq = DNASeq("AAA")
        bandwidth = 10
        log_p = [-2.0, -3.0, -3.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        moves = Rifraf.align_moves(template, pseq)
        t, s = Rifraf.moves_to_aligned_seqs(moves, template, seq)
        @test t == DNASequence("ATAA")
        @test s == DNASequence("A-AA")
    end

    @testset "align 2" begin
        template = DNASeq("AACCTT")
        seq = DNASeq("AAACCCTT")
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        moves = Rifraf.align_moves(template, pseq)
        t, s = Rifraf.moves_to_aligned_seqs(moves, template, seq)
        @test t[end-1:end] == DNASequence("TT")
    end

    @testset "moves_to_indices 1" begin
        template = DNASeq("AAA")
        seq = DNASeq("AAA")
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        moves = Rifraf.align_moves(template, pseq)
        indices = Rifraf.moves_to_indices(moves, length(template), length(pseq))
        @test indices == collect(1:3)
    end

    @testset "moves_to_indices 2" begin
        template = DNASeq("AAA")
        seq = DNASeq("AAAT")
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        moves = Rifraf.align_moves(template, pseq)
        indices = Rifraf.moves_to_indices(moves, length(template), length(pseq))
        @test indices == collect(1:3)
    end

    @testset "moves_to_indices 3" begin
        template = DNASeq("AAAT")
        seq = DNASeq("AAA")
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        moves = Rifraf.align_moves(template, pseq)
        indices = Rifraf.moves_to_indices(moves, length(template), length(pseq))
        @test indices == [1, 2, 3, 3]
    end

    @testset "moves_to_indices 4" begin
        template = DNASeq("TAAA")
        seq = DNASeq("AAA")
        bandwidth = 10
        log_p = fill(log10(0.1), length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        moves = Rifraf.align_moves(template, pseq)
        indices = Rifraf.moves_to_indices(moves, length(template), length(pseq))
        @test indices == [0, 1, 2, 3]
    end

    @testset "align and skew" begin
        ref_errors = ErrorModel(10.0, 1e-10, 1e-10, 1.0, 1.0)
        ref_scores = Scores(ref_errors)

        ref = DNASeq("CGGCGATTT")
        consensus_errors = LogProb[-8., -8., -8., -1., -8., -10., -10.]
        consensus = RifrafSequence(DNASeq("CTGCCGA"), consensus_errors, 10, ref_scores)
        a, b = Rifraf.align(ref, consensus, skew_matches=true)
        expected_a = dna"CGG-CGATTT"
        expected_b = dna"CTGCCGA---"
        @test a == expected_a
        @test b == expected_b
    end
end

@testset "align with self" begin
    # make sequence with long insertion
    seqstring = DNASequence("AAAGGGTTTCCC")
    seq = DNASeq(seqstring)
    errors = fill(0.1, length(seq))
    errors[1:6] = 0.3
    errors[end-3:end] = 0.45
    err_log_p = log10.(errors)
    bandwidth = 3
    scores = Scores(ErrorModel(1.0, 10.0, 10.0, 0.0, 0.0))
    rseq = RifrafSequence(seq, err_log_p, bandwidth, scores)

    a, b = Rifraf.align(seq, rseq)
    @test a == b
    @test a == seqstring
end
