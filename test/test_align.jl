using Bio.Seq
using Base.Test

using Rifraf

include("test_utils.jl")

import Rifraf.forward_moves,
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
        A, _ = forward_moves(template, pseq)
        # transpose because of column-major order
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, 2 * match]],
                                     (3, 3)))
        @test full(A) ≈ expected
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
        A, _ = forward_moves(template, pseq)
        B = backward(template, pseq)
        @test check_all_cols(A, B, false)
        expected = transpose(reshape([[0.0, lp + scores.deletion, 0.0];
                                      [lp + scores.insertion, match, match + lp + scores.deletion];
                                      [0.0, match + lp + scores.insertion, match + lp + scores.mismatch]],
                                     (3, 3)))

        @test full(A) ≈ expected
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
        @test full(B) ≈ expected
    end

    @testset "forward/backward agreement" begin
        @testset "forward/backward agreement 1" begin
            template = DNASeq("TG")
            seq = DNASeq("GTCG")
            log_p = [-1.2, -0.8, -0.7, -1.0]
            bandwidth = 5
            local_scores = Scores(ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            pseq = RifrafSequence(seq, log_p, bandwidth, local_scores)
            A, _ = forward_moves(template, pseq)
            B = backward(template, pseq)
            @test check_all_cols(A, B, true)
            @test true
        end

        @testset "forward/backward agreement 2" begin
            template = DNASeq("GCACGGTC")
            seq = DNASeq("GACAC")
            log_p = [-1.1, -1.1, -0.4, -1.0, -0.7]
            bandwidth = 5
            local_scores = Scores(ErrorModel(2.0, 1.0, 1.0, 3.0, 3.0))
            pseq = RifrafSequence(seq, log_p, bandwidth, local_scores)
            A, _ = forward_moves(template, pseq)
            B = backward(template, pseq)
            @test check_all_cols(A, B, true)
        end
    end

    @testset "insertion_agreement" begin
        template = DNASeq("AA")
        seq = DNASeq("ATA")
        bandwidth = 10
        log_p = [-5.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A, _ = forward_moves(template, pseq)
        B = backward(template, pseq)
        score = (inv_log10(log_p[1]) +
                 log_p[2] + scores.insertion +
                 inv_log10(log_p[3]))
        @test A[end, end] ≈ score
        @test check_all_cols(A, B, false)
    end

    @testset "deletion agreement 1" begin
        template = DNASeq("GATAG")
        seq = DNASeq("GAAG")
        bandwidth = 10
        log_p = [-5.0, -2.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A, _ = forward_moves(template, pseq)
        B = backward(template, pseq)
        score = (pseq.match_scores[1] +
                 pseq.match_scores[2] +
                 pseq.del_scores[3] +
                 pseq.match_scores[3] +
                 pseq.match_scores[4])
        @test A[end, end] ≈ score
        @test check_all_cols(A, B, false)
    end

    @testset "deletion agreement 2" begin
        template = DNASeq("ATA")
        seq = DNASeq("AA")
        bandwidth = 10
        log_p = [-2.0, -3.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A, _ = forward_moves(template, pseq)
        B = backward(template, pseq)
        score = (pseq.match_scores[1] +
                 pseq.del_scores[2] +
                 pseq.match_scores[2])
        @test A[end, end] ≈ score
        @test check_all_cols(A, B, false)
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
