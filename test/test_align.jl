using Bio.Seq
using Base.Test

using Rifraf

include("test_utils.jl")

import Rifraf.forward_moves,
       Rifraf.backward


# TODO: some codon tests
@testset "test forward and backward" begin
    const scores = Scores(-1.0, -1.0, -1.0, -Inf, -Inf)

    @testset "perfect_forward" begin
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

    @testset "imperfect_backward" begin
        bandwidth = 1
        template = DNASeq("AA")
        seq = DNASeq("AT")
        lp = -3.0
        match = inv_log10(lp)
        log_p = fill(lp, length(seq))
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        B = backward(template, pseq)
        expected = transpose(reshape([[lp + scores.mismatch + match,
                                       lp + scores.insertion + match, 0.0];
                                      [2*lp + scores.deletion + scores.mismatch,
                                       lp + scores.mismatch,
                                       lp + scores.insertion];
                                      [0.0, lp + scores.deletion, 0.0]],
                                     (3, 3)))
        @test full(B) ≈ expected
    end

    @testset "imperfect_forward" begin
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

        @testset "forward/backward agreement" begin
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

    @testset "deletion agreement" begin
        template = DNASeq("GATAG")
        seq = DNASeq("GAAG")
        bandwidth = 10
        log_p = [-5.0, -2.0, -1.0, -6.0]
        pseq = RifrafSequence(seq, log_p, bandwidth, scores)
        A, _ = forward_moves(template, pseq)
        B = backward(template, pseq)
        score = (inv_log10(log_p[1]) +
                 inv_log10(log_p[2]) +
                 maximum(log_p[2:3]) + scores.deletion +
                 inv_log10(log_p[3]) +
                 inv_log10(log_p[4]))
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
        score = (inv_log10(log_p[1]) +
                 maximum(log_p[1:2]) + scores.deletion +
                 inv_log10(log_p[2]))
        @test A[end, end] ≈ score
        @test check_all_cols(A, B, false)
    end
end
