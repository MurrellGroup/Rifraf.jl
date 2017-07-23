using Base.Test

using Rifraf


@testset "RifrafSequence" begin
    @testset "empty" begin
        bandwidth = 10
        scores = Scores(-1., -2., -3., -4., -5.)
        rseq = RifrafSequence(DNASeq(), LogProb[], bandwidth, scores)
        @test length(rseq) == 0
        @test length(RifrafSequence()) == 0
    end

    @testset "scores" begin
        seq = DNASeq("ACGT")
        error_log_p = LogProb[-1., -2., -3., -4.]
        bandwidth = 10
        scores = Scores(-1., -2., -3., -4., -5.)
        rseq = RifrafSequence(seq, error_log_p, bandwidth, scores)

        @test rseq.match_scores == log10.(1.0 - exp10.(error_log_p))
        @test rseq.mismatch_scores == error_log_p + scores.mismatch
        @test rseq.ins_scores == error_log_p + scores.insertion
        @test rseq.del_scores == LogProb[-1., -1., -2., -3., -4.] + scores.deletion
        @test rseq.codon_ins_scores == LogProb[-1., -2.] + scores.codon_insertion
        @test rseq.codon_del_scores == LogProb[-1., -1., -2., -3., -4.] + scores.codon_deletion
    end

    @testset "phred" begin
        seq = DNASeq("ACGT")
        phreds = Int8[3, 50, 10, 70]
        bandwidth = 10
        scores = Scores(-1., -2., -3., -4., -5.)

        # error log p
        rseq = RifrafSequence(seq, phreds, bandwidth, scores)
        # TODO: test
    end

    @testset "update scores" begin
        seq = DNASeq("ACGT")
        error_log_p = LogProb[-1., -2., -3., -4.]
        bandwidth = 10
        scores = Scores(-1., -2., -3., -4., -5.)
        rseq = RifrafSequence(seq, error_log_p, bandwidth, scores)

        new_scores = Scores(-1., -1., -1., -1., -1.)
        new_rseq = RifrafSequence(rseq, new_scores)
        @test new_rseq.ins_scores == new_rseq.mismatch_scores
    end
end
