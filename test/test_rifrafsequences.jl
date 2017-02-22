using Base.Test

using Rifraf


@testset "RifrafSequence" begin
    @testset "scores" begin
        seq = DNASeq("ACGT")
        error_log_p = LogProb[-1., -2., -3., -4.]
        bandwidth = 10
        scores = Scores(-1., -2., -3., -4., -5.)
        rseq = RifrafSequence(seq, error_log_p, bandwidth, scores)
        
        @test rseq.match_scores == LogProb[-Inf; log10(1.0 - exp10(error_log_p)); -Inf]
        @test rseq.mismatch_scores == LogProb[-Inf; error_log_p; -Inf] + scores.mismatch
        @test rseq.ins_scores == LogProb[-Inf; error_log_p; -Inf] + scores.insertion
        @test rseq.del_scores == LogProb[-1., -1., -2., -3., -4.] + scores.deletion
        @test rseq.codon_ins_scores == LogProb[-Inf, -Inf, -Inf, -1., -2., -Inf, -Inf, -Inf] + scores.codon_insertion
        @test rseq.codon_del_scores == LogProb[-1., -1., -2., -3., -4.] + scores.codon_deletion
    end
end
