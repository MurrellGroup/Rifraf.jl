using Bio.Seq
using Iterators
using Base.Test

using Rifraf

include("test_utils.jl")

import Rifraf.LogProb,
       Rifraf.initial_state,
       Rifraf.resample!,
       Rifraf.realign_rescore!,
       Rifraf.estimate_consensus_error_probs,
       Rifraf.estimate_point_error_probs

@testset "estimate_consensus_probs_sub" begin
    consensus = DNASeq("CGAC")
    seqs = [DNASeq("CGAC"),
            DNASeq("CGTC"),
            DNASeq("CGTC")]
    lps = Vector{LogProb}[[-9.0, -9.0, -5.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0]]

    params = RifrafParams(bandwidth=6)
    errors = Rifraf.normalize(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    scores = Scores(errors)
    pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                           for (s, p) in zip(seqs, lps)]
    state = initial_state(consensus, pseqs, DNASeq(), params)
    resample!(state, params)
    realign_rescore!(state, RifrafParams())
    probs = estimate_consensus_error_probs(state, false)
    @test probs.sub[3, 4] > 0.9
    @test probs.sub[3, 1] > probs.sub[3, 2]

    point_probs = estimate_point_error_probs(probs)
    @test point_probs[1] == point_probs[2] == point_probs[4]
    @test point_probs[1] < 0.1
    @test point_probs[3] < 0.9
end


@testset "estimate_consensus_probs_del" begin
    consensus = DNASeq("CGTAC")
    seqs = [DNASeq("CGAC"),
            DNASeq("CGAC"),
            DNASeq("CGAC")]
    lps = Vector{LogProb}[[-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0]]

    params = RifrafParams(bandwidth=6)
    errors = Rifraf.normalize(ErrorModel(1.0, 3.0, 3.0, 0.0, 0.0))
    scores = Scores(errors)
    pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                           for (s, p) in zip(seqs, lps)]
    state = initial_state(consensus, pseqs, DNASeq(), params)
    resample!(state, params)
    realign_rescore!(state, RifrafParams())
    probs = estimate_consensus_error_probs(state, false)
    @test probs.sub[1, 2] > 0.9
    @test probs.del[1] < 1e-9
    @test probs.del[3] > 0.9
    point_probs = estimate_point_error_probs(probs)
    @test point_probs[1] == point_probs[5]
    @test point_probs[2] == point_probs[4]
    @test point_probs[1] < 0.1
    @test point_probs[2] < 0.1
    @test point_probs[3] > 0.9
end

@testset "estimate_consensus_probs_ins" begin
    consensus = DNASeq("CGAT")
    seqs = [DNASeq("CGTAT"),
            DNASeq("CGTAT"),
            DNASeq("CGTAT")]
    params = RifrafParams(bandwidth=6)
    errors = Rifraf.normalize(ErrorModel(1.0, 3.0, 3.0, 0.0, 0.0))
    scores = Scores(errors)

    lps = Vector{LogProb}[[-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0],
                          [-9.0, -9.0, -9.0, -9.0, -9.0]]
    pseqs = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                           for (s, p) in zip(seqs, lps)]
    state = initial_state(consensus, pseqs, DNASeq(), params)
    resample!(state, params)
    realign_rescore!(state, RifrafParams())
    probs = estimate_consensus_error_probs(state, false)
    @test sum(probs.ins[1, :]) <= 1e-3
    @test probs.ins[3, 4] > 0.9

    point_probs = estimate_point_error_probs(probs)
    @test point_probs[1] == point_probs[4]
    @test point_probs[2] == point_probs[3]
    @test point_probs[1] < 0.1
    @test point_probs[2] > 0.9
end
