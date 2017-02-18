using Quiver2.Proposals

using Bio.Seq
using Base.Test

srand(1)

@testset "Proposals" begin
    @testset "test_apply_proposals" begin
        template = dna"ACG"
        proposals = Proposal[Proposals.Insertion(0, DNA_T),
                             Proposals.Insertion(3, DNA_C),
                             Proposals.Deletion(3),
                             Proposals.Substitution(2, DNA_T)]
        expected = dna"TATC"
        result = Proposals.apply_proposals(template, proposals)
        @test result == expected
    end

    @testset "test_apply_ambiguous_insertions" begin
        template = dna"ACG"
        proposals = Proposal[Proposals.Insertion(0, DNA_C),
                             Proposals.Insertion(0, DNA_T)]
        @test_throws AmbiguousProposalsError Proposals.apply_proposals(template, proposals)
    end

    @testset "test_apply_ambiguous_proposals" begin
        template = dna"ACG"
        proposals = Proposal[Proposals.Deletion(1),
                             Proposals.Substitution(1, DNA_T)]
        @test_throws AmbiguousProposalsError Proposals.apply_proposals(template, proposals)
    end
end
