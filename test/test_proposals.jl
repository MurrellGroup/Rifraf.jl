using Quiver2.Proposals

using Base.Test

srand(1)

@testset "Proposals" begin
    @testset "test_apply_proposals" begin
        template = "ACG"
        proposals = Proposal[Proposals.Insertion(0, 'T'),
                             Proposals.Insertion(3, 'C'),
                             Proposals.Deletion(3),
                             Proposals.Substitution(2, 'T')]
        expected = "TATC"
        result = Proposals.apply_proposals(template, proposals)
        @test result == expected
    end

    @testset "test_apply_ambiguous_insertions" begin
        template = "ACG"
        proposals = Proposal[Proposals.Insertion(0, 'C'),
                             Proposals.Insertion(0, 'T')]
        @test_throws AmbiguousProposalsError Proposals.apply_proposals(template, proposals)
    end

    @testset "test_apply_ambiguous_proposals" begin
        template = "ACG"
        proposals = Proposal[Proposals.Deletion(1),
                             Proposals.Substitution(1, 'T')]
        @test_throws AmbiguousProposalsError Proposals.apply_proposals(template, proposals)
    end
end
