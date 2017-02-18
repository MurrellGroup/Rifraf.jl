using Bio.Seq
using Base.Test

using Rifraf

import Rifraf.Proposal,
       Rifraf.Substitution,
       Rifraf.Insertion,
       Rifraf.Deletion,
       Rifraf.apply_proposals

srand(1)

@testset "Proposals" begin
    @testset "test_apply_proposals" begin
        template = dna"ACG"
        proposals = Proposal[Insertion(0, DNA_T),
                             Insertion(3, DNA_C),
                             Deletion(3),
                             Substitution(2, DNA_T)]
        expected = dna"TATC"
        result = apply_proposals(template, proposals)
        @test result == expected
    end

    @testset "test_apply_ambiguous_insertions" begin
        template = dna"ACG"
        proposals = Proposal[Insertion(0, DNA_C),
                             Insertion(0, DNA_T)]
        @test_throws Rifraf.AmbiguousProposalsError apply_proposals(template, proposals)
    end

    @testset "test_apply_ambiguous_proposals" begin
        template = dna"ACG"
        proposals = Proposal[Deletion(1),
                             Substitution(1, DNA_T)]
        @test_throws Rifraf.AmbiguousProposalsError apply_proposals(template, proposals)
    end
end
