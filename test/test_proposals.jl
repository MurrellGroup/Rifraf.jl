using BioSymbols
using BioSequences
using Base.Test

using Rifraf

import Rifraf.Proposal,
       Rifraf.Substitution,
       Rifraf.Insertion,
       Rifraf.Deletion,
       Rifraf.apply_proposals

srand(1)

@testset "Proposals" begin
    @testset "apply proposals" begin
        template = DNASeq("ACG")
        proposals = Proposal[Insertion(0, DNA_T),
                             Substitution(2, DNA_T),
                             Insertion(3, DNA_A),
                             Deletion(3)]
        expected = DNASeq("TATA")
        result = apply_proposals(template, proposals)
        @test result == expected
    end

    @testset "apply insertion" begin
        template = DNASeq("ATT")
        proposals = Proposal[Rifraf.Insertion(2,DNA_G)]
        expected = DNASeq("ATGT")
        result = apply_proposals(template, proposals)
        @test result == expected
    end

    @testset "apply ambiguous insertions" begin
        template = DNASeq("ACG")
        proposals = Proposal[Insertion(0, DNA_C),
                             Insertion(0, DNA_T)]
        @test_throws Rifraf.AmbiguousProposalsError apply_proposals(template, proposals)
    end

    @testset "apply ambiguous proposals" begin
        template = DNASeq("ACG")
        proposals = Proposal[Deletion(1),
                             Substitution(1, DNA_T)]
        @test_throws Rifraf.AmbiguousProposalsError apply_proposals(template, proposals)
    end
end
