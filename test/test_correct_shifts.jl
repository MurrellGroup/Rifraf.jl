using BioSymbols
using BioSequences
using Base.Test

using Rifraf

@testset "correct shifts" begin
    @testset "correct shifts 1" begin
        consensus = DNASeq("TTTT")
        reference = DNASeq("TTT")
        expected = DNASeq("TTT")
        result = Rifraf.correct_shifts(consensus, reference)
        @test result == expected
    end
    @testset "correct shifts 2" begin
        consensus = DNASeq("TT")
        reference = DNASeq("TTT")
        expected = DNASeq("TTT")
        result = Rifraf.correct_shifts(consensus, reference)
        @test result == expected
    end
    @testset "correct shifts 3" begin
        consensus = DNASeq("TTTACCC")
        reference = DNASeq("TTTCGC")
        expected = DNASeq("TTTCCC")
        result = Rifraf.correct_shifts(consensus, reference)
        @test result == expected
    end
    @testset "correct shifts 4" begin
        consensus = DNASeq("TTTAAACCC")
        reference = DNASeq("TTTCGC")
        expected = DNASeq("TTTAAACCC")
        result = Rifraf.correct_shifts(consensus, reference)
        @test result == expected
    end

    # TODO: this test fails because it currently prefers one insertion
    # to two deletions.

    # @testset "correct shifts 5" begin
    #     consensus = DNASeq("TTTAACCC")
    #     reference = DNASeq("TTTCCC")
    #     expected = DNASeq("TTTCCC")
    #     result = Rifraf.correct_shifts(consensus, reference)
    #     @test result == expected
    # end
end
