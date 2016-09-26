using Quiver2.Proposals

using Base.Test

function test_apply_proposals()
    template = "ACG"
    proposals = Proposal[Proposals.Insertion(0, 'T'),
                         Proposals.Insertion(3, 'C'),
                         Proposals.Deletion(3),
                         Proposals.Substitution(2, 'T')]
    expected = "TATC"
    result = Proposals.apply_proposals(template, proposals)
    @test result == expected
end

function test_apply_codon_proposals()
    template = "ACG"
    proposals = Proposal[Proposals.CodonInsertion(0, ('T', 'C', 'C')),
                         Proposals.CodonDeletion(1)]
    expected = "TCC"
    result = Proposals.apply_proposals(template, proposals)
    @test result == expected
end

function test_apply_ambiguous_insertions()
    template = "ACG"
    proposals = Proposal[Proposals.CodonInsertion(0, ('T', 'C', 'C')),
                         Proposals.Insertion(0, 'T')]
    @test_throws AmbiguousProposalsError Proposals.apply_proposals(template, proposals)
end

function test_apply_ambiguous_proposals()
    template = "ACG"
    proposals = Proposal[Proposals.Deletion(1),
                         Proposals.Substitution(1, 'T')]
    @test_throws AmbiguousProposalsError Proposals.apply_proposals(template, proposals)
end

srand(1)

test_apply_proposals()
test_apply_codon_proposals()
test_apply_ambiguous_insertions()
test_apply_ambiguous_proposals()

