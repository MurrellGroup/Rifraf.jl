using Quiver2.Mutations

using Base.Test

function test_apply_mutations()
    template = "ACG"
    mutations = Mutation[Mutations.Insertion(0, 'T'),
                         Mutations.Insertion(3, 'C'),
                         Mutations.Deletion(3),
                         Mutations.Substitution(2, 'T')]
    expected = "TATC"
    result = Mutations.apply_mutations(template, mutations)
    @test result == expected
end

function test_apply_codon_mutations()
    template = "ACG"
    mutations = Mutation[Mutations.CodonInsertion(0, ('T', 'C', 'C')),
                         Mutations.CodonDeletion(1)]
    expected = "TCC"
    result = Mutations.apply_mutations(template, mutations)
    @test result == expected
end

function test_apply_ambiguous_insertions()
    template = "ACG"
    mutations = Mutation[Mutations.CodonInsertion(0, ('T', 'C', 'C')),
                         Mutations.Insertion(0, 'T')]
    @test_throws AmbiguousMutationsError Mutations.apply_mutations(template, mutations)
end

function test_apply_ambiguous_mutations()
    template = "ACG"
    mutations = Mutation[Mutations.Deletion(1),
                         Mutations.Substitution(1, 'T')]
    @test_throws AmbiguousMutationsError Mutations.apply_mutations(template, mutations)
end

srand(1)

test_apply_mutations()
test_apply_codon_mutations()
test_apply_ambiguous_mutations()

