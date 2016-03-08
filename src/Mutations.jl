__precompile__()
module Mutations

export Mutation, Substitution, Insertion, Deletion, CandMutation

abstract Mutation

immutable Substitution <: Mutation
    pos::Int
    base::Char
end

immutable Insertion <: Mutation
    pos::Int
    base::Char
end

immutable Deletion <: Mutation
    pos::Int
end

immutable CandMutation
    mutation::Mutation
    score::Float64
end

end
