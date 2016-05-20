module Mutations

export Mutation, Substitution, Insertion, CodonInsertion,
       Deletion, CodonDeletion, CandMutation,
       are_unambiguous, base_shift, affected_positions, apply_mutations,
       AmbiguousMutationsError, random_codon

abstract Mutation
abstract SingleMutation <: Mutation
abstract CodonMutation <: Mutation

type Substitution <: SingleMutation
    pos::Int
    base::Char
end

type Insertion <: SingleMutation
    pos::Int  # insert after this position
    base::Char
end

type CodonInsertion <: CodonMutation
    pos::Int  # insert after this position
    bases::Tuple{Char,Char,Char}
end

type Deletion <: SingleMutation
    pos::Int
end

type CodonDeletion <: CodonMutation
    pos::Int  # position of the first base to delete
end

type CandMutation
    mutation::Mutation
    score::Float64
end

function affected_positions(m::SingleMutation)
    return [m.pos]
end

function affected_positions(m::CodonInsertion)
    return [m.pos]
end

function affected_positions(m::CodonDeletion)
    return [m.pos, m.pos + 1, m.pos + 2]
end

"""
Whether mutations can be applied in any order.

- ensure only one substitution or deletion at every position
- ensure only one insertion at every position

"""
function are_unambiguous(ms::Vector{Mutation})
    ins_positions = []
    other_positions = []
    for i in 1:length(ms)
        m = ms[i]
        t = typeof(m)
        if t == Insertion || t == CodonInsertion
            append!(ins_positions, affected_positions(m))
        else
            append!(other_positions, affected_positions(m))
        end
    end
    ins_good = length(Set(ins_positions)) == length(ins_positions)
    others_good = length(Set(other_positions)) == length(other_positions)
    return ins_good && others_good
end

function update_template(template::AbstractString,
                         mutation::Substitution)
    return string(template[1:(mutation.pos - 1)],
                  mutation.base,
                  template[(mutation.pos + 1):end])
end

function update_template(template::AbstractString,
                         mutation::Insertion)
    return string(template[1:(mutation.pos)],
                  mutation.base,
                  template[(mutation.pos+1):end])
end

function update_template(template::AbstractString,
                         mutation::CodonInsertion)
    return string(template[1:(mutation.pos)],
                  join(mutation.bases),
                  template[(mutation.pos+1):end])
end

function update_template(template::AbstractString,
                         mutation::Deletion)
    return string(template[1:(mutation.pos - 1)],
                  template[(mutation.pos + 1):end])
end

function update_template(template::AbstractString,
                         mutation::CodonDeletion)
    return string(template[1:(mutation.pos - 1)],
                  template[(mutation.pos + 3):end])
end

base_shift_dict = Dict(Insertion => 1,
                       CodonInsertion => 3,
                       Deletion => -1,
                       CodonDeletion => -3,
                       Substitution => 0)

"""
How many bases to shift mutations after this one.
"""
function base_shift(m::Mutation)
    return base_shift_dict[typeof(m)]
end


type AmbiguousMutationsError <: Exception end


function apply_mutations(template::AbstractString,
                         mutations::Vector{Mutation})
    if !are_unambiguous(mutations)
        throw(AmbiguousMutationsError())
    end
    remaining = deepcopy(mutations)
    while length(remaining) > 0
        m = pop!(remaining)
        template = update_template(template, m)
        shift = base_shift(m)
        for i in 1:length(remaining)
            m2 = remaining[i]
            if m2.pos >= m.pos
                m2.pos += shift
            end
        end
    end
    return template
end

end
