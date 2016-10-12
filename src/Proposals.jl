module Proposals

export Proposal, Substitution, Insertion, CodonInsertion,
       Deletion, CodonDeletion, CandProposal,
       are_unambiguous, base_shift, affected_positions, apply_proposals,
       AmbiguousProposalsError

abstract Proposal
abstract SingleProposal <: Proposal
abstract CodonProposal <: Proposal

type Substitution <: SingleProposal
    pos::Int
    base::Char
end

type Insertion <: SingleProposal
    pos::Int  # insert after this position
    base::Char
end

type CodonInsertion <: CodonProposal
    pos::Int  # insert after this position
    bases::Tuple{Char,Char,Char}
end

type Deletion <: SingleProposal
    pos::Int
end

type CodonDeletion <: CodonProposal
    pos::Int  # position of the first base to delete
end

type CandProposal
    proposal::Proposal
    score::Float64
end

function affected_positions(m::SingleProposal)
    return [m.pos]
end

function affected_positions(m::CodonInsertion)
    return [m.pos]
end

function affected_positions(m::CodonDeletion)
    return [m.pos, m.pos + 1, m.pos + 2]
end

"""
Whether proposals can be applied in any order.

- ensure only one substitution or deletion at every position
- ensure only one insertion at every position

"""
function are_unambiguous(ms::Vector{Proposal})
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

function update_template(template::String,
                         proposal::Substitution)
    return string(template[1:(proposal.pos - 1)],
                  proposal.base,
                  template[(proposal.pos + 1):end])
end

function update_template(template::String,
                         proposal::Insertion)
    return string(template[1:(proposal.pos)],
                  proposal.base,
                  template[(proposal.pos+1):end])
end

function update_template(template::String,
                         proposal::CodonInsertion)
    return string(template[1:(proposal.pos)],
                  join(proposal.bases),
                  template[(proposal.pos+1):end])
end

function update_template(template::String,
                         proposal::Deletion)
    return string(template[1:(proposal.pos - 1)],
                  template[(proposal.pos + 1):end])
end

function update_template(template::String,
                         proposal::CodonDeletion)
    return string(template[1:(proposal.pos - 1)],
                  template[(proposal.pos + 3):end])
end

base_shift_dict = Dict(Insertion => 1,
                       CodonInsertion => 3,
                       Deletion => -1,
                       CodonDeletion => -3,
                       Substitution => 0)

"""
How many bases to shift proposals after this one.
"""
function base_shift(m::Proposal)
    return base_shift_dict[typeof(m)]
end


type AmbiguousProposalsError <: Exception end


function apply_proposals(template::String,
                         proposals::Vector{Proposal})
    if !are_unambiguous(proposals)
        throw(AmbiguousProposalsError())
    end
    remaining = deepcopy(proposals)
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
