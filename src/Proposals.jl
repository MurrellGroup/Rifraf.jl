module Proposals

export Proposal, Substitution, Insertion,
       Deletion, CandProposal,
       are_unambiguous, base_shift, apply_proposals,
       AmbiguousProposalsError

abstract Proposal

immutable Substitution <: Proposal
    pos::Int
    base::Char
end

immutable Insertion <: Proposal
    pos::Int  # insert after this position
    base::Char
end

immutable Deletion <: Proposal
    pos::Int
end

function update_pos(p::Substitution, pos)
    return Substitution(pos, p.base)
end

function update_pos(p::Insertion, pos)
    return Insertion(pos, p.base)
end

function update_pos(p::Deletion, pos)
    return Deletion(pos)
end

immutable CandProposal
    proposal::Proposal
    score::Float64
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
        if t == Insertion
            push!(ins_positions, m.pos)
        else
            push!(other_positions, m.pos)
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
                         proposal::Deletion)
    return string(template[1:(proposal.pos - 1)],
                  template[(proposal.pos + 1):end])
end

base_shift_dict = Dict(Insertion => 1,
                       Deletion => -1,
                       Substitution => 0)
"""
How many bases to shift proposals after this one.
"""
function base_shift(m::Proposal)
    return base_shift_dict[typeof(m)]
end


immutable AmbiguousProposalsError <: Exception end


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
                remaining[i] = update_pos(m2, m2.pos + shift)
            end
        end
    end
    return template
end

end
