abstract Proposal

immutable Substitution <: Proposal
    pos::Int
    base::DNANucleotide
end

immutable Insertion <: Proposal
    pos::Int  # insert after this position
    base::DNANucleotide
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
function are_ambiguous(ms::Vector{Proposal})
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
    return !ins_good || !others_good
end

function apply_proposal(seq::DNASequence,
                         proposal::Substitution)
    return DNASequence(seq[1:(proposal.pos - 1)],
                       DNASequence([proposal.base]),
                       seq[(proposal.pos + 1):end])
end

function apply_proposal(seq::DNASequence,
                         proposal::Insertion)
    return DNASequence(seq[1:(proposal.pos)],
                       DNASequence([proposal.base]),
                       seq[(proposal.pos+1):end])
end

function apply_proposal(seq::DNASequence,
                         proposal::Deletion)
    return DNASequence(seq[1:(proposal.pos - 1)],
                       seq[(proposal.pos + 1):end])
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


function apply_proposals(seq::DNASequence,
                         proposals::Vector{Proposal})
    if are_ambiguous(proposals)
        throw(AmbiguousProposalsError())
    end
    remaining = deepcopy(proposals)
    while length(remaining) > 0
        m = pop!(remaining)
        seq = apply_proposal(seq, m)
        shift = base_shift(m)
        for i in 1:length(remaining)
            m2 = remaining[i]
            if m2.pos >= m.pos
                remaining[i] = update_pos(m2, m2.pos + shift)
            end
        end
    end
    return seq
end

function choose_candidates(candidates::Vector{CandProposal}, min_dist::Int)
    final_cands = CandProposal[]
    posns = Set{Int}()
    for c in sort(candidates, by=(c) -> c.score, rev=true)
        if any(Bool[(abs(c.proposal.pos - p) < min_dist) for p in posns])
            continue
        end
        push!(posns, c.proposal.pos)
        push!(final_cands, c)
    end
    return final_cands
end
