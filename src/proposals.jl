abstract type Proposal end

struct Substitution <: Proposal
    pos::Int
    base::DNA
end

struct Insertion <: Proposal
    pos::Int  # insert after this position
    base::DNA
end

struct Deletion <: Proposal
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

struct ScoredProposal
    proposal::Proposal
    score::Score
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

function get_proposal_base(seq::DNASeq, proposal::Substitution,
                           last_del_pos::Int)
    return DNASeq([proposal.base])
end

function get_proposal_base(seq::DNASeq, proposal::Insertion,
                           last_del_pos::Int)
    if proposal.pos > 0 && last_del_pos != proposal.pos
        return DNASeq([seq[proposal.pos], proposal.base])
    end
    return DNASeq([proposal.base])
end

function get_proposal_base(seq::DNASeq, proposal::Deletion,
                           last_del_pos)
    return DNASeq([])
end


struct AmbiguousProposalsError <: Exception end


function apply_proposals(seq::DNASeq,
                         proposals::Vector{Proposal})
    if are_ambiguous(proposals)
        throw(AmbiguousProposalsError())
    end
    result = DNASeq[]
    next = 1
    # keep track of deletions, so an insertion following a deletion
    # knows not to include the deleted base
    last_del_pos = 0
    # sort by position, making sure deletions go before insertions
    proposals = sort(proposals, by = p -> (p.pos, typeof(p) == Deletion ? 0 : 1))
    for p in proposals
        push!(result, seq[next:p.pos - 1])
        push!(result, get_proposal_base(seq, p, last_del_pos))
        next = p.pos + 1
        if typeof(p) == Deletion
            last_del_pos = p.pos
        end
    end
    push!(result, seq[next:end])
    return DNASeq(result...)
end

function choose_candidates(candidates::Vector{ScoredProposal}, min_dist::Int)
    final_cands = ScoredProposal[]
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
