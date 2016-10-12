module PartialOrderAligner

using LightGraphs
using Quiver2.Model

import Base.length

export PartialOrderGraph, length, align, add_alignment!, add_sequence!, consensus

immutable NodeLabel
    # which sequence this node came from
    seq::Int
    # position in that sequence
    position::Int
end


immutable PartialOrderGraph
    graph::DiGraph
    # has all same nodes as `graph`. edges link aligned nodes.
    aligned_nodes::Graph
    node_chars::Vector{Char}
    node_labels::Vector{Vector{NodeLabel}}
    n_seqs::Vector{Int}  # needs to be a vector so we can change its contents
    scores::Scores
end


function PartialOrderGraph(scores::Scores)
    if scores.codon_insertion > -Inf || scores.codon_deletion > -Inf
        error("codon indels not allowed")
    end
    g = DiGraph()
    aligned_nodes = Graph()
    node_chars = Char[]
    node_labels = Vector{NodeLabel}[]
    return PartialOrderGraph(g, aligned_nodes,
                             node_chars, node_labels, Int[0], scores)
end


function PartialOrderGraph(sequences::Vector{String},
                           scores::Scores)
    g = PartialOrderGraph(scores)
    for s in sequences
        add_sequence!(g, s)
    end
    return g
end


immutable GraphAlignment
    sequence::String
    # 0 for insertion
    nodes::Vector{Int}

    function GraphAlignment(s::String, n::Vector{Int})
        if length(s) != length(n)
            error("sequence and node vector lengths differ")
        end
        if minimum(n) < 0
            error("nodes must be non-negative")
        end
        return new(s, n)
    end
end


function length(a::GraphAlignment)
    return length(a.sequence)
end


function align(g::PartialOrderGraph, sequence::String)
    if nv(g.graph) == 0
        # trivial alignment
        return GraphAlignment(sequence, zeros(Int, length(sequence)))
    end
    sorted_nodes = topological_sort_by_dfs(g.graph)

    node_sort_posn = zeros(Int, length(sorted_nodes))
    for i in 1:length(sorted_nodes)
        node_sort_posn[sorted_nodes[i]] = i
    end

    # initialize array
    u = length(sorted_nodes) + 1
    v = length(sequence) + 1
    result = zeros(Float64, (u, v))
    moves = Array(Tuple{Int, Int}, (u, v))
    moves[1, 1] = (0, 0)
    for i = 2:u
        result[i, 1] = result[i-1, 1] + g.scores.deletion
        moves[i, 1] = (i-1, 1)
    end
    for j = 2:v
        result[1, j] = result[1, j-1] + g.scores.insertion
        moves[1, j] = (1, j-1)
    end

    # fill in matrices; column-major order
    for j = 2:v
        for i = 2:u
            node = sorted_nodes[i - 1]
            best_score = typemin(Float64)
            move = (0, 0)

            # implied start node at first row
            prev_row = 1
            # align
            is_match = false
            score = result[prev_row, j-1] + (is_match ? 0.0 : g.scores.mismatch)
            if score > best_score
                best_score = score
                move = (prev_row, j-1)
            end
            # skip node (vertical move)
            score = result[prev_row, j] + g.scores.deletion
            if score > best_score
                best_score = score
                move = (prev_row, j)
            end

            # consider all other predecessor nodes
            for prev in in_neighbors(g.graph, node)
                prev_row = node_sort_posn[prev] + 1
                if prev_row >= i
                    error("prev_row >= i")
                end
                # align
                is_match = (g.node_chars[prev] == sequence[j-1])
                score = result[prev_row, j-1] + (is_match ? 0.0 : g.scores.mismatch)
                if score > best_score
                    best_score = score
                    move = (prev_row, j-1)
                end
                # skip node (vertical move)
                score = result[prev_row, j] + g.scores.deletion
                if score > best_score
                    best_score = score
                    move = (prev_row, j)
                end
            end
            # insert new node (horizontal move)
            score = result[i, j-1] + g.scores.insertion
            if score > best_score
                best_score = score
                move = (i, j-1)
            end
            result[i, j] = best_score
            moves[i, j] = move
        end
    end

    # backtrace
    aligned_nodes = Int[]
    i = u
    j = v
    while i > 1 || j > 1
        prev_i, prev_j = moves[i, j]
        # horizontal move: insert new node
        if prev_i == i && prev_j != j
            push!(aligned_nodes, 0)
        # match: align to node
        elseif prev_i < i && prev_j < j
            push!(aligned_nodes, sorted_nodes[i-1])
        end
        # ignore vertical moves
        i = prev_i
        j = prev_j
    end
    aligned_nodes = reverse(aligned_nodes)
    return GraphAlignment(sequence, aligned_nodes)
end


function add_alignment!(g::PartialOrderGraph, alignment::GraphAlignment)
    seq_id = g.n_seqs[1] + 1

    prev_node = 0
    for i in 1:length(alignment)
        aligned_node = alignment.nodes[i]
        char = alignment.sequence[i]
        if has_vertex(g.graph, aligned_node) && g.node_chars[aligned_node] == char
            # add to node label
            push!(g.node_labels[aligned_node], NodeLabel(seq_id, i))
            prev_node = aligned_node
        else
            # create new node; maybe add alignment edge
            add_vertex!(g.graph)
            new_node = nv(g.graph)
            push!(g.node_chars, char)
            if has_vertex(g.graph, prev_node)
                add_edge!(g.graph, prev_node, new_node)
            end
            push!(g.node_labels, NodeLabel[NodeLabel(seq_id, i)])
            add_vertex!(g.aligned_nodes)
            if aligned_node != 0
                # TODO: check if any aligned nodes have same char
                add_edge!(g.aligned_nodes, aligned_node, new_node)
            end
            prev_node = new_node
        end
    end
    g.n_seqs[1] += 1
end


function add_sequence!(g::PartialOrderGraph, sequence::String)
    alignment = align(g, sequence)
    add_alignment!(g, alignment)
end


function consensus(g::PartialOrderGraph)
    n = nv(g.graph)
    sorted_nodes = topological_sort_by_dfs(g.graph)
    prev_nodes = zeros(Int, n)
    node_scores = zeros(Float64, n)

    # compute edge weights
    edge_weights = zeros(Int, (n, n))
    for e in edges(g.graph)
        start, stop = e
        start_labels = g.node_labels[start]
        stop_labels = g.node_labels[stop]
        s1 = Set(label.seq for label in start_labels)
        s2 = Set(label.seq for label in stop_labels)
        edge_weights[start, stop] = length(intersect(s1, s2))
    end

    # compute node scores
    for i = 1:n
        node = sorted_nodes[i]
        # set best predecessor
        best_score = typemin(Float64)
        best_prev = 0

        for prev_node in in_neighbors(g.graph, node)
            score = node_scores[prev_node] + edge_weights[prev_node, node]
            if score > best_score
                best_score = score
                best_prev = prev_node
            end
        end
        if best_prev > 0
            node_scores[node] = best_score
            prev_nodes[node] = best_prev
        end
    end

    # backtrace
    chars = Char[]
    score, node = findmax(node_scores)
    while node > 0
        push!(chars, g.node_chars[node])
        node = prev_nodes[node]
    end
    return string(reverse(chars)...)
end

end
