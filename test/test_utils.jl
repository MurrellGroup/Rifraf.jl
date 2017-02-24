function inv_log10(logvals)
    log10(1.0 - exp10(logvals))
end


function check_all_cols(A, B, codon_moves)
    expected = A[end, end]
    if !(A[end, end] ≈ B[1, 1])
        return false
    end
    ncols = size(A)[2]
    # if codon_moves is true, we cannot expect every column to contain
    # the correct score
    # TODO: every three columns should
    if !codon_moves
        for j in 1:ncols
            Acol = sparsecol(A, j)
            Bcol = sparsecol(B, j)
            score = maximum(Acol + Bcol)
            if !(expected ≈ score)
                return false
            end
        end
    end
    return A[end, end] ≈ B[1, 1]
end
