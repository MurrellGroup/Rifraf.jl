function inv_log10(logvals)
    log10(1.0 - exp10(logvals))
end


function check_all_cols(A, B, codon_moves)
    if !(A[end, end] ≈ B[1, 1])
        println()
        display(A)
        println()
        display(B)
        println()
    end
    @test A[end, end] ≈ B[1, 1]
    expected = A[end, end]
    ncols = size(A)[2]
    # if codon_moves is true, we cannot expect every column to contain
    # the correct score. every three columns should, though.
    if codon_moves
        for j in 1:ncols-2
            Acol = sparsecol(A, j)
            Bcol = sparsecol(B, j)
            score = maximum(Acol + Bcol)

            Acol = sparsecol(A, j+1)
            Bcol = sparsecol(B, j+1)
            score = max(score, maximum(Acol + Bcol))

            Acol = sparsecol(A, j+2)
            Bcol = sparsecol(B, j+2)
            score = max(score, maximum(Acol + Bcol))

            @test expected ≈ score
        end
    else
        for j in 1:ncols
            Acol = sparsecol(A, j)
            Bcol = sparsecol(B, j)
            score = maximum(Acol + Bcol)
            @test expected ≈ score
        end
    end
end
