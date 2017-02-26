function inv_log10(logvals)
    log10(1.0 - exp10(logvals))
end


function check_all_cols(A, B, codon_moves)
    @test A[end, end] ≈ B[1, 1]
    expected = A[end, end]
    ncols = size(A)[2]
    # if codon_moves is true, we cannot expect every column to contain
    # the correct score. every three columns should, though.
    if codon_moves
        for j in 1:ncols-2
            score = maximum(A[:, j:j+2] + B[:, j:j+2])
            @test score ≈ expected
        end
    else
        for j in 1:ncols
            score = maximum(A[:, j] + B[:, j])
            @test score ≈ expected
        end
    end
end
