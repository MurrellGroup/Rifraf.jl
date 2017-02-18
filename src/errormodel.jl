immutable ErrorModel
    mismatch::ErrorProb
    insertion::ErrorProb
    deletion::ErrorProb
    codon_insertion::ErrorProb
    codon_deletion::ErrorProb
end

function ErrorModel(mismatch::ErrorProb,
                    insertion::ErrorProb,
                    deletion::ErrorProb)
    return ErrorModel(mismatch, insertion, deletion, 0.0, 0.0)
end

function normalize(errors::ErrorModel)
    args = ErrorProb[errors.mismatch,
                   errors.insertion,
                   errors.deletion,
                   errors.codon_insertion,
                   errors.codon_deletion]
    m, i, d, ci, cd = normalize(args)
    return ErrorModel(m, i, d, ci, cd)
end


immutable Scores
    mismatch::Score
    insertion::Score
    deletion::Score
    codon_insertion::Score
    codon_deletion::Score
end

function Scores(errors::ErrorModel;
                mismatch::Score=0.0,
                insertion::Score=0.0,
                deletion::Score=0.0)
    args = Score[errors.mismatch,
                 errors.insertion,
                 errors.deletion,
                 errors.codon_insertion,
                 errors.codon_deletion]
    m, i, d, ci, cd = log10(normalize(args))
    return Scores(m + mismatch,
                  i + insertion,
                  d + deletion,
                  ci + 3 * insertion,
                  cd + 3 * deletion)
end
