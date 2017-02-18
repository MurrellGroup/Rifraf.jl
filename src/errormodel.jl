immutable ErrorModel
    mismatch::Float64
    insertion::Float64
    deletion::Float64
    codon_insertion::Float64
    codon_deletion::Float64
end

function ErrorModel(mismatch::Float64,
                    insertion::Float64,
                    deletion::Float64)
    return ErrorModel(mismatch, insertion, deletion, 0.0, 0.0)
end

function normalize(errors::ErrorModel)
    args = Float64[errors.mismatch,
                   errors.insertion,
                   errors.deletion,
                   errors.codon_insertion,
                   errors.codon_deletion]
    m, i, d, ci, cd = normalize(args)
    return ErrorModel(m, i, d, ci, cd)
end

immutable Scores
    mismatch::Float64
    insertion::Float64
    deletion::Float64
    codon_insertion::Float64
    codon_deletion::Float64
end

function Scores(errors::ErrorModel;
                mismatch::Float64=0.0,
                insertion::Float64=0.0,
                deletion::Float64=0.0)
    args = Float64[errors.mismatch,
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
