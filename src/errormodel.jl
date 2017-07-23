struct ErrorModel
    mismatch::Real
    insertion::Real
    deletion::Real
    codon_insertion::Real
    codon_deletion::Real
end

function ErrorModel(mismatch, insertion, deletion)
    return ErrorModel(mismatch, insertion, deletion, 0, 0)
end

function normalize(errors::ErrorModel)
    args = [errors.mismatch,
            errors.insertion,
            errors.deletion,
            errors.codon_insertion,
            errors.codon_deletion]
    m, i, d, ci, cd = normalize(args)
    return ErrorModel(m, i, d, ci, cd)
end

struct Scores
    mismatch::Score
    insertion::Score
    deletion::Score
    codon_insertion::Score
    codon_deletion::Score
end

function Scores(errors::ErrorModel;
                mismatch=0.0,
                insertion=0.0,
                deletion=0.0)
    args = [errors.mismatch,
            errors.insertion,
            errors.deletion,
            errors.codon_insertion,
            errors.codon_deletion]
    m, i, d, ci, cd = log10.(normalize(args))
    return Scores(Score(m + mismatch),
                  Score(i + insertion),
                  Score(d + deletion),
                  Score(ci + 3 * insertion),
                  Score(cd + 3 * deletion))
end
