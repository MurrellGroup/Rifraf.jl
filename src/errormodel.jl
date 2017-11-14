"""
    ErrorModel(mismatch, insertion, deletion, codon_insertion, codon_deletion)

Error model for sequencing.

Each field contains the relative rate of of that kind of error. For
instance, this model breaks the error rate into 80% mismatches, 10%
codon insertions, and 10% codon deletions:
`ErrorModel(8, 0, 0, 1, 1)`.

# Fields:
- `mismatch::Real`
- `insertion::Real`
- `deletion::Real`
- `codon_insertion::Real`
- `codon_deletion::Real`

"""
struct ErrorModel
    mismatch::Real
    insertion::Real
    deletion::Real
    codon_insertion::Real
    codon_deletion::Real
end

"""error model without codon errors"""
function ErrorModel(mismatch, insertion, deletion)
    return ErrorModel(mismatch, insertion, deletion, 0, 0)
end

"""turn error rates into probabilities"""
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

"""
    Scores(errors; mismatch, insertion, deletion)

Derive alignment scores from an error model.

Takes extra penalties to add to the mismatch, insertion, and deletion
scores.

# Arguments:
- `errors::ErrorModel`:
- `mismatch::Real`: substitution
- `insertion::Real`: insertion
- `deletion::Real`: deletion

"""
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
