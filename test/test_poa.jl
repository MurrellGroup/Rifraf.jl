using Quiver2.Model
using Quiver2.Sample
using Quiver2.PartialOrderAligner

using Base.Test

srand(1)

function test_poa()
    n_seqs=10
    ref_len=27
    ref_error_rate = 0.03
    template_error_mean = 0.01
    template_error_std = 0.001
    log_seq_actual_std = 0.01
    log_seq_reported_std = 0.01

    ref_errors = ErrorModel(1.0, 0.0, 0.0, 1.0, 1.0)
    errors = ErrorModel(3.0, 1.0, 1.0, 0.0, 0.0)
    scores = Quiver2.modified_emissions(Scores(errors))

    (ref, template, template_error_p,
     seqs, actual, phreds) = sample(n_seqs, ref_len,
                                    ref_error_rate,
                                    ref_errors,
                                    template_error_mean,
                                    template_error_std,
                                    log_seq_actual_std,
                                    log_seq_reported_std,
                                    errors)
    seqstrings = String[convert(String, s) for s in seqs]
    g = PartialOrderGraph(seqstrings, scores)
    c = consensus(g)
    @test c == convert(String, template)
end

test_poa()
