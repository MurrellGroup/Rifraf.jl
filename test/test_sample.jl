using Quiver2.Sample
using Quiver2.Model

using Base.Test

srand(1)

const seq_errors = ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0)
const ref_errors = ErrorModel(1.0, 0.0, 0.0, 1.0, 1.0)

function test_sample_reference()
    template = random_seq(102)
    reference = sample_reference(template, 0.1, ref_errors)
end

function test_sample_from_template()
    template = random_seq(102)
    template_error_p = 0.05 * ones(102)
    phred_scale = 5.0
    actual_std = 3.0
    reported_std = 0.3
    seq, actual, reported = sample_from_template(template,
                                                 template_error_p,
                                                 seq_errors,
                                                 phred_scale,
                                                 actual_std,
                                                 reported_std)
end

function test_sample()
    nseqs = 3
    len = 900
    ref_error_rate = 0.05
    error_rate = 0.01
    alpha = 1.0
    phred_scale = 5.0
    actual_std = 3.0
    reported_std = 1.0
    (ref, template, template_error_p,
     seqs, actual, phreds) = sample(nseqs, len,
                                    ref_error_rate, ref_errors,
                                    error_rate, alpha,
                                    phred_scale,
                                    actual_std, reported_std,
                                    seq_errors)
    @test_approx_eq_eps mean(template_error_p) error_rate 0.1
end

function test_sample_mixture()
    nseqs = (3, 3)
    len = 900
    n_diffs = 3
    ref_error_rate = 0.05
    error_rate = 0.01
    alpha = 1.0
    phred_scale = 5.0
    actual_std = 3.0
    reported_std = 1.0
    (ref, template, template_error_p,
     seqs, actual, phreds) = sample_mixture(nseqs, len, n_diffs,
                                            ref_error_rate, ref_errors,
                                            error_rate, alpha,
                                            phred_scale,
                                            actual_std, reported_std,
                                            seq_errors)
    @test length(seqs) == sum(nseqs)
end

test_sample_reference()
test_sample_from_template()
test_sample()
test_sample_mixture()
