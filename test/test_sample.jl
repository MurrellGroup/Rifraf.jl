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
    seq, actual, reported = sample_from_template(template,
                                                 template_error_p,
                                                 seq_errors,
                                                 0.5, 0.5)
end

function test_sample()
    len = 900
    error_rate = 0.01
    (ref, template, template_error_p,
     seqs, actual, phreds) = sample(3, len,
                                    0.05, ref_errors,
                                    error_rate, 1e-1,
                                    0.1, 0.01,
                                    seq_errors)
    @test_approx_eq_eps mean(template_error_p) error_rate 0.1
end

function test_sample_mixture()
    error_rate = 0.01
    (ref, template, template_error_p,
     seqs, actual, phreds) = sample_mixture((10, 10), 99, 1,
                                            0.05, ref_errors,
                                            error_rate, 1e-1,
                                            0.1, 0.01,
                                            seq_errors)
    @test length(seqs) == 20
end

test_sample_reference()
test_sample_from_template()
test_sample()
test_sample_mixture()
