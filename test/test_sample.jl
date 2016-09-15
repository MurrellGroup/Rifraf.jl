using Quiver2.Sample

using Base.Test

srand(1)

function test_sample_from_reference()
    reference = random_seq(102)
    template = sample_from_reference(reference,
                                     0.1,
                                     (1.0, 1.0, 1.0))
end

function test_sample_from_template()
    template = random_seq(102)
    template_error_p = 0.05 * ones(102)
    seq, actual, reported = sample_from_template(template,
                                                 template_error_p,
                                                 (1.0, 1.0, 1.0),
                                                 0.5,
                                                 0.5)
end

function test_sample()
    (ref, template, template_error_p,
     seqs, actual, reported) = sample(10, 99,
                                      0.05, (1.0, 1.0, 1.0),
                                      0.05, 0.01,
                                      0.5, 0.5,
                                      (1.0, 1.0, 1.0))
end

test_sample_from_reference()
test_sample_from_template()
test_sample()
