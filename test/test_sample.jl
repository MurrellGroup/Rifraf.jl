using Quiver2.Sample

using Base.Test

srand(1)

function test_do_substitutions()
    seq = random_seq(100)
    actual_error_p = zeros(100)
    actual_error_p[1] = 1.0
    sub_ratio = 1.0
    result = Quiver2.Sample.do_substitutions(seq, actual_error_p, sub_ratio)
    @test result[1] != seq[1]
    @test result[2:end] == seq[2:end]
end

function test_do_deletions()
    seq = random_seq(100)
    actual_error_p = zeros(100)
    actual_error_p[1] = 1.0
    del_ratio = 1.0
    result, new_error_p = Quiver2.Sample.do_deletions(seq, actual_error_p, del_ratio, false)
    @test length(result) == 99
    @test result == seq[2:end]
end

function test_do_codon_deletions()
    seq = random_seq(102)
    actual_error_p = zeros(102)
    actual_error_p[1] = 1.0
    del_ratio = 1.0
    result, new_error_p = Quiver2.Sample.do_deletions(seq, actual_error_p, del_ratio, true)
    @test length(result) == 99
    @test result == seq[4:end]
end

function test_do_insertions()
    seq = random_seq(100)
    actual_error_p = zeros(100)
    actual_error_p[1] = 1.0
    ins_ratio = 1.0
    result, new_error_p = Quiver2.Sample.do_insertions(seq, actual_error_p, ins_ratio, false)
    @test length(result) == 101
    @test result[2:end] == seq
end

function test_do_codon_insertions()
    seq = random_seq(102)
    actual_error_p = zeros(102)
    # must be 3.0 to ensure insertion, because codon insertion divides all probs by 3
    actual_error_p[1] = 3.0
    ins_ratio = 1.0
    result, new_error_p = Quiver2.Sample.do_insertions(seq, actual_error_p, ins_ratio, true)
    @test length(result) == 105
    @test result[4:end] == seq
end

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
    ref, template, template_error_p, seqs, actual, reported = sample(10, 99,
                                                                     0.05, (1.0, 1.0, 1.0),
                                                                     0.05, 0.01,
                                                                     0.5, 0.5,
                                                                     (1.0, 1.0, 1.0))
end

test_do_substitutions()
test_do_deletions()
test_do_codon_deletions()
test_do_insertions()
test_do_codon_insertions()
test_sample_from_reference()
test_sample_from_template()
test_sample()
