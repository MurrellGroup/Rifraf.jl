using Quiver2.Sample

using Base.Test

srand(1)

# test that codon sampling is always a multiple of 3
begin
    for i in 1:100
        template = random_seq(3 * 100)
        error_rate = 0.1
        sub_ratio = 1
        ins_ratio = 3
        del_ratio = 3
        error_std = 0.001
        seq, logp = sample_from_template(template, error_rate,
                                         sub_ratio, ins_ratio, del_ratio,
                                         error_std, codon=true)
        @test length(seq) % 3 == 0
    end
end
