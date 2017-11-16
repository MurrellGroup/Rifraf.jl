# Examples

Some examples showing how to use RIFRAF.

## On simulated data

First, we generate a random 1,200 bp template, along with a reference
and twenty simulated reads. We run RIFRAF both without and with the
reference and compare the result to the expected template.

```@repl
using Rifraf

sampled = sample_sequences(20, 1200);
(reference, template, _, sequences, _, phreds, _, _) = sampled;

result = rifraf(sequences, phreds;
                params=RifrafParams(batch_fixed_size=3, batch_size=5,
                                    verbose=1, max_iters=20));
result.consensus == template

result = rifraf(sequences, phreds; reference=reference,
                params=RifrafParams(verbose=1, max_iters=20));

result.consensus == template
```

## Reading data from FASTQ files

Rifraf.jl also provides a set of convenience functions for reading and
writing FASTA and FASTQ files. For instance, here is one way to run
RIFRAF on sequences from a file:

```
sequences, phreds, names = Rifraf.read_fastq("/path/to/sequences.fastq")
reference = Rifraf.read_fasta("/path/to/reference.fasta")[1]
result = rifraf(sequences, phreds; reference=reference)
```

There are also convenience functions for reading and writing the
output of `sample_sequences`.

```
reference, template, t_error, sequences, _, phreds, _, _ = sample_sequences()
Rifraf.write_samples("/path/to/basename", reference, template, t_error, sequences, phreds)
reference, template, t_error, sequences, phreds = Rifraf.read_samples("/path/to/basename")
```

## Command-line script

`scripts/rifraf.jl` is a command-line script for processing many sets
of reads at once. Julia takes some time to start up, so this script is
only recommended for long sequences or large numbers of reads.

The `data` directory includes some example data for testing this
script. The following command finds a consensus for each FASTQ file
that matches the glob `input-reads-*.fastq` and writes them to
`results.fasta`.

```
julia ./scripts/rifraf.jl \
    --reference ./data/references.fasta \
    --reference-map ./data/ref-map.tsv \
    --phred-cap 30 \
    --ref-errors 8,0.1,0.1,1,1 \
    1,2,2 \
    "./data/input-reads-*.fastq" \
    ./data/consensus-results.fasta
```

## Allowing frameshifts during frame correction

RIFRAF's default parameters penalize frameshift-causing indels
extremely heavily. If the template really does contain a frame shift
mutation, it will likely be removed from the result. To detect real
frameshifts, with a small risk of allowing some spurious ones, the
penalties for single insertions or deletions must be tuned. For
example:


```@repl
using Rifraf

sampled = Rifraf.sample_sequences(5, 3001; error_rate=.005);
(reference, template, _, sequences, _, phreds, _, _) = sampled;

result = rifraf(sequences, phreds; reference=reference);
length(result.consensus) % 3 == 0

result = rifraf(sequences, phreds; reference=reference,
                params=RifrafParams(ref_scores=Scores(ErrorModel(10, 1, 1, 1, 1)),
                                    ref_indel_mult=1.2, max_ref_indel_mults=3));
length(result.consensus) % 3 == 0

```