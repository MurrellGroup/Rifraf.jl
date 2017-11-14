# API

The main functions and types in Rifraf.jl.

## Functions

```@docs
rifraf
```

### Sequence simulations

```@docs
sample_sequences
write_samples
read_samples
```

### Utility IO functions

Rifraf.jl provides some utility functions for reading and writing
FASTQ and FASTA files. This functionality uses BioSequences.jl.


```@docs
read_fastq_records
read_fastq
write_fastq
read_fasta_records
read_fasta
write_fasta
```


## Types

```@docs
RifrafParams
RifrafResult
ErrorModel
Scores
```
