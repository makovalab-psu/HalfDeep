### Brief Tutorial, Identifying a misassembled sex chromosome.

This directory contains a toy example for a small 1Mbp genome,
fake_genome.fa.gz, and simulated pacbio-like reads, fake_reads.fa.gz.

The fake genome consists of three scaffolds -- FAKE1, FAKE2, and FAKE3.

Read depth is roughly 15X. The ground truth of the reads is that FAKE2 and
FAKE3 are covered at full depth. FAKE1 is mostly covered at half depth, but the
interval (120137,167385) is full depth. This mimics a common assembly error for
heterogametic individuals, where the two sex chromosomes have been misassembled
into a single scaffold, with the PAR between them.

## Directory layout

Inputs are the assemblies and the reads files (*.fasta.gz in the directory
layout shown below).

```
.
└── genomic_data
   ├── assembly
   │   ├── fake_genome.fasta.gz
   ├── pacbio
   │   ├── fake_reads.fasta.gz
```

## (1) Index the assembly

Run minimap2 to create an index for the assembly. This must use the option
"-x map-pb" to match what is used later when we map reads to the assembly. 

```
cd genomic_data/assembly
minimap2 -x map-pb -d fake_genome.idx fake_genome.fasta.gz
```

We now have fake_genome.idx.

```
.
└── genomic_data
   ├── assembly
   │   ├── fake_genome.fasta.gz
   │   ├── fake_genome.idx
   ├── pacbio
   │   ├── fake_reads.fasta.gz
```

_TBD_

(be sure to include the expected plot).

