### Brief Tutorial, Identifying a misassembled sex chromosome.

This directory contains a toy example for a small 1Mbp genome,
fake_genome.fa.gz, and simulated pacbio-like reads, fake_reads.fa.gz.

The fake genome consists of three scaffolds -- FAKE1, FAKE2, and FAKE3.

Read depth is roughly 15X. The ground truth of the reads is that FAKE2 and
FAKE3 are covered at full depth. FAKE1 is mostly covered at half depth, but the
interval (120137,167385) is full depth. This mimics a common assembly error for
heterogametic individuals, where the two sex chromosomes have been misassembled
into a single scaffold, with the PAR between them.

## Expected directory layout

Inputs are the assemblies and the reads files (*.fasta.gz and *.fastq.gz in the
directory layout shown below).

```
.
└── genomic_data
   ├── pacbio
   │   ├── m54178_170623_204539.subreads.fastq.gz
   ├── assembly
   │   ├── mBalMus1.pri.cur.20190618.fasta.gz
```

_TBD_

(be sure to include the expected plot).

