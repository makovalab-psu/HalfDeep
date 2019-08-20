# HalfDeep

## Dendencies

minimap2
samtools
https://github.com/rsharris/genodsp

Those should all be installed an in your `$PATH`. Additionally, the scripts in
this directory should be in your `$PATH`.


## Expected directory layout

Inputs are the assemblies and the fastq files.

The final output is a halfdeep.dat file. This is a list of
<scaffold> <start> <end> intervals (origin 1, closed) that the process has\
called as 'covered at half depth'.

```
.
└── genomic_data
   ├── pacbio
   │   ├── m54178_170623_204539.subreads.fastq.gz
   │   ├── m54178_170624_064412.subreads.fastq.gz
   │   ├──  ...
   ├── assembly_curated
   │   ├── mBalMus1.pri.cur.20190618.fasta.gz
   │   ├── mBalMus1.pri.cur.20190618.idx
   │   ├── mBalMus1.alt.cur.20190618.fasta.gz
   │   ├── mBalMus1.alt.cur.20190618.idx
   ├── halfdeep
   │   ├── mBalMus1.pri.cur.20190618
   │   |   ├── mapped_reads
   │   |   │   ├── m54178_170623_204539.subreads.bam
   │   |   │   ├── m54178_170623_204539.subreads.sort.bam
   │   |   │   ├── m54178_170623_204539.subreads.sort.bam.bai
   │   │   |   ├── m54178_170623_204539.subreads.depth.dat.gz
   │   |   │   ├── m54178_170624_064412.subreads.bam
   │   |   │   ├── m54178_170624_064412.subreads.sort.bam
   │   |   │   ├── m54178_170624_064412.subreads.sort.bam.bai
   │   │   |   ├── m54178_170624_064412.subreads.depth.dat.gz
   │   |   │   ├──  ...
   │   │   ├── scaffold_lengths.dat
   │   │   ├── depth.dat.gz
   │   │   ├── percentile_commands.sh
   │   │   ├── halfdeep.dat
   │   ├── mBalMus1.alt.cur.20190618
   │   │   ├──  ...
```

## Running

```
cd genomic_data
bam_depth.sh assembly_curated/mBalMus1.pri.cur.20190618.fasta.gz 1
bam_depth.sh assembly_curated/mBalMus1.pri.cur.20190618.fasta.gz 2
 ... (and so on, for 1..n where n is the number of fastq files)
halfdeep.sh assembly_curated/mBalMus1.pri.cur.20190618.fasta.gz
```
