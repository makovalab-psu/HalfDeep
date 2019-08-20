# HalfDeep

## Dendencies

https://github.com/rsharris/genodsp

## Expected directory layout

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


