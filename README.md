# HalfDeep

## Dependencies

minimap2

samtools

genodsp (https://github.com/rsharris/genodsp)

Those must all be installed and in your `$PATH`.

## Installation

Copy the files from the repository:
```
git clone https://github.com/makovalab-psu/halfdeep
```

The .sh and .py scripts in the halfdeep directory must be in your `$PATH`. One
way to accomplish that is to add this to your shell startup file (e.g. .bashrc):
```
export PATH=path_to_halfdeep:${PATH}
```

## Expected directory layout

I've tried to set this all up to be similar to the VGP assembly pipeline, given
the fArcCen1 Assembly Tutorial and the pipeline scripts for minimap2.

Inputs are the assemblies and the fastq files (*.fasta.gz and *.fastq.gz in the
directory layout shown below). The assembly should have one of these file
extensions: .fa, .fasta, .fsa_nt, .fa.gz, .fasta.gz, or .fsa_nt.gz.

The final output is a halfdeep.dat file. This is a list of
<scaffold> <start> <end> intervals (origin 1, closed) that the process has
called as 'covered at half depth'.

```
.
└── genomic_data
   │ input.fofn
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
ls pacbio/*.fastq.gz > input.fofn
bam_depth.sh assembly_curated/mBalMus1.pri.cur.20190618.fasta.gz 1
bam_depth.sh assembly_curated/mBalMus1.pri.cur.20190618.fasta.gz 2
 ... (and so on, for 1..n where n is the number of fastq files)
halfdeep.sh assembly_curated/mBalMus1.pri.cur.20190618.fasta.gz
```

Bam_depth determines which fastq file by reading line i of input.fofn, where i
is the number given on the command line.

## Plotting

Open R with the working directory at, e.g.
genomic_data/halfdeep/mBalMus1.pri.cur.20190618.

Then, in R:
```
source("path_to_halfdeep/halfdeep.r")
scaffolds         = read_scaffold_lengths("scaffold_lengths.dat")
scaffoldToOffset  = linearized_scaffolds(scaffolds)
depth             = read_depth("depth.dat.gz",scaffoldToOffset)
halfDeep          = read_halfdeep("halfdeep.dat",scaffoldToOffset)
percentileToValue = read_percentiles("percentile_commands.sh")
```

To plot to the screen:
```
halfdeep_plot(scaffolds,depth,halfDeep,percentileToValue,assembly)
```

To plot to a file (currently only supports pdf):
```
halfdeep_plot(scaffolds,depth,halfDeep,percentileToValue,assembly,plotFilename="half_deep.dat.pdf")
```

To plot just a few scaffolds:
```
scaffoldsToPlot=c("Super_scaffold_8","Super_scaffold_17")
halfdeep_plot(scaffolds,depth,halfDeep,percentileToValue,assembly,scaffoldsToPlot=scaffoldsToPlot)
```


