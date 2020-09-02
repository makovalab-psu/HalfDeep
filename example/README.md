### Brief Tutorial -- Identifying a misassembled sex chromosome.

This directory contains a toy example for a small 1Mbp genome,
fake_genome.fa.gz, and simulated pacbio-like reads, fake_reads_*.fa.gz.

The fake genome consists of three scaffolds -- FAKE1, FAKE2, and FAKE3.

Read depth is roughly 15X. The ground truth of the reads is that FAKE2 and
FAKE3 are covered at full depth. FAKE1 is mostly covered at half depth, but the
interval (120137,167385) is full depth. This mimics a common assembly error for
heterogametic individuals, where the two sex chromosomes have been misassembled
into a single scaffold, with the PAR between them.

## Initial directory layout

Inputs are an assembly and five reads files (*.fasta.gz in the directory
layout shown below).

```
.
└── genomic_data
   ├── assembly
   │   ├── fake_genome.fasta.gz
   ├── pacbio
   │   ├── fake_reads_001.fasta.gz
   │   ├── fake_reads_002.fasta.gz
   │   ├── fake_reads_003.fasta.gz
   │   ├── fake_reads_004.fasta.gz
   │   ├── fake_reads_005.fasta.gz
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
   │   ├── >>> fake_genome.idx <<<
   ├── pacbio
   │   ├── fake_reads_001.fasta.gz
   │   ├── fake_reads_002.fasta.gz
   │   ├── fake_reads_003.fasta.gz
   │   ├── fake_reads_004.fasta.gz
   │   ├── fake_reads_005.fasta.gz
```

## (2) Create the 'file of filenames'.

The pipeline needs a file containing a list of all the reads filenames. We
could create this manually with a text editor. But in this case we can use the
ls command to create it.

```
cd genomic_data
ls pacbio/*.fastq.gz > input.fofn
```

We now have input.fofn.

```
.
└── genomic_data
   │ >>> input.fofn <<<
   ├── assembly
   │   ├── fake_genome.fasta.gz
   │   ├── fake_genome.idx
   ├── pacbio
   │   ├── fake_reads_001.fasta.gz
   │   ├── fake_reads_002.fasta.gz
   │   ├── fake_reads_003.fasta.gz
   │   ├── fake_reads_004.fasta.gz
   │   ├── fake_reads_005.fasta.gz
```

This is what input.fofn should look like:

```
pacbio/fake_reads_001.fasta.gz
pacbio/fake_reads_002.fasta.gz
pacbio/fake_reads_003.fasta.gz
pacbio/fake_reads_004.fasta.gz
pacbio/fake_reads_005.fasta.gz
```

## (3) Map the reads.

For each reads file, we need to run 'bam_depth.sh'. The first argument is the
assembly file. The second argument is the number of the read file (this is the
line number of the reads file in input.fofn). We need to run bam_depth.sh five
times. For real data we might submit these as five separate jobs on a compute
cluster, but these files are small enough that all five should finish in less
than a minute.

```
cd genomic_data
bam_depth.sh assembly/fake_genome.fasta.gz 1
bam_depth.sh assembly/fake_genome.fasta.gz 2
bam_depth.sh assembly/fake_genome.fasta.gz 3
bam_depth.sh assembly/fake_genome.fasta.gz 4
bam_depth.sh assembly/fake_genome.fasta.gz 5
```

We now have a new subtree 'halfdeep', and the mapped read depth files
fake_reads_*.depth.dat.gz. Each run of bam_depth.sh will create a read depth
file for the corresponding reads file. Note that the alignment files themselves
are not saved, only the files containing coverage depth.

```
.
└── genomic_data
   │ input.fofn
   ├── assembly
   │   ├── fake_genome.fasta.gz
   │   ├── fake_genome.idx
   ├── pacbio
   │   ├── fake_reads_001.fasta.gz
   │   ├── fake_reads_002.fasta.gz
   │   ├── fake_reads_003.fasta.gz
   │   ├── fake_reads_004.fasta.gz
   │   ├── fake_reads_005.fasta.gz
   ├── >>> halfdeep <<<
   │   ├── >>> fake_genome <<<
   │   |   ├── >>> mapped_reads <<<
   │   |   │   ├── >>> fake_reads_001.depth.dat.gz <<<
   │   |   │   ├── >>> fake_reads_002.depth.dat.gz <<<
   │   |   │   ├── >>> fake_reads_003.depth.dat.gz <<<
   │   |   │   ├── >>> fake_reads_004.depth.dat.gz <<<
   │   |   │   ├── >>> fake_reads_005.depth.dat.gz <<<
```

## (4) Combine the coverage depth files and identify half-deep intervals

```
cd genomic_data
halfdeep.sh assembly/fake_genome.fasta.gz
```

This has produced four files -- scaffold_lengths.dat, depth.dat.gz,
percentile_commands.sh, and halfdeep.dat. These will be used in the next step
to produce a plot.

```
.
└── genomic_data
   │ input.fofn
   ├── assembly
   │   ├── fake_genome.fasta.gz
   │   ├── fake_genome.idx
   ├── pacbio
   │   ├── fake_reads_001.fasta.gz
   │   ├── fake_reads_002.fasta.gz
   │   ├── fake_reads_003.fasta.gz
   │   ├── fake_reads_004.fasta.gz
   │   ├── fake_reads_005.fasta.gz
   ├── halfdeep
   │   ├── fake_genome
   │   |   ├── mapped_reads
   │   |   │   ├── fake_reads_001.depth.dat.gz
   │   |   │   ├── fake_reads_002.depth.dat.gz
   │   |   │   ├── fake_reads_003.depth.dat.gz
   │   |   │   ├── fake_reads_004.depth.dat.gz
   │   |   │   ├── fake_reads_005.depth.dat.gz
   │   │   ├── >>> scaffold_lengths.dat <<<
   │   │   ├── >>> depth.dat.gz <<<
   │   │   ├── >>> percentile_commands.sh <<<
   │   │   ├── >>> halfdeep.dat <<<
```

## (5) Plotting

Open R with the working directory at genomic_data/halfdeep/fake_genome.

Load the HalfDeep functions in R:
```
source("path_to_halfdeep/halfdeep.r")
```

Then read in the data and prepare for plotting:
```
scaffolds         = read_scaffold_lengths("scaffold_lengths.dat")
scaffoldToOffset  = linearized_scaffolds(scaffolds)
depth             = read_depth("depth.dat.gz",scaffoldToOffset)
halfDeep          = read_halfdeep("halfdeep.dat",scaffoldToOffset)
percentileToValue = read_percentiles("percentile_commands.sh")

assembly = "fake_genome"
```

To plot to the screen:
```
halfdeep_plot(scaffolds,depth,halfDeep,percentileToValue,assembly)
```

To plot to a file (currently only supports pdf):
```
halfdeep_plot(scaffolds,depth,halfDeep,percentileToValue,assembly,plotFilename="half_deep.dat.pdf")
```

This has created the file half_deep.dat.pdf.

```
.
└── genomic_data
   │ input.fofn
   ├── assembly
   │   ├── fake_genome.fasta.gz
   │   ├── fake_genome.idx
   ├── pacbio
   │   ├── fake_reads_001.fasta.gz
   │   ├── fake_reads_002.fasta.gz
   │   ├── fake_reads_003.fasta.gz
   │   ├── fake_reads_004.fasta.gz
   │   ├── fake_reads_005.fasta.gz
   ├── halfdeep
   │   ├── fake_genome
   │   |   ├── mapped_reads
   │   |   │   ├── fake_reads_001.depth.dat.gz
   │   |   │   ├── fake_reads_002.depth.dat.gz
   │   |   │   ├── fake_reads_003.depth.dat.gz
   │   |   │   ├── fake_reads_004.depth.dat.gz
   │   |   │   ├── fake_reads_005.depth.dat.gz
   │   │   ├── scaffold_lengths.dat
   │   │   ├── depth.dat.gz
   │   │   ├── percentile_commands.sh
   │   │   ├── halfdeep.dat
   │   │   ├── >>> half_deep.dat.pdf <<<
```

The plot is shown below. The horizontal axis is positions along the assembly.
The vertical axis is coverage depth. Each dot is the average depth for a 1Kbp
window (non-overlapping). Median depth is 14.6X, as indicated in the second
line of the caption and as a dashed blue line. The plot is clipped at 1.5 times
the median, 21.8X (for real data we would have a few very deeply covered
intervals where the assembly has collapsed repeats; since we aren't interested
in those we clip the plot).

The 40th and 60th percentiles of depth are 10.4X and 17.0X. These are divided
by two and reported as 5.2 and 8.5 and indicated by dashed blue lines. These
two lines -- half 40th and half 60th -- form a band. Points within that band
are evidence of a half-depth interval.

Ignoring the red area, visually we can see that most of the FAKE1 scaffold
falls in or near the half-depth band. But there is an interval in the left
half of the scaffold that appears to be full depth.

HalfDeep attempts to automatically detect half-depth intervals. It shows these
intervals in red. For this dataset it has incorrectly called all of the FAKE1
scaffold as half-deep. This detection is still a  work in progress.

<img src="https://github.com/makovalab-psu/HalfDeep/blob/master/example/genomic_data/halfdeep/fake_genome/half_deep.dat.jpg?raw=true">

