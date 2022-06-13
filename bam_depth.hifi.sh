#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./bam_depth.hifi <ref> [<number>]"
	echo "    Assumes we have <ref> and input.fofn in the current dir"
	echo "    If <number> is not given, SLURM_ARRAY_TASK_ID is used"
	echo "    <number> is the line number in input.fofn of the file to process"
	exit -1
fi

#module load minimap2/2.11
#module load samtools

ref=$1
if [ ! -e $ref ]; then
	echo "reference $ref does not exist. Exit."
	exit -1
fi

ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fsa_nt$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' | sed 's/.fsa_nt.gz$//g'`
refbase=`basename $ref`

cpus=3  # fix this!
#cpus=$SLURM_CPUS_PER_TASK

# Unless specified, use slurm array task id for input line num.
if [ -z $2 ]; then
	i=$SLURM_ARRAY_TASK_ID
else
	i=$2
fi

if [ ! -e $ref.idx ]; then
	echo "reference index $ref.idx does not exist. Exit."
	exit -1
fi

if [ ! -e input.fofn ]; then
	echo "input.fofn does not exist. Exit."
	exit -1
fi

qry=`sed -n ${i}p input.fofn`

out=`basename $qry`
out=`echo $out | sed 's/.fasta.$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' | sed 's/.fastq.gz$//g'`
out=halfdeep/$refbase/mapped_reads/$out

if [ ! -d halfdeep ]; then
	mkdir halfdeep
fi

if [ ! -d halfdeep/$refbase ]; then
	mkdir halfdeep/$refbase
fi

if [ ! -d halfdeep/$refbase/mapped_reads ]; then
	mkdir halfdeep/$refbase/mapped_reads
fi


if [ -e $out.depth.dat.gz ]; then
	echo "$out.depth.dat.gz found. Exit."
	exit 0
fi

if [ -e $out.bam ]; then
	echo "$out.bam found. Skip alignment."
else
	echo "Start aligning $qry to $ref.idx"

	echo "\
	minimap2 -x map-hifi -a -t $cpus --split-prefix temp $ref.idx $qry | samtools view -hb - > $out.bam"
	minimap2 -x map-hifi -a -t $cpus --split-prefix temp $ref.idx $qry | samtools view -hb - > $out.bam
fi

if [ ! -e $out.bam ]; then
	echo "Something went wrong with minimap2. Intermediate file $out.bam does not exist. Exit."
	exit -1
fi

if [ -e $out.sort.bam ]; then
	echo "$out.sort.bam found. Skip sort."
else
	echo "Sort $out.bam"

	echo "\
	samtools sort -T $out.tmp -O bam -o $out.sort.bam $out.bam"
	samtools sort -T $out.tmp -O bam -o $out.sort.bam $out.bam
	samtools index $out.sort.bam
fi

if [ ! -e $out.sort.bam ]; then
	echo "Something went wrong with samtools sort. Intermediate file $out.sort.bam does not exist. Exit."
	exit -1
fi

echo "Compute coverage depth for $out.bam"

echo "\
samtools depth -Q 1 $out.sort.bam | gzip > $out.depth.dat.gz"
samtools depth -Q 1 $out.sort.bam | gzip > $out.depth.dat.gz

if [ ! -e $out.depth.dat.gz ]; then
	echo "Something went wrong with samtools depth. $out.depth.dat.gz does not exist. Exit."
	exit -1
fi

#rm $out.bam
#rm $out.sort.bam
#rm $out.sort.bam.bai
