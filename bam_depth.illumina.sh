#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./bam_depth <ref>"
	echo "    Assumes we have <ref> and <input.illumina.fofn> in the same dir"
	exit -1
fi

#module load bwa
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

if [ ! -e $ref.fa ]; then
	echo "reference index $ref.fa does not exist. Exit."
	exit -1
fi

if [ ! -e $ref.fa.bwt ]; then
	echo "reference index $ref.fa.bwt does not exist. Exit."
	exit -1
fi

if [ ! -e input.illumina.fofn ]; then
	echo "input.illumina.fofn does not exist. Exit."
	exit -1
fi

qry=`sed -n ${i}p input.illumina.fofn`
numQry=`echo ${qry} | awk '{ print NF }'`

if [ "$numQry" -eq 1 ]; then
	out=`basename $qry`
else
	qry1=`echo ${qry} | awk '{ print $1 }'`
	out=`basename $qry1`
fi

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
	bwa mem -t $cpus $ref.fa $qry | samtools view -hb - > $out.bam"
	bwa mem -t $cpus $ref.fa $qry | samtools view -hb - > $out.bam
fi

if [ ! -e $out.bam ]; then
	echo "Something went wrong with bwa mem. Intermediate file $out.bam does not exist. Exit."
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

rm $out.bam
rm $out.sort.bam
rm $out.sort.bam.bai
