#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./bam_depth <ref>"
	echo "    Assumes we have <ref> and <input.fofn> in the same dir"
	exit -1
fi

#module load minimap2/2.11
#module load samtools

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

cpus=$SLURM_CPUS_PER_TASK

# Unless specified, use slurm array task id for input line num.
if [ -z $2 ]; then
	i=$SLURM_ARRAY_TASK_ID
else
	i=$2
fi

if [ ! -d $ref ]; then
	mkdir $ref
fi

qry=`sed -n ${i}p input.fofn`

out=`basename $qry`
out=`echo $out | sed 's/.fasta.$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
out=halfdeep/$ref/mapped_reads/$out

if [ ! -d halfdeep/$ref ]; then
	mkdir halfdeep/$ref
fi

if [ ! -d halfdeep/$ref/mapped_reads ]; then
	mkdir halfdeep/$ref/mapped_reads
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
	minimap2 -x map-pb -a -t $cpus $ref.idx $qry | samtools view -hb - > $out.bam"
	minimap2 -x map-pb -a -t $cpus $ref.idx $qry | samtools view -hb - > $out.bam
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


echo "Compute coverage depth for $out.bam"

echo "\
samtools depth -Q 1 $out.sort.bam | gzip > $out.depth.dat.gz"
samtools depth -Q 1 $out.sort.bam | gzip > $out.depth.dat.gz

rm $out.bam
