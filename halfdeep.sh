#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./halfdeep <ref>"
	echo "    Assumes we have <ref>.lengths and <input.fofn> in the same dir"
	exit -1
fi


shortWindowSize=1K
longWindowSize=100K
detectionThreshold=0.10
gapFill=500K

refin=$1
ref=$refin
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
refbase=`basename $ref`

if [ ! -d halfdeep ]; then
	echo "    halfdeep dir not found, has bam_depth not been run?"
	exit -1
fi

if [ ! -d halfdeep/$refbase ]; then
	echo "    halfdeep/$refbase dir not found, has bam_depth not been run?"
	exit -1
fi

if [ ! -d halfdeep/$refbase/mapped_reads ]; then
	echo "    halfdeep/$refbase/mapped_reads dir not found, has bam_depth not been run?"
	exit -1
fi


if [ -e halfdeep/$refbase/scaffold_lengths.dat ]; then
	echo "halfdeep/$refbase/scaffold_lengths.dat found. Skip scaffold lengths step."
else
	echo "Collecting scaffold lengths from $refin"
	echo "\
    scaffold_lengths $refin > halfdeep/$refbase/scaffold_lengths.dat"
    scaffold_lengths $refin > halfdeep/$refbase/scaffold_lengths.dat
fi


echo "Combining individual reads-file depths to single track, in ${shortWindowSize} windows"

cat input.fofn \
  | while read readsname ; do
      readsname=`basename $readsname`
      readsname=`echo $out | sed 's/.fasta.$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' | sed 's/.fastq.gz$//g'`
      gzip -dc $readsname.depth.dat.gz
	  done \
  | awk '{ print $1,$2,$2,$3 }' \
  | genodsp --origin=one --uncovered:hide --precision=3 \
	  --chromosomes=halfdeep/$refbase/scaffold_lengths.dat \
	  = sum --window=$shortWindowSize --denom=actual \
  | gzip \
  > halfdeep/$refbase/depth.dat.gz

echo "Computing percentiles of depth distribution"

gzip -dc halfdeep/$refbase/depth.dat.gz \
  | genodsp --origin=one --uncovered:hide --precision=3 --nooutput \
	  --chromosomes=halfdeep/$refbase/scaffold_lengths.dat \
	  = percentile 40..60by10 --min=1/inf --report:bash \
  > halfdeep/$refbase/temp.percentile_commands

cat halfdeep/$refbase/temp.percentile_commands \
  | grep "# bash command" \
  | tr "=" " " \
  | awk '{ print "export",$1"="$2 }' \
  > halfdeep/$refbase/percentile_commands.sh

cat halfdeep/$refbase/temp.percentile_commands \
  | grep "# bash command" \
  | tr "=" " " \
  | awk '{ print "export","half"$1"="($2/2) }' \
  | sed "s/halfpercentile/halfPercentile/" \
  >> halfdeep/$refbase/percentile_commands.sh

rm halfdeep/$refbase/temp.percentile_commands

echo "Identifying half-deep intervals"

source halfdeep/$refbase/percentile_commands.sh
gzip -dc halfdeep/$refbase/depth.dat.gz \
  | genodsp --origin=one --uncovered:hide --nooutputvalue \
	  --chromosomes=halfdeep/$refbase/scaffold_lengths.dat \
	  = erase --keep:inside --min=${halfPercentile40} --max=${halfPercentile60} \
	  = binarize \
	  = dilate --right=$shortWindowSize-1 \
	  = sum --window=$longWindowSize --denom=actual \
	  = erase --max=$detectionThreshold \
	  = binarize \
	  = dilate --right=$longWindowSize-1 \
	  = close $gapFill \
  > halfdeep/$refbase/halfdeep.dat

