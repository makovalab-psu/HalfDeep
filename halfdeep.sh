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

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

if [ ! -d $ref ]; then
	mkdir $ref
fi


if [ -e $ref/scaffold_lengths.dat ]; then
	echo "$ref/scaffold_lengths.dat found. Skip scaffold lengths step."
else
	echo "Collecting scaffold lengths"
	echo "\
    fasta_lengths $1 > $ref/scaffold_lengths.dat"
    fasta_lengths $1 > $ref/scaffold_lengths.dat
fi


echo "Combining individual bam depths to single track, in short windows"

cat input.fofn \
  | while read bamName ; do
	  gzip -dc $bamName.depth.dat.gz
	  done \
  | awk '{ print $1,$2,$2,$3 }' \
  | genodsp --origin=one --uncovered:hide --precision=3 \
	  --chromosomes=$ref/scaffold_lengths.dat \
	  = sum --window=$shortWindowSize --denom=actual \
  | gzip \
  > $ref/depth.dat.gz

echo "Computing percentiles of depth distribution"

gzip -dc $ref/depth.dat.gz \
  | genodsp --origin=one --uncovered:hide --precision=3 --nooutput \
	  --chromosomes=$ref/scaffold_lengths.dat \
	  = percentile 40..60by10 --min=1/inf --report:bash \
  > $ref/temp.percentile_commands

cat $ref/temp.percentile_commands \
  | grep "# bash command" \
  | tr "=" " " \
  | awk '{ print "export",$1"="$2 }' \
  > $ref/percentile_commands.sh

cat $ref/temp.percentile_commands \
  | grep "# bash command" \
  | tr "=" " " \
  | awk '{ print "export","half"$1"="($2/2) }' \
  | sed "s/halfpercentile/halfPercentile/" \
  >> $ref/percentile_commands.sh

rm $ref/temp.percentile_commands

echo "Identifying half-deep intervals"

source $ref/percentile_commands.sh
gzip -dc $ref/depth.dat.gz \
  | genodsp --origin=one --uncovered:hide --nooutputvalue \
	  --chromosomes=$ref/scaffold_lengths.dat \
	  = erase --keep:inside --min=${halfPercentile40} --max=${halfPercentile60} \
	  = binarize \
	  = dilate --right=$shortWindowSize-1 \
	  = sum --window=$longWindowSize --denom=actual \
	  = erase --max=$detectionThreshold \
	  = binarize \
	  = dilate --right=$longWindowSize-1 \
	  = close $gapFill \
  > $ref/halfdeep.dat

