#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./halfdeep.kmers <ref> <numSamples>"
	echo "    Assumes we have <ref>.lengths and input.fofn in the current"
	echo "    directory"
	exit -1
fi

echo WORK IN PROGRESS
……… yank that

shortWindowSize=1K           # allow user to set this?
longWindowSize=100K          # allow user to set this?
detectionThreshold=0.10      # allow user to set this?
gapFill=500K                 # allow user to set this?

refIn=$1
numSamples=$2


ref=`echo ${refIn} | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fsa_nt$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' | sed 's/.fsa_nt.gz$//g'`
refBase=`basename ${ref}`

if [ ! -d halfdeep ]; then
	echo "    halfdeep dir not found, has kmer_depth.hifi not been run?"
	exit -1
fi

if [ ! -d halfdeep/${refBase} ]; then
	echo "    halfdeep/${refBase} dir not found, has kmer_depth.hifi not been run?"
	exit -1
fi

if [ ! -e halfdeep/${refBase}/kmer_profile.depth.dat.gz ]; then
	echo "    halfdeep/${refBase}/kmer_profile.depth.dat.gz not found, has kmer_depth.hifi not been run?"
	exit -1
fi


if [ -e halfdeep/${refBase}/scaffold_lengths.dat ]; then
	echo "halfdeep/${refBase}/scaffold_lengths.dat found. Skip scaffold lengths step."
else
	echo "Collecting scaffold lengths from ${refIn}"
	echo "\
    scaffold_lengths.py ${refIn} > halfdeep/${refBase}/scaffold_lengths.dat"
    scaffold_lengths.py ${refIn} > halfdeep/${refBase}/scaffold_lengths.dat
fi

if [ ! -e halfdeep/${refBase}/scaffold_lengths.dat ]; then
	echo "Something went wrong with scaffold_lengths. halfdeep/${refBase}/scaffold_lengths.dat does not exist. Exit."
	exit -1
fi

……… this needs to loop
………

    sampleNum=0
    time while [ ${sampleNum} -lt ${numSamples} ]; do
        sampleNum=$((sampleNum+1))
        subsampleFofn=input.${sampleNum}_of_${numSamples}.fofn
………
        outName=`head -n 1 ${subsampleFofn}`
        outName=`echo ${outName} | sed 's/.fastq$//g' | sed 's/.fastq.gz$//g'`
………
        done


if [ -e halfdeep/${refBase}/depth.dat.gz ]; then
	echo "halfdeep/${refBase}/depth.dat.gz found. Skip depth-windowing step."
else
	echo "Reduce the reads-file depths to ${shortWindowSize} windows ...

	# nota bene: the awk command converts samtools' scaffold,position,depth to
	#            scaffold,start,end,depth that genodsp wants
	echo "\
    gzip -dc halfdeep/${refBase}/kmer_profile.depth.dat.gz | ... > halfdeep/${refBase}/depth.dat.gz"
    time gzip -dc halfdeep/${refBase}/kmer_profile.depth.dat.gz \
      | awk '/^#/ { print $0; }
            !/^#/ { print $1,$2,$2,$3 }' \
	  | genodsp --origin=one --uncovered:hide --precision=3 \
	      --report=comments --progress=input:10M \
		  --chromosomes=halfdeep/${refBase}/scaffold_lengths.dat \
		  = sum --window=$shortWindowSize --denom=actual \
	  | gzip \
	  > halfdeep/${refBase}/depth.dat.gz
fi

if [ ! -e halfdeep/${refBase}/depth.dat.gz ]; then
	echo "Something went wrong with genodsp. halfdeep/${refBase}/depth.dat.gz does not exist. Exit."
	exit -1
fi


echo "Computing percentiles of depth distribution"

rm -f halfdeep/${refBase}/percentile_commands.sh
time gzip -dc halfdeep/${refBase}/depth.dat.gz \
  | genodsp --origin=one --uncovered:hide --precision=3 --nooutput \
	  --chromosomes=halfdeep/${refBase}/scaffold_lengths.dat \
	  = percentile 40..60by10 --min=1/inf --report:bash \
  > halfdeep/${refBase}/temp.percentile_commands

if [ ! -e halfdeep/${refBase}/temp.percentile_commands ]; then
	echo "Something went wrong with genodsp. halfdeep/${refBase}/temp.percentile_commands does not exist. Exit."
	exit -1
fi

cat halfdeep/${refBase}/temp.percentile_commands \
  | grep "# bash command" \
  | tr "=" " " \
  | awk '{ print "export",$1"="$2 }' \
  > halfdeep/${refBase}/percentile_commands.sh

cat halfdeep/${refBase}/temp.percentile_commands \
  | grep "# bash command" \
  | tr "=" " " \
  | awk '{ print "export","half"$1"="($2/2) }' \
  | sed "s/halfpercentile/halfPercentile/" \
  >> halfdeep/${refBase}/percentile_commands.sh

if [ ! -e halfdeep/${refBase}/percentile_commands.sh ]; then
	echo "Something went wrong. halfdeep/${refBase}/percentile_commands.sh does not exist. Exit."
	exit -1
fi

rm halfdeep/${refBase}/temp.percentile_commands

echo "Identifying half-deep intervals"

source halfdeep/${refBase}/percentile_commands.sh
rm -f halfdeep/${refBase}/halfdeep.dat
time gzip -dc halfdeep/${refBase}/depth.dat.gz \
  | genodsp --origin=one --uncovered:hide --nooutputvalue \
	  --chromosomes=halfdeep/${refBase}/scaffold_lengths.dat \
	  = erase --keep:inside --min=${halfPercentile40} --max=${halfPercentile60} \
	  = binarize \
	  = dilate --right=$shortWindowSize-1 \
	  = sum --window=$longWindowSize --denom=actual \
	  = erase --max=$detectionThreshold \
	  = binarize \
	  = dilate --right=$longWindowSize-1 \
	  = close $gapFill \
  > halfdeep/${refBase}/halfdeep.dat

if [ ! -e halfdeep/${refBase}/halfdeep.dat ]; then
	echo "Something went wrong with genodsp. halfdeep/${refBase}/halfdeep.dat does not exist. Exit."
	exit -1
fi
