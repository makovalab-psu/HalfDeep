#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./kmer_depth.hifi.sh <ref> <sampleNum>/<numSamples>"
	echo "    Assumes we have <ref> and input.fofn in the current directory,"
	echo "    and that <ref> has an associated file listing its contig names,"
	echo "    in order"
	echo "    input.fofn will be partitioned into <numSamples> different"
	echo "    subsets, and processing will be performed on the <sampleNum>th"
	echo "    subset"
	exit -1
fi

# FastK, subsample, and map_names are prerequisites



kmerSize=40     # allow user to set this?
minAbundance=3  # allow user to set this?
cpus=5          # allow user to set this?  OR cpus=$SLURM_CPUS_PER_TASK


refFull=$1
subsampleSpec=$2

if [ ! -e ${refFull} ]; then
	echo "reference ${refFull} does not exist. Exit."
	exit -1
fi

ref=`echo ${refFull} | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
refBase=`basename ${ref}`
refExt=`echo ${refFull} | sed 's/^.*\.fasta$/.fasta/g' | sed 's/^.*\.fa$/.fa/g' | sed 's/^.*\.fasta.gz$/.fasta.gz/g' | sed 's/^.*\.fa.gz$/.fa.gz/g'`

if [ ! -e input.fofn ]; then
	echo "input.fofn does not exist. Exit."
	exit -1
fi

outDir=halfdeep/fastk
tmpDir=halfdeep/fastk/temp

if [ ! -d halfdeep ]; then
	mkdir halfdeep
fi

if [ ! -d halfdeep/fastk ]; then
	mkdir halfdeep/fastk
fi

if [ ! -d halfdeep/fastk/temp ]; then
	mkdir halfdeep/fastk/temp
fi

if [ ! -d halfdeep/${refBase} ]; then
	mkdir halfdeep/${refBase}
fi

if [ ! -d halfdeep/${refBase}/kmer_profiles ]; then
	mkdir halfdeep/${refBase}/kmer_profiles
fi


# extract the <sampleNum>th subset of the fastq files

sampleNum=` echo ${subsampleSpec} | tr "/" " " | awk '{ print $1 }'`
numSamples=`echo ${subsampleSpec} | tr "/" " " | awk '{ print $2 }'`
subsampleName=${sampleNum}_of_${numSamples}
subsampleFofn=input.${subsampleName}.fofn

echo "\
cat input.fofn | subsample ${subsampleSpec} > ${subsampleFofn}"
cat input.fofn | subsample ${subsampleSpec} > ${subsampleFofn}

# the intended directory for the k-mers table is ${outDir}; we will be creating
# symlinks in that directory to the actual fastqs, because we don't want FastK
# to litter our actual fastq directory with all its invisible files; here, we
# create a list of all the symlinks (which may not yet exist); and the first of
# these will be used by FastK as the name of its results

fastqFiles=`cat ${subsampleFofn} \
  | while read f ; do
      fBase=$(basename ${f})
      echo ${outDir}/${fBase}
      done`
out=`echo ${fastqFiles} | awk '{ print $1 }'`
out=`echo ${out} | sed 's/.fastq$//g' | sed 's/.fastq.gz$//g'`
outBase=`basename ${out}`

# make sure the final product has not already been computed

finalProduct=halfdeep/${refBase}/kmer_profiles/${subsampleName}.depth.dat.gz

if [ -e ${finalProduct} ]; then
	echo "${finalProduct} found. Exit."
	exit 0
fi

# build the k-mers table

if [ -e ${out}.ktab ]; then
	echo "${out}.ktab found. Skip FastK table build."
else
	echo "Using FastK to build k-mers table"

    # create the symlinks mentioned above
	cat ${subsampleFofn} \
	  | while read f ; do
	      fBase=`basename ${f}`
	      ln -s ../../${f} ${outDir}/${fBase}
	      done

	echo "\
	FastK -v -k${kmerSize} -t${minAbundance} -T${cpus} -P${tmpDir} ${fastqFiles}"
    time FastK -v -k${kmerSize} -t${minAbundance} -T${cpus} -P${tmpDir} ${fastqFiles}
fi

if [ ! -e ${out}.ktab ]; then
	echo "Something went wrong with FastK table build. ${out}.ktab does not exist. Exit."
	exit -1
fi

# the intended directory for the k-mers table is also ${outDir}; we will be
# creating a symlink in that directory to the actual assembly, because we don't
# want FastK to litter our actual assembly directory with its files

# create k-mer contig profiles

if [ -e ${outDir}/${refBase}.${subsampleName}.prof ]; then
	echo "${outDir}/${refBase}.${subsampleName}.prof found. Skip FastK k-mer profile build."
else
	echo "Using FastK to build k-mer profiles"

    # create the symlink mentioned above
	ln -s ../../${refFull} ${outDir}/${refBase}.${subsampleName}${refExt}

	echo "\
    FastK -v -k${kmerSize} -p:${out}.ktab -P${tmpDir} ${outDir}/${refBaseWithExt}"
    time FastK -v -k${kmerSize} -p:${out}.ktab -P${tmpDir} ${outDir}/${refBase}.${subsampleName}${refExt}
fi

if [ ! -e ${outDir}/${refBase}.${subsampleName}.prof ]; then
	echo "Something went wrong with FastK k-mer profile build. ${outDir}/${refBase}.${subsampleName}.prof does not exist. Exit."
	exit -1
fi

# convert k-mer contig profiles to the same format as a bam depth file

cat ${ref}.names \
 | awk '{ print ++n,$1; }' \
 > ${outDir}/${refBase}.${subsampleName}.num_to_name

echo "\
Profex ${outDir}/${refBase}.${subsampleName}.prof 1-# ..."
time Profex ${outDir}/${refBase}.${subsampleName}.prof 1-# \
  | awk '/^Read/ { readNum = $2 } \
        !/^Read/ { if ((NF>0)&&(readNum != "")&&($2 > 0)) print readNum,$1,$2 }' \
  | tr -d ":" \
  | map_names --key=1 ${outDir}/${refBase}.${subsampleName}.num_to_name \
  | gzip \
  > ${finalProduct}

if [ ! -e ${finalProduct} ]; then
	echo "Something went wrong with k-mer profile conversion. ${finalProduct} does not exist. Exit."
	exit -1
fi

# $$$ there are a lot of intermediate files that should be deleted

