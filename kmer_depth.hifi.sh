#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./kmer_depth.hifi.sh <ref>"
	echo "    Assumes we have <ref> and input.fofn in the current directory,"
	echo "    and that <ref> has an associated file listing its contig names,"
	echo "    in order"
	exit -1
fi

# FastK and map_names are prerequisites

refFull=$1
if [ ! -e ${refFull} ]; then
	echo "reference ${refFull} does not exist. Exit."
	exit -1
fi

ref=`echo ${refFull} | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fsa_nt$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' | sed 's/.fsa_nt.gz$//g'`
refBase=`basename ${ref}`

#cpus=3  # allow user to set this?
#cpus=$SLURM_CPUS_PER_TASK

kmerSize=40     # allow user to set this?
minAbundance=3  # allow user to set this?

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

# the intended directory for the k-mers table is ${outDir}; we will be creating
# symlinks in that directory to the actual fastqs, because we don't want FastK
# to litter our actual fastq directory with all its invisible files; here, we
# create a list of all the symlinks (which may not yet exist); and the first of
# these will be used by FastK as the name of its results

……… can't do two `s here
fastqFiles=`cat input.fofn \
  | while read f ; do
      fBase=$(basename ${f})
      echo ${outDir}/${fBase}
      done`
out=`echo ${fastqFiles} | awk '{ print $1 }'`
out=`echo ${out} | sed 's/.fastq$//g' | sed 's/.fastq.gz$//g'`

# make sure the final product has not already been computed

finalProduct=halfdeep/${refBase}/kmer_profile.depth.dat.gz

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
	cat input.fofn \
	  | while read f ; do
	      fBase=`basename ${f}`
	      ln -s ../../${f} ${outDir}/${fBase}
	      done

	echo "\
	FastK -v -k${kmerSize} -t${minAbundance} -P${tmpDir} ${fastqFiles}"
    time FastK -v -k${kmerSize} -t${minAbundance} -P${tmpDir} ${fastqFiles}
fi

if [ ! -e ${out}.ktab ]; then
	echo "Something went wrong with FastK table build. ${out}.ktab does not exist. Exit."
	exit -1
fi

# the intended directory for the k-mers table is also ${outDir}; we will be
# creating a symlink in that directory to the actual assembly, because we don't
# want FastK to litter our actual assembly directory with its files

# create k-mer contig profiles

if [ -e ${outDir}/${refBase}.prof ]; then
	echo "${outDir}/${refBase}.prof found. Skip FastK k-mer profile build."
else
	echo "Using FastK to build k-mer profiles"

    # create the symlink mentioned above
	refBaseWithExt=`basename ${refFull}`
	ln -s ../../${refFull} ${outDir}/${refBaseWithExt}

	echo "\
    FastK -v -k${kmerSize} -p:${out}.ktab -P${tmpDir} ${outDir}/${refBaseWithExt}"
    time FastK -v -k${kmerSize} -p:${out}.ktab -P${tmpDir} ${outDir}/${refBaseWithExt}
fi

if [ ! -e ${outDir}/${refBase}.prof ]; then
	echo "Something went wrong with FastK k-mer profile build. ${outDir}/${refBase}.prof does not exist. Exit."
	exit -1
fi

# convert k-mer contig profiles to same format as a bam depth file

cat ${ref}.names \
 | awk '{ print ++n,$1; }' \
 > ${outDir}/${refBase}.num_to_name

echo "\
Profex ${outDir}/${refBase}.prof 1-# ..."
time Profex ${outDir}/${refBase}.prof 1-# \
  | awk '/^Read/ { readNum = $2 } \
        !/^Read/ { if ((NF>0)&&(readNum != "")&&($2 > 0)) print readNum,$1,$2 }' \
  | tr -d ":" \
  | map_names --key=1 ${outDir}/${refBase}.num_to_name \
  | gzip \
  > ${finalProduct}

if [ ! -e ${finalProduct} ]; then
	echo "Something went wrong with k-mer profile conversion. ${finalProduct} does not exist. Exit."
	exit -1
fi

# $$$ there are a lot of intermediate files that should be deleted

