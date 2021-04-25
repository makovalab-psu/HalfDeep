#!/usr/bin/env python
"""
Read a fasta file and report the names and lengths of all the sequences.
"""

from sys   import argv,stdin,stdout,stderr,exit
from gzip import open as gzip_open

def main():

	assert (len(argv) in [1,2])

	(inFilename,inFile) = (None,stdin)
	if (len(argv) == 2):
		inFilename = argv[1]
		if (inFilename.endswith(".gz")) or (inFilename.endswith(".gzip")):
			inFile = gzip_open(inFilename,"rt")
		else:
			inFile = file(inFilename,"rt")

	for (name,sequence) in fasta_sequences(inFile):
		print "%s\t%d" % (name,len(sequence))

	if (inFilename != None): inFile.close()


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = sequence_name(line)
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


# sequence_name--
#	Extract the sequence name from a fasta header.

def sequence_name(s):
	s = s[1:].strip()
	if (s == ""): return ""
	else:         return s.split()[0]


if __name__ == "__main__": main()

