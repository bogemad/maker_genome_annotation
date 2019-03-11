#!/usr/bin/env python

import sys
from Bio import SeqIO

def main():
	zff = sys.argv[1]
	fasta = sys.argv[2]
	new_fasta = sys.argv[3]
	keep_contigs = [ line[1:].strip() for line in open(zff) if line.startswith('>') ]
	recs = SeqIO.parse(fasta,'fasta')
	new_recs = [ rec for rec in recs if rec.id in keep_contigs ]
	SeqIO.write(new_recs,new_fasta,'fasta')

main()