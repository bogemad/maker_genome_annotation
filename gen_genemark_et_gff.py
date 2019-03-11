#!/usr/bin/env python

import sys

zff = open(sys.argv[1])
out_gff = open(sys.argv[2], 'w')

for line in zff:
	if line.startswith('>'):
		contig = line.strip()[1:]
		continue
	data = line.strip().split('\t')
	if int(data[1]) < int(data[2]):
		strand = '+'
	elif int(data[1]) > int(data[2]):
		strand = '-'
	else:
		print("Strand error: start = {}, end = {}".format(data[1],data[2]))
		sys.exit(1)
	if data[0] == 'Einit':
		if strand == '+':
			intron_start = int(data[2]) + 1
		elif strand == '-':
			intron_end = int(data[2]) - 1
	elif data[0] == 'Exon':
		if strand == '+':
			intron_end = int(data[1]) - 1
		elif strand == '-':
			intron_start = int(data[1]) + 1
		out_gff.write("{0}\tPASA\tintron\t{1}\t{2}\t10\t{3}\t.\t{4}\n".format(contig,intron_start,intron_end, strand, data[3]))
		if strand == '+':
			intron_start = int(data[2]) + 1
		elif strand == '-':
			intron_end = int(data[2]) - 1
	elif data[0] == 'Eterm':
		if strand == '+':
			intron_end = int(data[1]) - 1
		elif strand == '-':
			intron_start = int(data[1]) + 1
		out_gff.write("{0}\tPASA\tintron\t{1}\t{2}\t10\t{3}\t.\t{4}\n".format(contig,intron_start,intron_end, strand, data[3]))
	elif data[0] == 'Esngl':
		continue
