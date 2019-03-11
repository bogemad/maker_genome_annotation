#!/usr/bin/env python

import sys

file = open(sys.argv[1])
outfile = open(sys.argv[2],'w')
lst_asmbl = ""
asmbl_list = []
first_line = True

def output_fixed_zff(asmbl_list, outfile):
	for (i, (start, end, asmbl)) in enumerate(asmbl_list):
		if i == 0:
			outfile.write("Einit\t{0}\t{1}\t{2}\n".format(start, end, asmbl))
		elif i == (len(asmbl_list)-1):
			outfile.write("Eterm\t{0}\t{1}\t{2}\n".format(start, end, asmbl))
		else:
			outfile.write("Exon\t{0}\t{1}\t{2}\n".format(start, end, asmbl))


for line in file:
	if line.startswith(">"):
		if first_line == False:
			output_fixed_zff(asmbl_list, outfile)
			asmbl_list = []
		else:
			first_line = False
		outfile.write(line)
		lst_asmbl = ""
		continue
	data = line.strip().split('\t')
	if data[0] == "Esngl":
		outfile.write(line)
		lst_asmbl = ""
		continue
	start = data[1]
	end = data[2]
	asmbl = data[3]
	if asmbl == lst_asmbl or lst_asmbl == "":
		asmbl_list.append((start, end, asmbl))
	else:
		output_fixed_zff(asmbl_list, outfile)
		asmbl_list = [(start, end, asmbl)]
	lst_asmbl = asmbl
output_fixed_zff(asmbl_list, outfile)

file.close()
outfile.close()