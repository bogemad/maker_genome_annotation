#!/usr/bin/env python

import sys, re

def test_gm_in_line(gene_model_list, line):
	for gm in gene_model_list:
		if gm in line:
			return False
	return True

def read_error_file(error_file):
	gene_model_list = []
	with open(error_file) as error_data:
		for line in error_data:
			sO = re.search(r'(contig_[0-9]+_[0-9]+\-[0-9]+)', line)
			if sO:
				gene_model_list.append(sO.group(1))
	return gene_model_list

def write_filtered_ts(gene_model_list, filtered_training_set, training_set):
	with open(filtered_training_set, 'w') as ts2:
		with open(training_set) as ts:
			for line in ts:
				if line.startswith("LOCUS"):
					print_switch = test_gm_in_line(gene_model_list, line)
				if print_switch == True:
					ts2.write(line)

def main():
	error_file = sys.argv[1]
	training_set = sys.argv[2]
	filtered_training_set = sys.argv[3]
	gene_model_list = read_error_file(error_file)
	print(gene_model_list)
	write_filtered_ts(gene_model_list, filtered_training_set, training_set)


main()