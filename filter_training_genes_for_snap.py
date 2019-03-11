#!/usr/bin/env python

import sys, re
from collections import defaultdict

def tryint(s):
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

def test_gm_in_line(error_list, line):
	for gm in error_list:
		if gm in line:
			return False
	return True

def read_error_file(error_file):
	error_list = []
	with open(error_file) as error_data:
		for line in error_data:
			sO = re.search(r'(contig_[0-9]+_[0-9]+\-[0-9]+)', line)
			if sO:
				error_list.append(sO.group(1))
	return error_list


def extract_group(field):
	x = field.split(';')
	y = x[0].split('|')
	prefix = y[0].split('=')[1]
	pre_split = prefix.split('.')
	prefix = pre_split[len(pre_split)-1]
	suffix = 'g.' + y[1].split('.')[1]
	group = prefix + '|' + suffix
	if not group.startswith('asmbl'):
		print("Group name error: {}".format(group))
		sys.exit(1)
	return group, x[0][3:]


def modify_label(field):
	if field == "gene":
		return "ignore"
	if field == "mRNA":
		return "ignore"
	if field == "five_prime_utr":
		return "ignore"
	if field == "exon":
		return "weird_Exon"
	if field == "CDS":
		return "Exon"
	if field == "three_prime_utr":
		return "ignore"


def read_gff(training_set):
	gff_d = defaultdict(lambda : defaultdict(list))
	exon_num = 0
	with open(training_set) as gff_handle:
		for line in gff_handle:
			if line.strip() == '':
				exon_num = 1
				continue
			data = line.strip().split('\t')
			contig = data[0]
			group, id = extract_group(data[8])
			label = modify_label(data[2])
			strand = data[6]
			if strand == '-':
				start = int(data[4])
				end = int(data[3])
			else:
				start = int(data[3])
				end = int(data[4])
			if label == "Exon":
				gff_d[contig][group].append((id,label,start,end,strand,exon_num))
				exon_num += 1
	return gff_d


def write_filtered_ts(error_list, filtered_training_set, training_set):
	with open(filtered_training_set, 'w') as ts2:
		with open(training_set) as ts:
			for line in ts:
				if line.startswith("LOCUS"):
					print_switch = test_gm_in_line(error_list, line)
				if print_switch == True:
					ts2.write(line)


def filter_gff_dict(gff_d, complete_list):
	new_gff_d = defaultdict(lambda : defaultdict(list))
	for contig in list(gff_d):
		for group in gff_d[contig]:
			if group in complete_list:
				new_gff_d[contig][group] = gff_d[contig][group]
	return new_gff_d


def write_zff(filtered_training_set,gff_d,complete_list):
	with open(filtered_training_set, 'w') as ts2:
		iter_list = list(gff_d)
		sort_nicely(iter_list)
		for contig in iter_list:
			ts2.write(">{}\n".format(contig))
			for group in list(gff_d[contig]):
				for item in gff_d[contig][group]:
					id,label,start,end,strand,exon_num = item
					num_exons = sum([1 for x in gff_d[contig][group] if x[1] == "Exon"])
					print("{}\t{}\t{}".format(group,exon_num,num_exons))
					if num_exons == 1:
						if label == "Exon":
							label = "Esngl"
					if num_exons > 1:
						if exon_num == 1:
							label = "Einit"
						if exon_num == num_exons:
							label = "Eterm"
					ts2.write("{0}\t{1}\t{2}\t{3}\n".format(label, start, end, group))

def main():
	error_file = sys.argv[1] #train.err from augustus run
	training_set = sys.argv[2] #{database_name}.assemblies.fasta.transdecoder.genome.gff3 from augustus run
	filtered_training_set = sys.argv[3] #filtered output file
	complete_genes = sys.argv[4] #pasa.complete.lst from augustus run
	print("Reading error file...")
	error_list = read_error_file(error_file)
	print("Reading gff file...")
	gff_d = read_gff(training_set)
	print("Reading complete genes file...")
	complete_list = [line.strip().replace('|m.','|g.')[:-1] for line in open(complete_genes)]
	# print(complete_list)
	# print([ x for sub in list(gff_d) for x in list(gff_d[sub])])
	print("Filtering data...")
	gff_d = filter_gff_dict(gff_d, complete_list)
	print("Writing zff...")
	write_zff(filtered_training_set,gff_d,complete_list)


main()