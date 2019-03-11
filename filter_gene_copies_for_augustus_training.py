#!/usr/bin/env python

import sys, os
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def build_blacklist(file, blacklist):
	file_handle = open(file)
	file = os.path.basename(file)
	blast_results = NCBIXML.read(file_handle)
	if len(blast_results.alignments) > 1: #check if multiple alignments
		for alignment in blast_results.alignments:
			if file[:-4] in alignment.title:
				continue #ignore if blast subject is same as query
			for hsp in alignment.hsps:
				if hsp.identities / float(hsp.align_length) >= 0.7:
					print("HSP identity = {1:.3f} to {2}: Adding {0} to blacklist...".format(file[:-4],(hsp.identities / float(hsp.align_length)),alignment.title))
					blacklist.append(file[:-4])
					return blacklist
	return blacklist


def gen_locus_blacklist(blacklist, genome):
	locus_blacklist = []
	locus_list = []
	for chrom in genome:
		for feat in chrom.features:
			if feat.type == 'CDS':
				if feat.qualifiers['protein_id'][0] in blacklist:
					locus_blacklist.append(feat.qualifiers['locus_tag'][0])
					print("Protein {0} associated with locus tag {1}. Adding {1} to locus_blacklist.".format(feat.qualifiers['protein_id'][0],feat.qualifiers['locus_tag'][0]))
				# if len(feat) < 1200:
					# locus_blacklist.append(feat.qualifiers['locus_tag'][0])
					# print("Protein {0} shorter than 400 amino acids. Adding {1} to locus_blacklist.".format(feat.qualifiers['protein_id'][0],feat.qualifiers['locus_tag'][0]))
	return locus_blacklist


def filter_genome(locus_blacklist, genome):
	log_dict = {}
	del_dict = {}
	for chrom in genome:
		new_feature_list = []
		for index, feat in enumerate(chrom.features):
			if feat.type in list(log_dict):
				log_dict[feat.type] += 1
			else:
				log_dict[feat.type] = 1
				del_dict[feat.type] = 0
			if "locus_tag" in feat.qualifiers:
				if feat.qualifiers['locus_tag'][0] in locus_blacklist:
					# chrom.features.pop(index)
					print("Filtered {0} associated with locus_tag {1}.".format(feat.type,feat.qualifiers['locus_tag'][0]))
					del_dict[feat.type] += 1
				else:
					new_feature_list.append(feat)
			else:
				del_dict[feat.type] += 1
				print("Filtered {0}. Not associated with recognised genes.".format(feat.type))
		chrom.features = new_feature_list
	print("Filtering completed.")
	for feat_type in list(log_dict):
		print("{0} {1} found, {2} {1} deleted.".format(log_dict[feat_type], feat_type, del_dict[feat_type]))
	return genome


def divide_genome(genome):
	output_genes = []
	new_record = ""
	for chrom in genome:
		for feat in chrom.features:
			if feat.type == 'gene':
				if new_record:
					output_genes.append(new_record)
				record_sequence = feat.extract(chrom.seq)
				new_record = SeqRecord(record_sequence)
				new_record.id = feat.qualifiers['locus_tag'][0]
				new_record.description = chrom.description
				new_record.name = feat.qualifiers['locus_tag'][0]
				new_location = FeatureLocation(feat.location.start - feat.location.start, feat.location.end - feat.location.start )
				new_feat = SeqFeature(location = new_location, type = feat.type, strand = feat.strand, ref = feat.ref, ref_db = feat.ref_db)
				new_feat.qualifiers = feat.qualifiers
				new_record.features = [new_feat]
				new_record.annotations['topology'] = chrom.annotations['topology']
				new_record.annotations['date'] = chrom.annotations['date']
				new_record.annotations['taxonomy'] = chrom.annotations['taxonomy']
				new_record.annotations['source'] = chrom.annotations['source']
				new_record.annotations['organism'] = chrom.annotations['organism']
				new_record.annotations['sequence_version'] = chrom.annotations['sequence_version']
				new_record.annotations['data_file_division'] = chrom.annotations['data_file_division']
				new_record.annotations['references'] = chrom.annotations['references']
			else:
				new_location = FeatureLocation(feat.location.start - feat.location.start, feat.location.end - feat.location.start )
				new_feat = SeqFeature(location = new_location, type = feat.type, strand = feat.strand, ref = feat.ref, ref_db = feat.ref_db)
				new_feat.qualifiers = feat.qualifiers
				new_record.features.append(new_feat)
	return output_genes


def main():
	xml_dir = sys.argv[1]
	gb_genome = sys.argv[2]
	outfile = sys.argv[3]
	blacklist = []
	for file in [ file for file in os.listdir(xml_dir) if file.endswith('.xml') ]:
		file = os.path.join(xml_dir,file)
		blacklist = build_blacklist(file, blacklist)
	genome = [ chrom for chrom in SeqIO.parse(gb_genome,"genbank") ]
	print("Number of features in genome = {}".format(sum(len(chrom.features) for chrom in genome)))
	locus_blacklist = gen_locus_blacklist(blacklist, genome)
	genome = filter_genome(locus_blacklist, genome)
	# output_genes = divide_genome(genome)
	print("Number of features in genome after filtering = {}".format(sum(len(chrom.features) for chrom in genome)))
	SeqIO.write(genome, outfile, "genbank")


main()