#!/usr/bin/python3

import os
import sys
import re
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import gzip


def get_parser():
	"""
	This function parses the arguments provided by the user
	:return: a dictionary having a key for each arguments
	:rtype: dict
	"""
	parser = argparse.ArgumentParser(
	    description = 'This script is used to calculate the scaffold length, N50, GC% and N% of the genome.\n',
	    usage = 'python3 Calculate_N50_gc_N_of_genome.py -i <genome.fa/genome.fa.gz> -o <output.txt> -m <N50|len>',
	    formatter_class = RawTextHelpFormatter, add_help = False)

	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')

	required.add_argument(
	    '-i', '--in', dest = 'in', required = True, metavar = 'FASTA FILE', help = 'Input sequence file in FASTA format')

	required.add_argument(
	    '-o', '--out', dest = 'out', required = True, metavar = 'OUTPUT',
	    help = 'Output folders and files will be labelled with this name. WARNING: do not provide a path')

	required.add_argument(
	    '-m', '--mode', dest = 'mode', required = True, metavar = 'MODE',
	    help = 'Specify which analysis mode to run.\n'
	  	       'There are two valid modes:\n- len for caculating the length, GC%%, and N%% of each scaffold or genome\n'
	           '- N50 for calcuting the scaffold N50 of the genome')

	optional.add_argument('-h', '--help', action = "help", help = "Show this help message and exit")

	return parser.parse_args()


def Read_fasta(filename):
	if filename.endswith('gz'):
		Fasta = gzip.open(filename, 'rt')
	else:
		Fasta = open(filename, 'r')

	seq = {}

	for line in Fasta:
		if line.startswith('>'):
			seq_id = line.replace('>', '').split()[0]
			seq[seq_id] = [0, 0, 0]
		else:
			seq[seq_id][0] += len(line.replace('\n', '').strip())
			seq[seq_id][1] += len(re.findall('[GCgc]', line))
			seq[seq_id][2] += len(re.findall('[Nn]', line))

	Fasta.close()

	return seq


def Count_len_gc_n(in_dict, outfile):
	o = open(outfile, "w")

	df = pd.DataFrame(data = in_dict).T
	df = df.sort_values(by = 0, ascending = False)
	df.columns = ['length', 'GC_length', 'N_length']

	tot_len = df['length'].sum()
	tot_gc = df['GC_length'].sum()
	tot_N_len = df['N_length'].sum()

	df.loc['Sum'] = df.apply(lambda x: x.sum())
	df['GC%'] = (df['GC_length'] / df['length']).round(2)
	df['N%'] = (df['N_length'] / df['length']).round(2)
	cols = df.columns.tolist()
	cols = ['length', 'GC%', 'N_length', 'N%']
	print(df[cols].to_string(), file = o)

	num_ge_100 = (df['length'] >= 100).sum()
	num_ge_200 = (df['length'] >= 200).sum()
	num_ge_1000 = (df['length'] >= 1000).sum()
	num_ge_2000 = (df['length'] >= 2000).sum()
	print("\nnumber (length >= 100bp): ", num_ge_100, file = o)
	print("number (length >= 200bp): ", num_ge_200, file = o)
	print("number (length >= 1000bp): ", num_ge_1000, file = o)
	print("number (length >= 2000bp): ", num_ge_2000, file = o)


def Count_N50(in_dict, outfile):
	o = open(outfile, "w")
	print("\tlength (bp)\tnumber", file = o)

	df = pd.DataFrame(data = in_dict).T
	df = df.sort_values(by = 0, ascending = False)
	df.columns = ['length', 'GC_length', 'N_length']

	tot_len = df['length'].sum()
	tot_GC = df['GC_length'].sum()
	tot_N = df['N_length'].sum()

	df['cum_sum'] = df['length'].cumsum()
	for key in range(10, 100, 10):
		label = 'n' + str(key)
		label = df[df['cum_sum'].gt(tot_len*key/100)].index[0]
		num = (df['cum_sum'] < (tot_len*key/100)).sum() + 1
		print("N%d\t%11d\t%6d" % (key, df['length'][label], num), file = o)

	num_ge_100 = (df['length'] >= 100).sum()
	num_ge_200 = (df['length'] >= 200).sum()
	num_ge_1000 = (df['length'] >= 1000).sum()
	num_ge_2000 = (df['length'] >= 2000).sum()
	print("\nnumber (length >= 100bp): %4d" % num_ge_100, file = o)
	print("number (length >= 200bp): %4d" % num_ge_200, file = o)
	print("number (length >= 1000bp): %5d" % num_ge_1000, file = o)
	print("number (length >= 2000bp): %5d" % num_ge_2000, file = o)

	max_len = df['length'].max()
	max_id = df.loc[df['length'].idxmax()].name
	print("\nmax_length: %d bp\t\tSeq_id: %s" % (max_len, max_id), file = o)
	fa_num = len(df.index)
	mean_len = round(tot_len / fa_num, 2)
	print("average_lenth: %6.2f bp" % (tot_len / fa_num), file = o)
	print("Total: %d bp\t\tGC%%: %2.2f%%\t\tN%%:%2.2f%%" % (tot_len, (100 * tot_GC / tot_len), (100 * tot_N / tot_len)), file = o)


def main():
	params = vars(get_parser())

	# print help if no parameters
	if len(sys.argv) == 1:
		params.print_help()
		sys.exit(1)

	genome = Read_fasta(params['in'])
	outfile = params['out']

	if params['mode'] == "N50":
		Count_N50(genome, outfile)

	if params['mode'] == "len":
		Count_len_gc_n(genome, outfile)


if __name__ == "__main__":
	main()
