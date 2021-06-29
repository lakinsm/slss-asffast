#!/usr/bin/env python3

import sys


def fasta_parse(infile):
	with open(infile, 'r') as fasta_file:
		# Skip whitespace
		while True:
			line = fasta_file.readline()
			if line == "":
				return  # Empty file or premature end of file?
			if line[0] == ">":
				break
		while True:
			if line[0] != ">":
				raise ValueError("Records in FASTA should begin with '>'")
			header = line[1:].rstrip()
			all_lines = []
			line = fasta_file.readline()
			while True:
				if not line:
					break
				if line[0] == ">":
					break
				all_lines.append(line.rstrip())
				line = fasta_file.readline()
			yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
			if not line:
				return  # Stop Iteration


def parse_final_observed_file(infile):
	ret = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')
	for line in data:
		if not line:
			continue
		entries = line.split()
		ret[entries[0]] = entries[1]
	return ret


if __name__ == '__main__':
	D = {k: v for k, v in fasta_parse(sys.argv[1])}
	best_genomes = parse_final_observed_file(sys.argv[2])
	barcode_id = sys.argv[3]

	selected_genome = best_genomes[barcode_id]
	sys.stdout.write('>{}\n{}\n'.format(
		selected_genome,
		D[selected_genome]
	))
