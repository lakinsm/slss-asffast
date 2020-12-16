#!/usr/bin/env python3

import sys


def fastq_to_fasta(infile):
	out_data = {}
	line_buffer = ''
	n_lines = 0
	with open(infile, 'r') as f:
		line = f.readline()
		while line:
			header = line[1:].rstrip('\n')
			line = f.readline()
			if not line:
				break
			while line[0] != '+':
				line_buffer += line.rstrip('\n')
				n_lines += 1
				line = f.readline()
				if not line:
					break
			if not line:
				break
			for _ in range(n_lines+1):
				line = f.readline()
			if not line:
				break
			out_data[header] = line_buffer
			line_buffer = ''
			n_lines = 0

	for header, seq in out_data.items():
		sys.stdout.write('>{}\n{}\n'.format(
			header,
			seq
		))


if __name__ == '__main__':
	fastq_to_fasta(sys.argv[1])
