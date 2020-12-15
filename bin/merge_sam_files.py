#!/usr/bin/env python3

import sys
import ast


def parse_sam(infile, headers=False):
	"""
	Screen input Sequence Alignment Map (SAM) file and output to stdout to merge all input SAM files.
	:param infile: STR, file path to input SAM file
	:param headers: BOOL, include header information in output
	:return: None
	"""
	with open(infile, 'r') as f:
		line = f.readline()
		while line:
			if line[0] == '@' and not headers:
				line = f.readline()
				continue
			sys.stdout.write(line)
			line = f.readline()


if __name__ == '__main__':
	# Args: space-separated list of SAM files to consider
	headers_written = False

	for sam_file in sys.argv[1:]:
		sam_file = sam_file.lstrip('[').rstrip(' ').rstrip(',').rstrip(']')
		if headers_written:
			parse_sam(sam_file)
		else:
			parse_sam(sam_file, headers=True)
			headers_written = True
