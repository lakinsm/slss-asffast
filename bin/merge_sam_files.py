#!/usr/bin/env python3

import sys


def parse_sam(infile, outfile, headers=False):
	"""
	Screen input Sequence Alignment Map (SAM) file and output to stdout to merge all input SAM files.
	:param infile: STR, file path to input SAM file
	:param outfile: STR, file path to output SAM file in append mode
	:param headers: BOOL, include header information in output
	:return: None
	"""
	with open(infile, 'r') as f, open(outfile, 'a') as out:
		line = f.readline()
		while line:
			if line[0] == '@':
				if headers:
					out.write(line)
			else:
				entries = line.split('\t')
				if int(entries[1]) & 4 == 0:
					out.write(line)
			line = f.readline()


if __name__ == '__main__':
	# Args: space-separated list of SAM files to consider
	headers_written = set()

	with open(sys.argv[1], 'r') as sfile:
		data = sfile.read().rstrip('\n')
	sam_file_list = data.lstrip("\'").rstrip("\'").lstrip('[').rstrip(']').split(', ')

	for sam_file in sam_file_list:
		barcode = sam_file.split('/')[-1].split('_')[0]
		out_file_path = '{}_aligned_reads.sam'.format(barcode)
		if barcode not in headers_written:
			parse_sam(sam_file, out_file_path, headers=True)
			headers_written.add(barcode)
		else:
			parse_sam(sam_file, out_file_path)
