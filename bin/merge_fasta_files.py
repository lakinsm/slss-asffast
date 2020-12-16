#!/usr/bin/env python3

import sys
import glob
import os


def merge_fastas(glob_dir):
	for f in glob.glob(glob_dir + '/*.fasta'):
		barcode = f.split('/')[-1].split('_')[1]
		with open('{}_aligned_reads.fasta'.format(barcode), 'a') as out, open(f, 'r') as infile:
			data = infile.read()
			out.write(data)
		os.remove(f)


if __name__ == '__main__':
	merge_fastas(sys.argv[1])
