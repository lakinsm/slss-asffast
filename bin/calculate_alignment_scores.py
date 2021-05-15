#!/usr/bin/env python3

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from samscore.samscore import SamParser
from samscore.samscore import ReadScoreCache
from samscore.samscore import score_cigar


def insert_read_scores(target_info_list, read_score_obj, zstart):
	# idx_scores
	for i, v in enumerate(read_score_obj[1]):
		try:
			target_info_list[0][zstart + i] += v
		except KeyError:
			target_info_list[0][zstart + i] = v

	# match_idxs
	for i in read_score_obj[2]:
		try:
			target_info_list[1][i] += 1
		except KeyError:
			target_info_list[1][i] = 1

	# mismatch_idxs
	for i in read_score_obj[3]:
		try:
			target_info_list[2][i] += 1
		except KeyError:
			target_info_list[2][i] = 1

	# insert_idxs
	for i in read_score_obj[4]:
		try:
			target_info_list[3][i] += 1
		except KeyError:
			target_info_list[3][i] = 1

	# delete_idxs
	for i in read_score_obj[5]:
		try:
			target_info_list[4][i] += 1
		except KeyError:
			target_info_list[4][i] = 1


if __name__ == '__main__':
	sam_parser = SamParser(sys.argv[1])
	read_cache = ReadScoreCache()
	target_info = {}
	for target in sam_parser.ref_len.keys():
		target_info[target] = [{} for _ in range(5)]

	for qheader, _, theader, tstart, cigar, _ in sam_parser:
		# read_score = tuple(score, [idx_scores], (match_idxs,), (mismatch_idxs,), (insert_idxs,), (delete_idxs,))
		read_score = score_cigar(cigar, tstart - 1)
		top_read = read_cache.smart_insert(qheader, read_score)
		if top_read:
			insert_read_scores(target_info[theader], top_read, tstart - 1)

	# final entry
	final_read = read_cache.finalize()
	insert_read_scores(sam_parser.values()[0], final_read, sam_parser.values()[3] - 1)

	for k, v in target_info.items():
		sys.stdout.write('{}\t{} ({})\t{} ({})\t{} ({})\t{} ({})\t{} ({})'.format(
			k,
			len(v[0]),
			sum(v[0].values()),
			len(v[1]),
			sum(v[1].values()),
			len(v[2]),
			sum(v[2].values()),
			len(v[3]),
			sum(v[3].values()),
			len(v[4]),
			sum(v[4].values()),
		))
