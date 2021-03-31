#!/usr/bin/env python3


import sys
import argparse
import queue
from queue import PriorityQueue
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import multiprocessing as mp


class SamParser(object):
	"""
	Line-by-line parsing of Sequence Alignment Map (SAM) UTF-8 encoded files.  Stores alignment information for
	genome lengths from the provided reference database using the SAM headers, if specified.  Outputs
	aligned reads if specified and only if the aligned reads do not match a provided filter file (headers of reference
	to exclude).  Aligning reads, if output, are output as they are parsed.
	"""

	def __init__(self, sam_path, screen_headers=None, output_reads=None, capture_ref_len=True):
		"""
		Initialize object and data structures for storing coverage and coverage over time, if specified.
		:param sam_path: STR, path to input SAM file
		:param screen_headers: LIST, reference headers to ignore if encountered in SAM file alignments
		:param output_reads: STR, output aligned reads passing filter/screen to this filepath in FASTA format
		:param capture_ref_len: BOOL, store reference genome lengths from SAM headers
		"""
		self.sam_path = sam_path
		self.ref_len = {}
		if screen_headers:
			self.filter = set(screen_headers)
		self.output_reads = output_reads
		self.output_handle = None
		if output_reads:
			self.seen_reads = set()
		self.handle = None

		# Read SAM header information to retrieve database headers and sequence lengths
		self._open()
		self.line = self.handle.readline()
		if self.line:
			while self.line[0] == '@':
				if capture_ref_len and self.line.startswith('@SQ'):
					entries = self.line.split('\t')
					db_header = entries[1].split(':')[-1]
					db_seq_len = int(entries[2].split(':')[-1])
					self.ref_len[db_header] = db_seq_len
				self.line = self.handle.readline()

	def __iter__(self):
		return self

	def __next__(self):
		"""
		Iterator for SAM file handle, yields aligned read entries. The start_idx yielded from the SAM file is
		one-indexed and must have 1 subtracted from it to translate it into zero-indexing.
		:yield: (query, query_reverse_bool, target, target_start_idx, CIGAR)
		"""
		if not self.line:
			self._close()
			raise StopIteration
		entries = self.line.split('\t')
		sam_flag = int(entries[1])
		while sam_flag & 4 != 0:
			self.line = self.handle.readline()
			if not self.line:
				self._close()
				raise StopIteration
			entries = self.line.split('\t')
			sam_flag = int(entries[1])
		query_header = entries[0]
		query_seq = entries[9]
		target_header = entries[2]
		target_start = int(entries[3])
		cigar = entries[5]
		if sam_flag & 16 != 0:  # if reverse
			reverse = True
		else:
			reverse = False
		if self.output_reads:
			if query_header not in self.seen_reads:
				self.seen_reads.add(query_header)
				self.output_handle.write('>{}\n{}\n'.format(
					query_header,
					query_seq
				))
		self.line = self.handle.readline()
		return query_header, reverse, target_header, target_start, cigar

	def _open(self):
		self.handle = open(self.sam_path, 'r')
		if self.output_reads:
			self.output_handle = open(self.output_reads, 'w')

	def _close(self):
		self.handle.close()
		if self.output_reads:
			self.output_handle.close()


def parse_cigar(s, t_idx, rev):
	"""
	Parse SAM CIGAR alignment string and return indices to which the read aligned.
	:param s: STR, CIGAR string
	:param t_idx: INT, zero-based index for target start position
	:param rev: BOOL, read aligned in reverse orientation if true
	:return: tuple of integers, zero-indexed indices to which the read aligned
	"""
	ret = ()
	num = ''
	c_idx = 0
	while c_idx < len(s):
		if s[c_idx].isdigit():
			num += s[c_idx]
		else:
			op = s[c_idx]
			if rev:
				if op == 'M' or op == '=':
					ret += tuple(range(t_idx - int(num) + 1, t_idx + 1))
					t_idx -= int(num)
				elif op == 'D':
					t_idx -= int(num)
				elif op == 'N':
					t_idx -= int(num)
				elif op == 'X':
					ret += tuple(range(t_idx - int(num) + 1, t_idx + 1))
					t_idx -= int(num)
			else:
				if op == 'M' or op == '=':
					ret += tuple(range(t_idx, t_idx + int(num)))
					t_idx += int(num)
				elif op == 'D':
					t_idx += int(num)
				elif op == 'N':
					t_idx += int(num)
				elif op == 'X':
					ret += tuple(range(t_idx, t_idx + int(num)))
					t_idx += int(num)
			num = ''
		c_idx += 1
	return ret


def find_top_targets(cov_dict, n=5):
	"""
	Using a min priority queue, store top 5 coverage results and return a list of genome names to target.
	:param cov_dict: DICT, dictionary of coverage over indices from SamParser {target: {idx: cov}}
	:param n: INT, number of top coverage genomes to return
	:return: tuple of genome names (target keys) in order of coverage, descending
	"""
	ret = ()
	tops = PriorityQueue()
	for k, v in cov_dict.items():
		tops.put((-sum(v.values()), k))
	for _ in range(n):
		try:
			ret += (tops.get(block=False)[1],)
		except queue.Empty:
			break
	return ret


def rolling_average(data, window_size, interval, max_len):
	"""
	Generate a vector of averages for each window of size window_size at fixed intervals interval.
	:param data: TUPLE, (x, y) pairs of integers, where x is the genome idx and y is the depth over that idx
	:param window_size: INT, size of each window over which to compute the averages
	:param interval: INT, fixed interval to shift the window
	:param max_len: INT, maximum index for the end window
	:return: tuple of rolling averages
	"""
	ret = ()
	local_data = np.array(data)
	local_data.sort(axis=0)
	if window_size > max_len:
		window_size = np.ceil(max_len / 50)
		interval = np.ceil(window_size / 2)
	window_start = 0
	window_end = window_size
	while window_end <= max_len:
		idxs = np.where((local_data[:, 0] < window_end) * (local_data[:, 0] >= window_start))
		ret += (np.sum(local_data[idxs, 1]) / float(window_size),)
		window_start += interval
		window_end += interval
	return ret


def plot_cov(cov_dict, target, target_len, window_size, pdf_handle):
	"""
	Using matplotlib.pyplot, produce coverage plots for the top n targets and output to PDF, one per page.
	:param cov_dict: DICT, dictionary of coverage over indices from SamParser {target: {idx: cov}}
	:param target: STR, reference genome name to plot, must be present in cov_dict
	:param target_len: INT, reference genome length of target
	:param window_size: INT, window size over which to produce rolling averages
	:param pdf_handle: FILE HANDLE, pdf file handle object for saving figures
	:return: None
	"""
	y = rolling_average([(k, v) for k, v in cov_dict[target].items()], window_size, np.ceil(window_size / 2),
	                    target_len)
	x = np.ceil(window_size / 2) * np.array(range(len(y)))
	plt.figure(figsize=(15, 10))
	plt.fill_between(x, y, 0,
	                 facecolor='blue',
	                 color='black',
	                 alpha=0.3)
	plt.xlabel('Genomic Location')
	plt.ylabel('Coverage (Rolling Average)')
	plt.title('Rolling Average Coverage, Reference Sequence {}'.format(target))
	plt.ticklabel_format(style='sci', axis='x')
	pdf_handle.savefig()
	plt.close()


def plot_timeseries(cov_dict, output_csv_path, n=5):
	"""
	Plot a multi-line matplotlib pyplot of genome coverage over time for the top n targets at the final timepoint.
	:param cov_dict: DICT, dictionary of coverage over indices from SamParser {target: {idx: cov}}
	:param output_csv_path: STR, filepath of the output .csv file, which will be modified to .pdf
	:param n: INT, number of targets to plot with highest genomic coverage at the final timepoint
	:return: None
	"""
	pdf_path = output_csv_path.split('.')[0] + '.pdf'
	ordered_timepoints = sorted(cov_dict.keys())

	# Extract top n targets at final timepoint
	tops = PriorityQueue()
	top_targets = ()
	for target, cov in cov_dict[ordered_timepoints[-1]].items():
		tops.put((-cov, target))
	for _ in range(n):
		try:
			top_targets += (tops.get(block=False)[1],)
		except queue.Empty:
			break

	# Generate color palette for n colors and setup figure
	colors = [plt.cm.Set1(i) for i in range(n)]
	plt.figure(figsize=(15, 10))
	axes = plt.gca()
	axes.set_xlim((0, ordered_timepoints[-1]))
	axes.set_ylim((0, 100))
	axes.yaxis.set_major_locator(MaxNLocator(integer=True))
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))

	# Extract timepoint data for each target in top_targets
	x = [0] + ordered_timepoints
	analytic_data = {}
	with PdfPages(pdf_path) as pdf_handle:
		for i, t in enumerate(top_targets):
			analytic_data[t] = [0.]
			for timepoint in ordered_timepoints:
				try:
					analytic_data[t].append(cov_dict[timepoint][t])
				except KeyError:
					analytic_data[t].append(0.)

			plt.plot(x, analytic_data[t], marker='', color=colors[i], label=t)
		plt.legend()
		plt.xlabel('Time Point')
		plt.ylabel('Percent of Genome Covered (>= 1x)')
		plt.title('Cumulative Coverage Over Time for Top {} Genomic Targets'.format(n))
		pdf_handle.savefig()
		plt.close()


def write_timeseries(cov_dict, ref_len_dict, output_csv_path):
	"""
	Write a long-format .csv file containing genomic coverage percentages at each timepoint for all targets.
	:param cov_dict: DICT, dictionary of coverage over indices from SamParser {target: {idx: cov}}
	:param ref_len_dict: DICT, dictionary of reference genome lengths {target: length}
	:param output_csv_path: STR, filepath of the output .csv file
	:return: None
	"""
	all_targets = set()
	for target in [vk for _, v in cov_dict.items() for vk, _ in v.items()]:
		all_targets.add(target)

	with open(output_csv_path, 'w') as out:
		out.write('timepoint,target,percent_coverage,genome_length\n')
		for timepoint, data in cov_dict.items():
			for target, cov in data.items():
				try:
					out.write('{},{},{},{}\n'.format(
						timepoint,
						target,
						cov,
						ref_len_dict[target]
					))
				except KeyError:
					out.write('{},{},0,{}\n'.format(
						timepoint,
						target,
						ref_len_dict[target]
					))


def write_coverage(cov_dict, ref_len_dict, barcode_to_sample_dict, output_csv_path, final=False):
	"""
	Write a long-format .csv file containing genomic coverage percentages at each timepoint for all targets.
	:param cov_dict: DICT, dictionary of coverage over indices from worker {barcode: {target: set(idxs)}
	:param ref_len_dict: DICT, dictionary of reference genome lengths {target: length}
	:param output_csv_path: STR, filepath of the output .csv file
	:param final: BOOL, if this is the final run, output avg coverage in addition to genomic coverage percent
	:return: None
	"""
	target_info = {}
	if final:
		for barcode, tdict in cov_dict.items():
			if barcode not in target_info:
				target_info[barcode] = ()
			for target, idx_dict in tdict.items():
				genome_len = ref_len_dict[target]
				bases_aligned = sum(idx_dict.values())
				bases_covered = len(set(range(genome_len)).intersection(idx_dict.keys()))
				target_info[barcode] += ((bases_aligned, bases_covered, genome_len, target),)
		with open(output_csv_path, 'w') as out:
			out.write('barcode,sample_id,target,genome_length,avg_coverage,percent_bases_covered\n')
			for barcode, data in target_info.items():
				sorted_targets = sorted(data, key=lambda x: x[0], reverse=True)
				for bases_aligned, bases_covered, genome_len, target in sorted_targets:
					out.write('{},{},{},{},{},{}\n'.format(
						barcode,
						barcode_to_sample_dict[barcode],
						target,
						genome_len,
						float(bases_aligned) / float(genome_len),
						100 * float(bases_covered) / float(genome_len)
					))
	else:
		for barcode, idx_dict in cov_dict.items():
			if barcode not in target_info:
				target_info[barcode] = ()
			for target, idxs in idx_dict.items():
				genome_len = ref_len_dict[target]
				bases_covered = len(set(range(genome_len)).intersection(idxs))
				target_info[barcode] += ((100 * float(bases_covered) / float(genome_len), genome_len, target),)
		with open(output_csv_path, 'w') as out:
			out.write('barcode,sample_id,target,genome_length,percent_bases_covered\n')
			for barcode, data in target_info.items():
				sorted_targets = sorted(data, key=lambda x: x[0], reverse=True)
				for bases_covered, genome_len, target in sorted_targets:
					out.write('{},{},{},{},{}\n'.format(
						barcode,
						barcode_to_sample_dict[barcode],
						target,
						genome_len,
						bases_covered
					))


def worker(infile):
	ref_cov = {}
	barcode_id = infile.split('/')[-1].split('_')[0]
	sample_id = infile.split('/')[-1].replace(barcode_id + '_', '').replace('.sam', '')
	sam_parser = SamParser(infile)
	for query, q_reverse, target, t_start, cigar in sam_parser:
		idxs = parse_cigar(cigar, t_start - 1, q_reverse)
		if target not in ref_cov:
			ref_cov[target] = set()
		for idx in idxs:
			ref_cov[target].add(idx)
	return barcode_id, sample_id, ref_cov


def final_worker(infile):
	ref_cov = {}
	barcode_id = infile.split('/')[-1].split('_')[0]
	sample_id = infile.split('/')[-1].replace(barcode_id + '_', '').replace('.sam', '')
	timepoint = infile.split('/')[-1].replace('.sam', '').split('_')[-1]
	sam_parser = SamParser(infile, output_reads='{}_{}_aligned_reads.fasta'.format(timepoint, barcode_id))
	for query, q_reverse, target, t_start, cigar in sam_parser:
		idxs = parse_cigar(cigar, t_start - 1, q_reverse)
		if target not in ref_cov:
			ref_cov[target] = {}
		for i in idxs:
			if i in ref_cov[target]:
				ref_cov[target][i] += 1
			else:
				ref_cov[target][i] = 1
	return barcode_id, sample_id, int(timepoint), ref_cov


parser = argparse.ArgumentParser('sam_parser.py')
parser.add_argument('-i', '--input', type=str, default=None, required=True,
                    help='Input SAM filepath (single file, or an array as a string if -ot enabled)')
parser.add_argument('-s', '--screen', type=str, default=None, help='Text file of reference database headers to exclude'
                                                                   'from SAM alignments, one per line')
parser.add_argument('-r', '--read_output', action='store_true', default=False, help='Output aligned reads to stdout')
parser.add_argument('-oc', '--output_cov', type=str, default=None, help='Filepath for primary coverage output file')
parser.add_argument('-op', '--output_pdf', type=str, default=None, help='Filepath for pdf genome coverage output file')
parser.add_argument('-ot', '--output_time', type=str, default=None,
                    help='Filepath for time series .csv output file, overrides/exclusive with -oc and -op')
parser.add_argument('-tmax', '--time_max', type=int, default=400, help='Maximum time point to analyze for time series')
parser.add_argument('--final', action='store_true', default=False,
                    help='Final run (not an intermediate): output all data as opposed to minimal')
parser.add_argument('--threads', type=int, default=1, help='Number of threads to utilize in parallel')

if __name__ == '__main__':
	mp.freeze_support()
	args = parser.parse_args()

	# Setup data structures
	barcode_to_ref_cov = {}
	barcode_to_sample_id = {}

	# Format input list
	with open(args.input, 'r') as sfile:
		data = sfile.read().rstrip('\n')
	sam_file_list = data.lstrip("\'").rstrip("\'").lstrip('[').rstrip(']').split(', ')

	if not args.final:
		pool = mp.Pool(processes=args.threads)
		res = pool.map(worker, sam_file_list)
		pool.close()
		pool.join()
		pool.terminate()
		del pool

		# Get initial reference lengths
		this_sam_parser = SamParser(sam_file_list[0])

		# Combine results
		for barcode, sample, cov_dict in res:
			if barcode not in barcode_to_ref_cov:
				barcode_to_ref_cov[barcode] = cov_dict
				barcode_to_sample_id[barcode] = sample.split('_')[0]
			else:
				for target, idxs in cov_dict.items():
					if target in barcode_to_ref_cov[barcode]:
						for idx in idxs:
							barcode_to_ref_cov[barcode][target].add(idx)
					else:
						barcode_to_ref_cov[barcode][target] = idxs

		write_coverage(barcode_to_ref_cov, this_sam_parser.ref_len, barcode_to_sample_id, args.output_cov)
	else:

		pool = mp.Pool(processes=args.threads)
		res = pool.map(final_worker, sam_file_list)
		pool.close()
		pool.join()
		pool.terminate()
		del pool

		# Get initial reference lengths
		this_sam_parser = SamParser(sam_file_list[0])

		# Dictionaries for time series reference database coverage and for storing observed idxs cumulatively
		timeseries_cov = {}
		timepoint_observed_idxs = {}  # target: set(idxs)

		# Combine timeseries results
		for barcode, sample, timepoint, cov_dict in res:
			if timepoint <= args.time_max:
				if barcode not in timepoint_observed_idxs:
					timepoint_observed_idxs[barcode] = {timepoint: {}}
					timeseries_cov[barcode] = {}
				if barcode not in barcode_to_sample_id:
					barcode_to_sample_id[barcode] = sample.split('_')[0]
				for target, idx_dict in cov_dict.items():
					if timepoint not in timepoint_observed_idxs[barcode]:
						timepoint_observed_idxs[barcode][timepoint] = {}
					timepoint_observed_idxs[barcode][timepoint][target] = set(idx_dict.keys())

		# Finalize timeseries data
		for barcode, tsdict in timepoint_observed_idxs.items():
			cumulative = {k: set() for k in this_sam_parser.ref_len.keys()}
			for timepoint in range(args.time_max):
				timeseries_cov.setdefault(barcode, {timepoint: {}})
				if timepoint not in timeseries_cov[barcode]:
					timeseries_cov[barcode][timepoint] = {}
				for target in cumulative.keys():
					if timepoint in tsdict:
						if target in tsdict[timepoint]:
							for idx in tsdict[timepoint][target]:
								cumulative[target].add(idx)
					timeseries_cov[barcode][timepoint][target] = np.min((100., 100 * float(len(cumulative[target])) / \
					                                                     float(this_sam_parser.ref_len[target])))

		for barcode, tsdict in timeseries_cov.items():
			csv_path = '{}_{}_timeseries_results.csv'.format(barcode, barcode_to_sample_id[barcode])
			plot_timeseries(tsdict, csv_path)
			write_timeseries(tsdict, this_sam_parser.ref_len, csv_path)

		# Combine coverage results
		for barcode, sample, timepoint, cov_dict in res:
			if barcode not in barcode_to_ref_cov:
				barcode_to_ref_cov[barcode] = cov_dict
			if barcode not in barcode_to_sample_id:
				barcode_to_sample_id[barcode] = sample.split('_')[0]
			else:
				for target, idx_dict in cov_dict.items():
					if target not in barcode_to_ref_cov[barcode]:
						barcode_to_ref_cov[barcode][target] = idx_dict
					else:
						for idx, count in idx_dict.items():
							if idx in barcode_to_ref_cov[barcode][target]:
								barcode_to_ref_cov[barcode][target][idx] += count
							else:
								barcode_to_ref_cov[barcode][target][idx] = count

		# Finalize coverage data
		write_coverage(barcode_to_ref_cov, this_sam_parser.ref_len, barcode_to_sample_id, args.output_cov, final=True)

		for barcode, cov_dict in barcode_to_ref_cov.items():
			top_targets = find_top_targets(cov_dict)
			n_windows = 500  # Number of points to graph rolling avg for coverage plots
			sample_id = barcode_to_sample_id[barcode]
			pdf_path = '{}_{}_coverage_results.pdf'.format(barcode, sample_id)
			with PdfPages(pdf_path) as pdf:
				for t in top_targets:
					# k = 2n / (m + 1)
					window_size = np.ceil(float(2 * this_sam_parser.ref_len[t]) / float(n_windows + 1))
					window_size = np.max((window_size, 1))
					plot_cov(cov_dict, t, this_sam_parser.ref_len[t], window_size, pdf)
