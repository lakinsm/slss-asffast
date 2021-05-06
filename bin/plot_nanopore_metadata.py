#!/usr/bin/env python3

import sys
import numpy as np
from queue import PriorityQueue
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages


def plot_data(fcount_data, thr_data, aln_data, output_pdf_dir):
	# Filecount data
	y = np.array(fcount_data)
	x = np.array(range(len(fcount_data)))

	plt.figure(figsize=(15, 10))
	axes = plt.gca()
	axes.set_xlim((x[0], x[-1]))
	axes.set_ylim(y[0], y[-1])
	axes.yaxis.set_major_locator(MaxNLocator(integer=True))
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	colors = [plt.cm.Set1(i) for i in range(5)]

	fcount_pdf = output_pdf_dir + '/filecount_timeseries_graph.pdf'
	with PdfPages(fcount_pdf) as pdf_handle:
		plt.plot(x, y, marker='')
		plt.xlabel('Minutes of Sequencing')
		plt.ylabel('Number of Files Produced')
		plt.title('Files produced per minute of Nanopore sequencing time')
		pdf_handle.savefig()
		plt.close()

	# Throughput data
	zipped_thr = list(zip(*thr_data))
	bcr_pass = np.array(zipped_thr[1])
	bcr_fail = np.array(zipped_thr[2])
	x = np.array(range(len(bcr_pass)))


	plt.figure(figsize=(15, 10))
	axes = plt.gca()
	axes.set_xlim((x[0], x[-1]))
	axes.set_ylim(0, max(bcr_pass[-1], bcr_fail[-1]))
	axes.yaxis.set_major_locator(MaxNLocator(integer=True))
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.ticklabel_format(useOffset=False)

	tcount_pdf = output_pdf_dir + '/throughput_timeseries_graph.pdf'
	with PdfPages(tcount_pdf) as pdf_handle2:
		plt.plot(x, bcr_pass, marker='', color=colors[0], label='Basecalled Reads Pass')
		plt.plot(x, bcr_fail, marker='', color=colors[1], label='Basecalled Reads Fail')
		plt.legend()
		plt.xlabel('Minutes of Sequencing')
		plt.ylabel('Number of Basecalled Reads Passing and Failing Filter')
		plt.title('Throughput of Nanopore read generation over minutes of sequencing time')
		pdf_handle2.savefig()
		plt.close()

	# Alignment data dict(target: [ydat])
	x = np.array(range(len(aln_data[list(aln_data.keys())[0]])))
	plt.figure(figsize=(15, 10))
	axes = plt.gca()
	axes.set_xlim((x[0], x[-1]))
	axes.set_ylim(0, 100)
	axes.yaxis.set_major_locator(MaxNLocator(integer=True))
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.ticklabel_format(useOffset=False)

	aln_pdf = output_pdf_dir + '/aligned_timeseries_graph.pdf'
	with PdfPages(aln_pdf) as pdf_handle3:
		counter = 0
		for target, ydat in aln_data.items():
			plt.plot(x, ydat, marker='', color=colors[counter], label=target)
			counter += 1
		plt.legend()
		plt.xlabel('Minutes of Sequencing')
		plt.ylabel('Percent ASFV Genomic Coverage')
		plt.title('Percent ASFV Genome Coverage by Minutes of Sequencing, Top 5 Targets')
		pdf_handle3.savefig()
		plt.close()


def write_timeseries(seq_dict, through_data, aln_data, output_csv_dir):
	q = PriorityQueue(maxsize=-1)
	for v in seq_dict.values():
		q.put(v, block=False)
	filecounts = []
	cur_sectime = 0.
	qval = q.get(block=False)
	while not q.empty():
		cur_sectime += 60.
		try:
			filecounts.append(filecounts[-1])
		except IndexError:
			filecounts.append(0)
		while qval < cur_sectime:
			filecounts[-1] += 1
			if q.empty():
				break
			qval = q.get(block=False)
	with open(output_csv_dir + '/nanopore_filecounts.csv', 'w') as fcsv_out:
		fcsv_out.write('time_minutes,file_count\n')
		for i, val in enumerate(filecounts):
			fcsv_out.write('{},{}\n'.format(
				i+1,
				val
			))
	with open(output_csv_dir + '/nanopore_throughput.csv', 'w') as tcsv_out:
		tcsv_out.write('time_minutes,basecalled_reads_pass,basecalled_reads_fail,basecalled_basepairs\n')
		for mintime, bcr_pass, bcr_fail, bcbp in through_data:
			tcsv_out.write('{},{},{},{}\n'.format(
				mintime,
				bcr_pass,
				bcr_fail,
				bcbp
			))

	with open(output_csv_dir + '/nanopore_stats_overall.txt', 'w') as out:
		out.write('Total Minutes: {}\nTotal Basecalled Reads Pass Filter: {}\nTotal Basecalled Reads Fail Filter: {}'
		          '\nPercent Basecalled Reads Pass Filter: {}%\nTotal Bases Called: {}\n'.format(
			through_data[-1][0],
			through_data[-1][1],
			through_data[-1][2],
			100 * float(through_data[-1][1]) / (float(through_data[-1][1]) + float(through_data[-1][2])),
			through_data[-1][3]
		))

	with open(output_csv_dir + '/alignment_timeseries.csv', 'w') as acsv_out:
		acsv_out.write('target,timepoint,percent_cov\n')
		for target, covdat in aln_data.items():
			for i, cov in enumerate(covdat):
				acsv_out.write('{},{},{}\n'.format(
					target,
					i,
					cov
				))

	plot_data(filecounts, through_data, aln_data, output_csv_dir)


def parse_seqfile(infile, sam_dict):
	ret = {}
	with open(infile, 'r') as f:
		f.readline()
		line = f.readline()
		while line:
			entries = line.split()
			fq_name = entries[0]
			if '_fail_' not in fq_name:
				sectime = float(entries[6]) + float(entries[7])
				if fq_name not in ret:
					ret[fq_name] = sectime
				elif ret[fq_name] < sectime:
					ret[fq_name] = sectime
				read_id = entries[2]
				try:
					sam_dict[read_id][0] = sectime
				except KeyError:
					pass
			line = f.readline()
	return ret, sam_dict


def parse_throughfile(infile):
	ret = tuple()
	with open(infile, 'r') as f:
		f.readline()
		line = f.readline()
		while line:
			entries = line.split(',')
			mintime = int(entries[0])
			bcr_pass = int(entries[2])
			bcr_fail = int(entries[3])
			bcbp = int(entries[8])
			ret += ((mintime, bcr_pass, bcr_fail, bcbp),)
			line = f.readline()
	return ret


class SamParser(object):
	"""
	Line-by-line parsing of Sequence Alignment Map (SAM) UTF-8 encoded files.  Stores alignment information for
	genome lengths from the provided reference database using the SAM headers, if specified.  Outputs
	aligned reads if specified and only if the aligned reads do not match a provided filter file (headers of reference
	to exclude).  Aligning reads, if output, are output as they are parsed.
	"""

	def __init__(self, sam_path, capture_ref_len=True):
		"""
		Initialize object and data structures for storing coverage and coverage over time, if specified.
		:param sam_path: STR, path to input SAM file
		:param capture_ref_len: BOOL, store reference genome lengths from SAM headers
		"""
		self.sam_path = sam_path
		self.ref_len = {}
		self.output_handle = None
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
		target_header = entries[2]
		target_start = int(entries[3])
		cigar = entries[5]
		if sam_flag & 16 != 0:  # if reverse
			reverse = True
		else:
			reverse = False
		self.line = self.handle.readline()
		return query_header, reverse, target_header, target_start, cigar

	def _open(self):
		self.handle = open(self.sam_path, 'r')

	def _close(self):
		self.handle.close()


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


def parse_sam(infile):
	ret = {}
	sam_parser = SamParser(infile)
	for header, rev, thead, tstart, cigar in sam_parser:
		if header not in ret:
			# sequencing sectime, target_list, start_idx_list, rev_list, cigar_list
			ret[header] = [None, (thead,), (tstart,), (rev,), (cigar,)]
		else:
			ret[header][1] += (thead,)
			ret[header][2] += (tstart,)
			ret[header][3] += (rev,)
			ret[header][4] += (cigar,)
	return ret, sam_parser.ref_len


def format_alignment_data(sam, ref):
	ret = {x: [0.] for x in ref.keys()}
	cov = {x: set() for x in ref.keys()}
	ref_sets = {k: set(range(v)) for k, v in ref.items()}
	cur_len = 1
	for sec, targets, starts, revs, cigars in sam:
		idx_min = int(np.ceil(sec / 60.))
		if idx_min > cur_len:
			for target, covlist in ret.items():
				cov_val = 100. * float(len(cov[target].intersection(ref_sets[target]))) / float(ref[target])
				covlist += [cov_val for _ in range(idx_min - cur_len)]
			cur_len = idx_min
		for i in range(len(targets)):
			cov[targets[i]].update(parse_cigar(cigars[i], starts[i] - 1, revs[i]))
	last_cov = [(v[-1], k) for k, v in ret.items()]
	print(last_cov)
	for k, v in cov.items():
		print(k, len(v), ref[k])
	tops = set(j for _, j in sorted(last_cov, reverse=True)[:5])
	return {k: v for k, v in ret.items() if k in tops}


if __name__ == '__main__':
	sam_data, ref_lens = parse_sam(sys.argv[1])
	seqfile_data, sam_data = parse_seqfile(sys.argv[2], sam_data)
	through_data = parse_throughfile(sys.argv[3])
	aln_data = format_alignment_data(sorted(sam_data.values()), ref_lens)
	write_timeseries(seqfile_data, through_data, aln_data, sys.argv[4])
