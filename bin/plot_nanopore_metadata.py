#!/usr/bin/env python3

import sys
import numpy as np
from queue import PriorityQueue
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages


def plot_data(fcount_data, thr_data, output_pdf_dir):
	# Filecount data
	y = np.array(fcount_data)
	x = np.array(range(len(fcount_data)))

	plt.figure(figsize=(15, 10))
	axes = plt.gca()
	axes.set_xlim((x[0], x[-1]))
	axes.set_ylim(y[0], y[-1])
	axes.yaxis.set_major_locator(MaxNLocator(integer=True))
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))

	fcount_pdf = output_pdf_dir + '/filecount_timeseries_graph.pdf'
	with PdfPages(fcount_pdf) as pdf_handle:
		plt.plot(x, y, marker='')
		plt.xlabel('Minutes of Sequencing')
		plt.ylabel('Number of Files Produced')
		plt.title('Files produced per minute of Nanopore sequencing time')
		pdf_handle.savefig()
		plt.close()

		# Filecount data
		zipped_thr = list(zip(*thr_data))
		bcr_pass = np.array(zipped_thr[1])
		bcr_fail = np.array(zipped_thr[2])
		x = np.array(range(len(bcr_pass)))
		colors = [plt.cm.Set1(i) for i in range(2)]

		plt.figure(figsize=(15, 10))
		axes = plt.gca()
		axes.set_xlim((x[0], x[-1]))
		axes.set_ylim(0, max(bcr_pass[-1], bcr_fail[-1]))
		axes.yaxis.set_major_locator(MaxNLocator(integer=True))
		axes.xaxis.set_major_locator(MaxNLocator(integer=True))
		axes.ticklabel_format(useOffset=False)

		tcount_pdf = output_pdf_dir + '/throughput_timeseries_graph.pdf'
		with PdfPages(tcount_pdf) as pdf_handle:
			plt.plot(x, bcr_pass, marker='', color=colors[0], label='Basecalled Reads Pass')
			plt.plot(x, bcr_fail, marker='', color=colors[1], label='Basecalled Reads Fail')
			plt.legend()
			plt.xlabel('Minutes of Sequencing')
			plt.ylabel('Number of Basecalled Reads Passing and Failing Filter')
			plt.title('Throughput of Nanopore read generation over minutes of sequencing time')
			pdf_handle.savefig()
			plt.close()


def write_timeseries(seq_dict, through_data, output_csv_dir):
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

	plot_data(filecounts, through_data, output_csv_dir)


def parse_seqfile(infile):
	ret = {}
	with open(infile, 'r') as f:
		f.readline()
		line = f.readline()
		while line:
			entries = line.split()
			fq_name = entries[0]
			if '_fail_' not in fq_name:
				sectime = float(entries[6])
				if fq_name not in ret:
					ret[fq_name] = sectime
				elif ret[fq_name] < sectime:
					ret[fq_name] = sectime
			line = f.readline()
	return ret


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


if __name__ == '__main__':
	seqfile_data = parse_seqfile(sys.argv[1])
	through_data = parse_throughfile(sys.argv[2])
	write_timeseries(seqfile_data, through_data, sys.argv[3])
