#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from samscore.samscore import SamParser
from samscore.samscore import parse_cigar


def plot_data(aln_data, aln_order, barcode, output_pdf_dir):
	colors = [plt.cm.Set1(i) for i in range(5)]

	# Alignment data dict(target: [ydat])
	x = np.array(range(len(aln_data[list(aln_data.keys())[0]])))
	plt.figure(figsize=(15, 10))
	axes = plt.gca()
	axes.set_xlim((x[0], x[-1]))
	axes.set_ylim(0, 100)
	axes.yaxis.set_major_locator(MaxNLocator(integer=True))
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.ticklabel_format(useOffset=False)

	aln_pdf = '{}/{}_alignment_timeseries_graph.pdf'.format(output_pdf_dir, barcode)
	with PdfPages(aln_pdf) as pdf_handle:
		counter = 0
		for target in aln_order:
			plt.plot(x, aln_data[target], marker='', color=colors[counter], label=target)
			counter += 1
		plt.legend()
		plt.xlabel('Minutes of Sequencing')
		plt.ylabel('Percent ASFV Genomic Coverage')
		plt.title('Percent ASFV Genome Coverage by Minutes of Sequencing, Top 5 Targets')
		pdf_handle.savefig()
		plt.close()


def write_timeseries(aln_data, barcode, output_csv_dir):
	target_order = [(v[-1], k) for k, v in aln_data.items()]
	target_order = [x for _, x in sorted(target_order, reverse=True)]
	with open('{}/{}_alignment_timeseries_data.csv'.format(output_csv_dir, barcode), 'w') as acsv_out:
		acsv_out.write('target,timepoint,percent_cov\n')
		for target in target_order:
			for i, cov in enumerate(aln_data[target]):
				acsv_out.write('{},{},{}\n'.format(
					target,
					i,
					cov
				))
	plot_data(aln_data, target_order, barcode, output_csv_dir)


def parse_seqfile(infile, sam_dict):
	with open(infile, 'r') as f:
		f.readline()
		line = f.readline()
		while line:
			entries = line.split()
			fq_name = entries[0]
			if '_fail_' not in fq_name:
				sectime = float(entries[6]) + float(entries[7])
				read_id = entries[2]
				try:
					sam_dict[read_id][0] = sectime
				except KeyError:
					pass
			line = f.readline()
	return sam_dict


def parse_sam(infile):
	ret = {}
	sam_parser = SamParser(infile)
	for header, _, thead, tstart, cigar in sam_parser:
		if header not in ret:
			# sequencing sectime, target_list, start_idx_list, cigar_list
			ret[header] = [None, (thead,), (tstart,), (cigar,)]
		else:
			ret[header][1] += (thead,)
			ret[header][2] += (tstart,)
			ret[header][3] += (cigar,)
	return ret, sam_parser.ref_len


def format_alignment_data(sam, ref):
	ret = {x: [0.] for x in ref.keys()}
	cov = {x: set() for x in ref.keys()}
	ref_sets = {k: set(range(v)) for k, v in ref.items()}
	cur_len = 1
	for sec, targets, starts, cigars in sam:
		idx_min = int(np.ceil(sec / 60.))
		if idx_min > cur_len:
			for target, covlist in ret.items():
				cov_val = 100. * float(len(cov[target].intersection(ref_sets[target]))) / float(ref[target])
				covlist += [cov_val for _ in range(idx_min - cur_len)]
			cur_len = idx_min
		for i in range(len(targets)):
			cov[targets[i]].update(parse_cigar(cigars[i], starts[i] - 1))
	last_cov = [(v[-1], k) for k, v in ret.items()]
	tops = set(j for _, j in sorted(last_cov, reverse=True)[:5])
	return {k: v for k, v in ret.items() if k in tops}


if __name__ == '__main__':
	sam_data, ref_data = parse_sam(sys.argv[1])
	sam_data = parse_seqfile(sys.argv[2], sam_data)
	aln_data = format_alignment_data(sorted(sam_data.values()), ref_data)
	barcode = sys.argv[1].split('/')[-1].split('_')[0]
	write_timeseries(aln_data, barcode, sys.argv[3])
