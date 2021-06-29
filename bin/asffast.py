#!/usr/bin/env python3

import sys
import os
import re
import subprocess
import math
import argparse
import glob
import time
from queue import PriorityQueue
from datetime import datetime


SINGULARITY_LIB_PATH = 'library://lakinsm/default/slss-asffast:alpha1'
SLURM_CONFIG = 'slss_hpc_slurm.config'
NEXTFLOW_MAX_FORKS = 16  # Maximum file parallelism to execute (should match maxForks in nextflow.config)
BARCODE_REGEX = re.compile(r'/(barcode[0-9]{1,2}[a-z]?)/')  # Regular expression to detect if barcode dirs are present
BAR_LENGTH = 35
TERMINAL_LENGTH = 80


def available_cpu_count():
	"""
	Determine the number of available CPUs on the system in a platform-independent manner.
	:return: INT, number of virtual cores/threads not in use
	"""

	# cpuset
	# cpuset may restrict the number of *available* processors
	try:
		m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
					  open('/proc/self/status').read())
		if m:
			res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
			if res > 0:
				return res
	except IOError:
		pass

	# Python 2.6+
	try:
		import multiprocessing
		return multiprocessing.cpu_count()
	except (ImportError, NotImplementedError):
		pass

	# https://github.com/giampaolo/psutil
	try:
		import psutil
		return psutil.cpu_count()   # psutil.NUM_CPUS on old versions
	except (ImportError, AttributeError):
		pass

	# POSIX
	try:
		res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

		if res > 0:
			return res
	except (AttributeError, ValueError):
		pass

	# Windows
	try:
		res = int(os.environ['NUMBER_OF_PROCESSORS'])

		if res > 0:
			return res
	except (KeyError, ValueError):
		pass

	# jython
	try:
		from java.lang import Runtime
		runtime = Runtime.getRuntime()
		res = runtime.availableProcessors()
		if res > 0:
			return res
	except ImportError:
		pass

	# BSD
	try:
		sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
								  stdout=subprocess.PIPE)
		scStdout = sysctl.communicate()[0]
		res = int(scStdout)

		if res > 0:
			return res
	except (OSError, ValueError):
		pass

	# Linux
	try:
		res = open('/proc/cpuinfo').read().count('processor\t:')

		if res > 0:
			return res
	except IOError:
		pass

	# Solaris
	try:
		pseudoDevices = os.listdir('/devices/pseudo/')
		res = 0
		for pd in pseudoDevices:
			if re.match(r'^cpuid@[0-9]+$', pd):
				res += 1

		if res > 0:
			return res
	except OSError:
		pass

	# Other UNIXes (heuristic)
	try:
		try:
			dmesg = open('/var/run/dmesg.boot').read()
		except IOError:
			dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
			dmesg = dmesgProcess.communicate()[0]

		res = 0
		while '\ncpu' + str(res) + ':' in dmesg:
			res += 1

		if res > 0:
			return res
	except OSError:
		pass

	raise Exception('Can not determine number of CPUs on this system')


def get_real_dir(infile):
	return os.path.dirname(os.path.realpath(infile))


def get_create_real_dir(inpath):
	real_work_dir = os.path.realpath(inpath)
	dirname = inpath.rstrip('/').split('/')[-1]
	if os.path.isdir(real_work_dir):
		return real_work_dir
	else:
		parent_dir = '/'.join(inpath.rstrip('/').split('/')[:-1])
		real_parent_dir = os.path.realpath(parent_dir)
		if not os.path.isdir(real_parent_dir):
			sys.stderr.write('Error: {} is not a directory. For the working directory, at least '
			                 'its parent directory must exist.\n')
			raise ValueError
		real_dirpath = real_parent_dir + '/' + dirname
		os.mkdir(real_dirpath)
	return real_dirpath


def report_intermediate_coverage(infile, n=3):
	intermediate_best_genomes = {}
	end_col = '\u001b[0m'
	sys.stdout.write('\n{:%Y-%m-%d %H:%M:%S} Coverage Results for Top {} Genomes Per Barcode:\n\n'.format(
		datetime.now(),
		n
	))
	data_dict = {}
	with open(infile, 'r') as f:
		data = f.read().split('\n')[1:]
		for line in data:
			if not line:
				continue
			barcode, sample, target, _, percent_cov = line.split(',')
			if barcode not in data_dict:
				data_dict[barcode] = (sample, PriorityQueue())
			data_dict[barcode][1].put((-float(percent_cov), target))

	ordered_keys = sorted(data_dict.keys())
	for barcode in ordered_keys:
		sys.stdout.write('{}: {}\n'.format(
			barcode,
			data_dict[barcode][0]
		))
		n_refs = n
		if data_dict[barcode][1].qsize() < n:
			n_refs = data_dict[barcode][1].qsize()
		for i in range(n_refs):
			percent_cov, target = data_dict[barcode][1].get(block=False)
			percent_cov = -percent_cov
			if 0. <= percent_cov < 20.:
				bar_col = '\u001b[41m'  # 41 = bg red
			elif 20. <= percent_cov < 80:
				bar_col = '\u001b[43m'  # 43 = bg yellow
			else:
				bar_col = '\u001b[42m'  # 42 = bg green
			if i == 0:
				intermediate_best_genomes[barcode] = target
			count = int(math.floor(BAR_LENGTH * percent_cov / 100))
			short_target = target.split('|')[-1].rstrip('_')
			second_gap = int(max(1, (TERMINAL_LENGTH - len(short_target) - BAR_LENGTH)))
			sys.stdout.write('\t{}:{}|{}{}{}{}| {:.2f}%\n'.format(
				short_target,
				' ' * second_gap,
				bar_col,
				' ' * count,
				end_col,
				' ' * (BAR_LENGTH - count),
				percent_cov
			))
		sys.stdout.write('\n')
	return intermediate_best_genomes


parser = argparse.ArgumentParser('asffast.py')
parser.add_argument('-i', '--input', type=str, default=None, required=True,
					help='Existing input root data directory for ONT sequencing run (must contain the desired sample name)')
parser.add_argument('-o', '--output', type=str, default=None, required=True,
					help='Path to desired result output directory')
parser.add_argument('-r', '--reference', type=str, default=None, required=True,
					help='Path to multi-FASTA of reference genomes to use as the alignment database')
parser.add_argument('-w', '--working_dir', type=str, default='/tmp',
					help='Path to working directory for the Python watcher script')
parser.add_argument('-s', '--singularity', type=str, default=None,
                    help='Path to Singularity container if other than default (pulls from cloud if this argument isn\'t used)')
parser.add_argument('-p', '--prefix', type=str, default=None,
                    help='Optional flowcell or sample ID, defaults to Flowcell1')
parser.add_argument('--slurm', action='store_true', default=False,
                    help='Flag for use of SLURM for HPC clusters.  Modify {} to change cluster options.'.format(
						SLURM_CONFIG
                    ))
parser.add_argument('--wait', type=int, default=1,
                    help='Number of seconds to watch for new data files before finalizing [1, inf)')


if __name__ == '__main__':
	args = parser.parse_args()
	nextflow_work_dir = get_create_real_dir(args.working_dir)
	nextflow_path = '/'.join(get_real_dir(sys.argv[0]).split('/')[:-1]) + '/asffast.nf'
	nextflow_config = nextflow_path.replace('asffast.nf', 'nextflow.config')
	globstar_input_path = None

	# Calculate thread number dynamically
	avail_threads = available_cpu_count()
	if args.slurm:  # BWA will be run with 1 thread but forked into max available CPUs here for efficiency on cluster
		threads = 1
	else:
		threads = str(int(max((1, math.ceil(float(avail_threads) / float(NEXTFLOW_MAX_FORKS))))))

	# Setup data structures
	samples = set()  # Either tuples of samplenames (no barcodes) or tuples of (samplename, barcode)
	files_present = 0
	watch_timer = 0
	cancel_flag = False  # Flag to cancel the run prematurely
	barcode_flag = False  # Flag indicating that barcodes are present in the run
	observed_best_genomes = None

	# Calculate input structure, find barcode ids
	sys.stdout.write('\n')
	while (watch_timer < args.wait) and not cancel_flag:
		this_file_counter = 0
		fastqs = glob.iglob(args.input + '/**/fastq_pass/**/*.fastq', recursive=True)
		if not fastqs:
			sys.stderr.write('\nNo folder \"fastq_pass\" or .fastq files detected in subdirectories of {}\n'.format(
				args.input
			))
			raise ValueError
		for f in fastqs:
			# Handles unclassified reads from barcoded runs
			samplename = '_'.join(f.split('/')[-1].split('.')[0].split('_')[:-1])
			barcode = BARCODE_REGEX.search(f)
			if barcode:
				barcode_flag = True
				samples.add((samplename, barcode.group(1)))
				if not globstar_input_path:
					globstar_input_path = get_real_dir(f).replace(barcode.group(1), 'barcode*') + '/*.fastq'
			else:
				parent_fastq_dir = f.split('/')[-2]
				if parent_fastq_dir == 'unclassified':
					continue
				samples.add((samplename,))
				if not globstar_input_path:
					globstar_input_path = get_real_dir(f) + '/*.fastq'
			this_file_counter += 1

		# Determine if Nextflow needs to be run, else start timer
		if this_file_counter > files_present:
			sys.stdout.write('{} new files found, Nextflow run triggered...\n\n'.format(
				this_file_counter - files_present
			))
			if barcode_flag:
				nextflow_arglist = [
					nextflow_path,
					'--reads',
					globstar_input_path,
					'--output',
					os.path.realpath(args.output),
					'--db',
					os.path.realpath(args.reference),
					'--threads',
					threads,
					'--barcodes',
					'-w',
					nextflow_work_dir,
					'-config',
					nextflow_config,
					'-resume'
				]
			else:
				nextflow_arglist = [
					nextflow_path,
					'--reads',
					globstar_input_path,
					'--output',
					os.path.realpath(args.output),
					'--threads',
					threads,
					'--db',
					os.path.realpath(args.reference),
					'-w',
					nextflow_work_dir,
					'-config',
					nextflow_config,
					'-resume'
				]
			if args.singularity:
				nextflow_arglist += ['-with-singularity', args.singularity]
			else:
				nextflow_arglist += ['-with-singularity', SINGULARITY_LIB_PATH]
			if args.slurm:
				config_file_idx = nextflow_arglist.index('-config') + 1
				nextflow_arglist[config_file_idx] = nextflow_path.replace('asffast.nf', SLURM_CONFIG)
			p = subprocess.Popen([str(x) for x in nextflow_arglist])
			exit_code = p.wait()
			sys.stdout.write('\n')
			sys.stdout.flush()

			observed_best_genomes = report_intermediate_coverage(args.output + '/CoverageAnalysis/coverage_results.csv')

			files_present = this_file_counter
			watch_timer = 0
		else:
			time.sleep(1)
			watch_timer += 1
			sys.stdout.write('Waiting for new files for {} seconds\r'.format(watch_timer))
			sys.stdout.flush()
	sys.stdout.write('\n\n')

	# Start final Nextflow run with genomic subset of reference
	sys.stdout.write('Beginning final Nextflow run...\n\n')

	# Write best observed genomes for each barcode to file
	obs_genome_filepath = nextflow_work_dir + '/observed_best_genomes.txt'
	with open(obs_genome_filepath, 'w') as obs:
		for barcode, genome in observed_best_genomes.items():
			obs.write('{}\t{}\n'.format(barcode, genome))

	# Find metadata files if they exist and optionally include in the final run
	fq_pass = list(glob.iglob(args.input + '/**/fastq_pass', recursive=True))
	if len(fq_pass) > 1:
		sys.stderr.write('Error: Multiple fastq_pass folders were found in data directory when only one'
		                 ' should be present.  Please re-run the script by pointing to the top-level '
		                 'folder output by MinKNOW for a single flowcell.')
		raise ValueError
	fq_pass_dir = fq_pass[0]
	metadata_dir = os.path.realpath(fq_pass_dir.split('fastq_pass')[0])
	sequencing_files = glob.glob(metadata_dir + '/sequencing_summary_*')
	throughput_files = glob.glob(metadata_dir + '/throughput_*')

	globstar_sequencing_path = None
	globstar_througput_path = None
	if len(sequencing_files) > 0:
		globstar_sequencing_path = metadata_dir + '/sequencing_summary_*'
	if len(throughput_files) > 0:
		globstar_througput_path = metadata_dir + '/throughput_*'

	if barcode_flag:
		nextflow_arglist = [
			nextflow_path,
			'--reads',
			globstar_input_path,
			'--output',
			os.path.realpath(args.output),
			'--db',
			os.path.realpath(args.reference),
			'--threads',
			threads,
			'--barcodes',
			'--final_info',
			obs_genome_filepath,
			'-w',
			nextflow_work_dir,
			'-config',
			nextflow_config,
			'-resume'
		]
	else:
		nextflow_arglist = [
			nextflow_path,
			'--reads',
			globstar_input_path,
			'--output',
			os.path.realpath(args.output),
			'--threads',
			threads,
			'--db',
			os.path.realpath(args.reference),
			'--final_info',
			obs_genome_filepath,
			'-w',
			nextflow_work_dir,
			'-config',
			nextflow_config,
			'-resume'
		]
	if args.singularity:
		nextflow_arglist += ['-with-singularity', args.singularity]
	else:
		nextflow_arglist += ['-with-singularity', SINGULARITY_LIB_PATH]
	if args.slurm:
		config_file_idx = nextflow_arglist.index('-config') + 1
		nextflow_arglist[config_file_idx] = nextflow_path.replace('asffast.nf', SLURM_CONFIG)
	if globstar_sequencing_path and globstar_througput_path:
		nextflow_arglist += ['--throughput', globstar_througput_path, '--sequencing', globstar_sequencing_path]
	if args.prefix:
		nextflow_arglist += ['--out_prefix', args.prefix]
	p = subprocess.Popen([str(x) for x in nextflow_arglist])
	exit_code = p.wait()
	sys.stdout.write('\n')
	sys.stdout.flush()
