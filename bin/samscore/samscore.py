# Author: Steven Lakin

import re
import sys


class SamParser(object):
	"""
	Line-by-line parsing of Sequence Alignment Map (SAM) UTF-8 encoded files.  Stores alignment information for
	genome lengths from the provided reference database using the SAM headers, if specified.  Outputs
	aligned reads if specified and only if the aligned reads do not match a provided filter file (headers of reference
	to exclude).  Aligning reads, if output, are output as they are parsed.
	"""

	def __init__(self, sam_path, capture_ref_len=True, aln_scores=True):
		"""
		Initialize object and data structures for storing coverage and coverage over time, if specified.
		:param sam_path: STR, path to input SAM file
		:param capture_ref_len: BOOL, store reference genome lengths from SAM headers
		:param aln_scores: BOOL, return alignment scores if present
		"""
		self.aln_regex = None
		if aln_scores:
			self.aln_regex = re.compile(r'AS:i:([0-9]+)\t')
			self.dp_regex = re.compile(r'ms:i:([0-9]+)\t')
		self.sam_path = sam_path
		self.ref_len = {}
		self.output_handle = None
		self.handle = None

		# Values
		self.query_header = None
		self.reverse = None
		self.target_header = None
		self.target_start = None
		self.cigar = None
		self.aln_score = None

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
		:yield: (query, query_reverse_bool, target, target_start_idx, CIGAR, aln_score)
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
		self.query_header = entries[0]
		self.target_header = entries[2]
		self.target_start = int(entries[3])
		self.cigar = entries[5]
		if sam_flag & 16 != 0:  # if reverse
			self.reverse = True
		else:
			self.reverse = False
		if self.aln_regex:
			try:
				self.aln_score = int(self.aln_regex.search(self.line).group(1))
			except AttributeError:
				self.aln_score = int(self.dp_regex.search(self.line).group(1))
		else:
			self.aln_score = None
		self.line = self.handle.readline()
		return self.values()

	def _open(self):
		self.handle = open(self.sam_path, 'r')

	def _close(self):
		self.handle.close()

	def values(self):
		return self.query_header, self.reverse, self.target_header, self.target_start, self.cigar, self.aln_score


class ReadScoreCache(object):
	"""

	"""
	def __init__(self):
		self.cache = {}
		self.current_read = ''

	def _put(self, target_name, read_score_obj, idxs):
		if target_name not in self.cache:
			self.cache[target_name] = [(read_score_obj, idxs)]
		else:
			self.cache[target_name].append((read_score_obj, idxs))

	def _clear(self):
		self.cache = {}

	def _top(self):
		finals = {}
		for target, read_cache in self.cache.items():
			finals[target] = sorted(read_cache)[-1]
		return True, finals

	def smart_insert(self, read_name, target_name, read_score_obj, read_idxs):
		top_val = (None, None)
		if not read_name:
			sys.stderr.write('Error: read_name is undefined\n')
			raise ValueError
		if not self.current_read:
			self.current_read = read_name
		if read_name != self.current_read:
			top_val = self._top()
			self._clear()
		self._put(target_name, read_score_obj, read_idxs)
		return top_val

	def finalize(self):
		return self._top()


def parse_cigar(s, t_idx):
	"""
	Parse SAM CIGAR alignment string and return indices to which the read aligned.
	:param s: STR, CIGAR string
	:param t_idx: INT, zero-based index for target start position
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


def score_cigar(s, t_idx, match=2, mismatch=-4, indel_start=-2, indel_extend=-1):
	"""
	Parse SAM CIGAR alignment string and return reference indices for matches, mismatches, insertions, and deletions.
	All indices are zero-indexed.  The score is calculated with respect to each position in the reference genome.
	The overall read-wise score is the sum of the first tuple.  This score should roughly correlate with alignment
	scores from seed-and-extend aligners.
	:param s: STR, CIGAR string
	:param t_idx: INT, zero-based index for target start position
	:param match: INT, score/penalty for an alignment match (match of sequence to reference)
	:param mismatch: INT, score/penalty for an alignment mismatch (mismatch of sequence to reference)
	:param indel_start: INT, score/penalty for starting an insertion/deletion event
	:param indel_extend: INT, score/penalty for extension of an existing insertion/deletion event
	:return: tuple of tuples/lists of integers, (score, [idx_scores], (match_idxs,), (mismatch_idxs,), (insert_idxs,), (delete_idxs,))
	"""
	score = 0
	idx_scores = []
	match_idxs = tuple()
	mismatch_idxs = tuple()
	insert_idxs = tuple()
	delete_idxs = tuple()
	num = ''
	c_idx = 0
	start_idx = t_idx
	while c_idx < len(s):
		if s[c_idx].isdigit():
			num += s[c_idx]
		else:
			op = s[c_idx]
			if op == 'M' or op == '=':
				match_idxs += tuple(range(t_idx, t_idx + int(num)))
				idx_scores += [match for _ in range(int(num))]
				score += match * int(num)
				t_idx += int(num)
			elif op == 'D':
				delete_idxs += tuple(range(t_idx, t_idx + int(num)))
				idx_scores += [indel_start]
				idx_scores += [indel_extend for _ in range(int(num) - 1)]
				score += (indel_extend * (int(num) - 1)) + indel_start
				t_idx += int(num)
			elif op == 'N' or op == 'X':
				mismatch_idxs += tuple(range(t_idx, t_idx + int(num)))
				idx_scores += [mismatch for _ in range(int(num))]
				score += mismatch + int(num)
				t_idx += int(num)
			elif op == 'I':
				insert_cost = (indel_extend * (int(num) - 1)) + indel_start
				insert_idxs += (t_idx,)
				idx_scores[-1] += insert_cost
				score += insert_cost
			num = ''
		c_idx += 1
	assert ((t_idx - start_idx) == len(idx_scores))
	return score, idx_scores, match_idxs, mismatch_idxs, insert_idxs, delete_idxs


if __name__ == '__main__':
	sam_parser = SamParser(sys.argv[1])
	read_cache = ReadScoreCache()
	for header, rev, thead, tstart, cigar, aln_score in sam_parser:
		cigar_score = score_cigar(cigar, tstart - 1)
		aln_idxs = parse_cigar(cigar, tstart - 1)
		top_read, top_idx_dict = read_cache.smart_insert(header, thead, cigar_score, aln_idxs)
		if top_read:
			for target, top_idxs in top_idx_dict.items():
				sys.stdout.write('{}\t{}\n'.format(target, len(top_idxs)))

	_, final_idx_dict = read_cache.finalize()
	for target, final_idxs in final_idx_dict.items():
		sys.stdout.write('{}\t{}\n'.format(target, final_idxs))
