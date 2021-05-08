# Author: Steven Lakin

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


def score_cigar(s, t_idx, match=1, mismatch=-1, indel_start=-2, indel_extend=-1):
	"""
	Parse SAM CIGAR alignment string and return reference indices for matches, mismatches, insertions, and deletions.
	All indices are zero-indexed.  The score is calculated with respect to each position in the reference genome,
	using +1 for match, -1 for mismatch, -2 for indel start, and -1 for indel extension (affine gap penalty).
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
	print((t_idx - start_idx), len(idx_scores))
	assert((t_idx - start_idx) == len(idx_scores))
	return score, idx_scores, match_idxs, mismatch_idxs, insert_idxs, delete_idxs
