#!/usr/bin/env python3

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from samscore.samscore import SamParser
from samscore.samscore import score_cigar


if __name__ == '__main__':
	sam_parser = SamParser(sys.argv[1])
	for qheader, _, theader, tstart, cigar in sam_parser:
		score, idx_scores, matches, mismatches, insertions, deletions = score_cigar(cigar, tstart - 1)
		print(qheader, theader, score, len(matches), len(mismatches), len(insertions), len(deletions))
