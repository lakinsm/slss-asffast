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
		print(score_cigar(cigar, tstart - 1))
