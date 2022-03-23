#! /usr/bin/python3
# -*- coding: utf-8 -*-

# VarHound - TSO500 - Coverage diagnostics launcher

import os
import sys
import gzip
import shutil
from subprocess import call

# VarHound absolute path
VarHound = "~/VarHound"

# Coverage (unzipped) file extension
suffix = "thresholds.bed"

if len(sys.argv) == 1:
	exit("VarHound software - v0.0.2\n" +
	      "USAGE:\n  vhLaunch.py COVERAGE_DIRECTORY [snv|cnv|rna|fusion]")

covpath = os.path.normpath(sys.argv[1])
vhrunx = ' ' + os.path.join(VarHound, "vhLaunch.R")
vhcore = ' ' + os.path.join(VarHound, "varhound.R")
vhcovr = ' ' + os.path.join(VarHound, "vhCoverage.R")
covdir = ' ' + covpath

F = []
for root, sub, files in os.walk(covpath):
	for x in files:
		if x.endswith(suffix + '.gz'):
			f = os.path.join(root, x)
			F += [f]
			with gzip.open(f, 'rb') as infile:
				with open(f.rstrip('.gz'), 'wb') as outfile:
					shutil.copyfileobj(infile, outfile)
F = ' ' + str(len(F))
suffix = ' ' + suffix

if len(sys.argv) > 2:
	run = ' ' + sys.argv[2]
else:
	run = ' snv'

call('Rscript' + vhrunx + vhcore + vhcovr + covdir + run + F + suffix,
     shell = True)
