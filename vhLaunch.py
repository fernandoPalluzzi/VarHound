#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
from subprocess import call

# VarHound absolute path
VarHound = "~/varhound"

# Coverage (unzipped) file extension
suffix = "thresholds.bed"


covpath = os.path.normpath(sys.argv[1])
vhrunx = ' ' + os.path.join(VarHound, "vhLaunch.R")
vhcore = ' ' + os.path.join(VarHound, "varhound.R")
vhcovr = ' ' + os.path.join(VarHound, "vhCoverage.R")
covdir = ' ' + covpath

call('gzip -rdkf' + covdir, shell = True)

f = []
for root, sub, files in os.walk(covpath):
	for x in files:
		if x.endswith(suffix):
			f += [os.path.join(root, x)]

f = ' ' + str(len(f))

if len(sys.argv) > 2:
	run = ' ' + sys.argv[2]
else:
	run = ' snv'

call('Rscript' + vhrunx + vhcore + vhcovr + covdir + run + f, shell = True)
