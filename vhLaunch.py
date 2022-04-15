#! /usr/bin/python3
# -*- coding: utf-8 -*-

# VarHound - TSO500 - Coverage diagnostics launcher

#  Copyright (C) 2022 Fernando Palluzzi
#  e-mail: <fernando.palluzzi@gmail.com>
#  Bioinformatics facility 
#  Gemelli Science and Technological Park (GSTeP)
#  Fondazione Policlinico Universitario Agostino Gemelli IRCCS,
#  Largo Agostino Gemelli 8, 00168 Roma, Italy

#  VarHound is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  VarHound is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys
import gzip
import shutil
from configparser import ConfigParser
from subprocess import call

if len(sys.argv) == 1:
	exit("\n          ## VarHound software ##\n" +\
	     "# Coverage Analysis (VarHound-CA) - " +\
	     "v0.1.0 #\n\n" +
	     "BASIC USAGE:\n  vhLaunch.py COVERAGE_DIRECTORY " +\
	     "[snv|cnv|rna|fusion]\n" +\
	     "Example:\n" +\
	     "  vhLaunch.py ~/TSO500 cnv\n\n" +\
	     "PER-BASE ANALYSIS:\n  vhLaunch.py COVERAGE_DIRECTORY " +\
	     "x<COVERAGE_THRESHOLD>\n" +\
	     "Example:\n" +\
	     "  vhLaunch.py ~/TSO500 x80\n")

VHCONF = os.path.expanduser('~') + "/VarHound/conf/paths.ini"
cp = ConfigParser()
CONF = cp.read(VHCONF)

VarHound = cp.get("path", "VarHound")
suffix = cp.get("path", "suffix")
pbcov = cp.get("path", "pbcov")
manifest = cp.get("path", "manifest")
manifest = os.path.join(VarHound, manifest)
reference = cp.get("path", "reference")
reference = os.path.join(VarHound, reference)

covpath = os.path.normpath(sys.argv[1])
vhrunx = ' ' + os.path.join(VarHound, "vhLaunch.R")
vhcore = ' ' + os.path.join(VarHound, "varhound.R")
vhcovr = ' ' + os.path.join(VarHound, "vhCoverage.R")
vhpbca = ' ' + os.path.join(VarHound, "vhPerbase.py")
covdir = ' ' + covpath

print("\n# VarHound-CA suite started: coverage data preparation ...")

perbase = False
if len(sys.argv) > 2:
	if sys.argv[2].startswith("x"):
		p = ' ' + sys.argv[2].split("x")[1]
		perbase = True
		suffix = pbcov
	else:
		run = ' ' + sys.argv[2]
else:
	run = ' snv'

F = []
for root, sub, files in os.walk(covpath):
	for x in files:
		if x.endswith(suffix + '.gz'):
			f = os.path.join(root, x)
			if not perbase:
				F += [f]
			with gzip.open(f, 'rb') as infile:
				with open(f.rstrip('.gz'), 'wb') as outfile:
					shutil.copyfileobj(infile, outfile)
suffix = ' ' + suffix

print ("# Done.")

if perbase:
	manifest = ' ' + manifest
	reference = ' ' + reference
	call(vhpbca + covdir + suffix + p + manifest + reference,
	     shell = True)
else:
	F = ' ' + str(len(F))
	call('Rscript' + vhrunx + vhcore + vhcovr + covdir + run + F + suffix,
		 shell = True)
