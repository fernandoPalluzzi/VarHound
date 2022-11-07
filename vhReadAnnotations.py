#! /usr/bin/python3
# -*- coding: utf-8 -*-

# VarHound - vhPlot library - VCF annotations descriptives

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
import argparse
import numpy as np
import pandas as pd
import vhPlot as vhp
from shutil import rmtree

print("\n# VarHound - vhPlot library - VCF annotations descriptives")


parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("indir", help = "input directory")

parser.add_argument("-e", "--extension",
                    help = "input file extension",
                    dest = "ext",
                    default = "txt",
                    type = str)

parser.add_argument("-c", "--vclass",
                    help = "variant class field name",
                    dest = "vclass",
                    default = "VARIANT_CLASS",
                    type = str)

parser.add_argument("-i", "--impact",
                    help = "variant impact field name",
                    dest = "impact",
                    default = "IMPACT",
                    type = str)

parser.add_argument("-o", "--order",
                    help = "order of the variant impact levels",
                    dest = "order",
                    default = ["HIGH", "MODERATE", "LOW", "MODIFIER"],
                    type = str,
                    nargs = "+")

parser.add_argument("-id", "--identifier",
                    help = "annotation element identifier (e.g., gene symbol)",
                    dest = "id",
                    default = ["SYMBOL", "labels"],
                    type = str,
                    nargs = 2)

parser.add_argument("-t", "--target",
                    help = "target files for variant frequency aggregates",
                    dest = "target",
                    default = ["geneCounts", "ClinVar"],
                    type = str,
                    nargs = 2)

parser.add_argument("-k", "--topk",
                    help = "top-k annotations (e.g., genes) shown in barplots",
                    dest = "top",
                    default = 50,
                    type = int)

parser.add_argument("-w", "--weight",
                    help = "score weights for high, moderate, and low impact variants",
                    dest = "weight",
                    default = [1.1, 1, 0.1, 0.01],
                    type = float,
                    nargs = 3)

parser.add_argument("-x", "--exclude",
                    help = "missing data value(s)",
                    dest = "exclude",
                    default = ["None", "NaN", "NA", "-"],
                    type = str,
                    nargs = "+")

parser.add_argument("-v", "--verbose",
                    help = "show execution information to screen",
                    dest = "verb",
                    default = "store_false")

args = parser.parse_args()


infiles = [os.path.join(sys.argv[1], x) for x in os.listdir(sys.argv[1]) \
           if os.path.isfile(os.path.join(sys.argv[1], x)) \
           and x.endswith(args.ext)]

for infile in infiles:
	
	if args.verb:
		print("\n[INFO]    Computation started for:",
		      os.path.basename(infile), "\n")
		print("[INFO]    Data acquisition ...")
	
	x = vhp.importTable(infile)
	
	if args.verb:
		print("[INFO]    Done.\n")
	
	tree = vhp.setdirs(infile)
	
	vhp.varCharts(x, varfreq = tree[0], clinsig = tree[1], clinvar = tree[2],
				  gcounts = tree[3], stackedbar = tree[4], freqbar = tree[5],
				  cake1 = tree[6], cake2 = tree[7], gvariants = tree[8],
				  rawvar = tree[9], verbose = args.verb, hue = args.impact,
				  order = args.order)
	
	vhp.classCharts(x, infile, attr_vclass = args.vclass, verbose = args.verb)

if args.verb:
	print("\n[INFO]    Generating overall descriptives ...")

tables = os.path.join(sys.argv[1], "overall/trash")
if not os.path.exists(tables):
	os.mkdir(tables)

if args.verb:
	print("[INFO]    Done.")
	print("\n[INFO]    Collecting variant frequencies ...")

f = vhp.fetchCountpaths(sys.argv[1], target = args.target[0])

f = vhp.combineCounts(f, ID = args.id[0], where = tables, verbose = True)

if args.verb:
	print("\n[INFO]    Gathering frequencies by annotation ...")

g = vhp.fetchCountpaths(sys.argv[1], target = args.target[1])

g = vhp.combineCounts(g, ID = args.id[1], suffix = '_impact.txt',
                      where = tables, verbose = True)

if args.verb:
	print("\n[INFO]    Generating variant frequency plots ...")


for k in f.keys():
	
	outdir = os.path.join(sys.argv[1], "overall", k)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	
	out = os.path.join(outdir, "Genewise_varFreqs_" + k + ".pdf")
	
	m = vhp.selectMutants(f[k], symbols = args.id[0], top = args.top,
	                      high_impact = args.order[0],
			    		  moderate_impact = args.order[1],
			    		  low_impact = args.order[2],
				    	  modifier = args.order[3],
				    	  w1 = args.weight[0],
				    	  w2 = args.weight[1],
				    	  w3 = args.weight[2],
					      exclude = args.exclude,
					      computeFreqs = False)
	
	vhp.countbars(m, where = out, vclass = args.id[0],
	              check_impact = True, topgenes = args.top, angle = 0,
	              xlab = 'log2 variant frequency', ylab = 'Gene',
	              stacked = True, vertical = False,
	              exclude = args.exclude)
	
	out = os.path.join(outdir, "Genewise_rawFreqs_" + k + ".pdf")
	
	raw = vhp.rawVarcounts(f[k], symbols = args.id[0], ids = args.id[0],
	                       top = args.top, high_impact = args.order[0],
	                       moderate_impact = args.order[1],
	                       low_impact = args.order[2],
	                       modifier = args.order[3],
	                       exclude = args.exclude,
	                       computeFreqs = False)
	
	vhp.countbars(raw, where = out, vclass = args.id[0], angle = 0,
	              xlab = 'log2 variant frequency', ylab = 'Gene',
	              stacked = True, vertical = False, gscores = True,
	              exclude = args.exclude, check_impact = True,
	              showgrid = True)
	
	out = os.path.join(outdir, "ClinVar_" + k + ".pdf")
	
	vhp.cake(g[k], where = out, exclude = args.exclude,
	         label = args.id[1], mincount = 1, maxcount = None,
	         dpi = 450, width = 20, height = 16,
	         title = 'Overall variant impact', titlesize = 22,
	         fontsize = 16, dec = 2, pctdist = 0.8, donut = True,
	         radius = 0.6, detach = 0.6, basicdonut = False,
	         silent = True, warning_title = "",
	         pal = ["#0080FF", "#66B2FF", "#99CCFF", "#CCCCFF",
	                "#00E100", "#FAE100", "#FFB266", "#FFCCFF",
	                "#FF9999", "#FF6666", "#FF0000", "#C0C0C0",
	                "#E0E0E0"])

if args.verb:
	print("[INFO]    Done.")
	print("\n[INFO]    Writing count tables ...")

countfiles = [os.path.join(tables, x) for x in os.listdir(tables) \
              if x.find('_counts.txt') > -1]

for countfile in countfiles:
	
	x = vhp.importTable(countfile)
	x = vhp.checkImpact(x, high = args.order[0],
	                    moderate = args.order[1],
	                    low = args.order[2],
	                    modifier = args.order[3])
	
	x = vhp.setScore(x, symbols = args.id[0],
	                 high_impact = args.order[0],
	                 moderate_impact = args.order[1],
	                 low_impact = args.order[2],
	                 modifier = args.order[3],
	                 w1 = args.weight[0],
	                 w2 = args.weight[1],
	                 w3 = args.weight[2],
	                 w4 = args.weight[3],
	                 exclude = args.exclude)
	
	x = x[[args.id[0], args.order[0], args.order[1], args.order[2],
	       args.order[3], 'rawscore', 'hiscore']]
	
	base = os.path.basename(countfile)
	dest = os.path.join(sys.argv[1], "overall", base.split('_')[0])
	
	x.to_csv(os.path.join(dest, base), index = False, sep = "\t")

if args.verb:
	print("[INFO]    All tasks completed.")
