#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys

print("# VarHound - VCF / annotation file cleaner")

outdir = os.path.join(os.path.dirname(sys.argv[1]), "clean")
if not os.path.exists(outdir):
	os.mkdir(outdir)

infiles = []
for root, sub, files in os.walk(sys.argv[1]):
	for x in files:
		infiles += [os.path.join(root, x)]

for infile in infiles:
	
	base = os.path.basename(infile)
	out = os.path.join(outdir, base)
	w = open(out, "w")
	
	with open(infile, "r") as han:
		for line in han:
			if line.startswith("##"):
				pass
			elif line.startswith("#"):
				w.write(line.lstrip("#"))
			else:
				w.write(line)
	
	w.close()
	print("#", base + ': cleaned')
