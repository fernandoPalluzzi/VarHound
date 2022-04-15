#! /usr/bin/python3
# -*- coding: utf-8 -*-

# VarHound - TSO500 - Per-base coverage diagnostics

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
from subprocess import call

def setHeader(f, header):
	tmp = f.split('.')[0] + '_setheader.tmp'
	header = '\t'.join(header.split()) + '\n'
	w = open(tmp, 'w')
	w.write(header)
	with open(f, 'r') as lines:
		for line in lines:
			w.write(line)
	w.close()
	os.rename(tmp, f)
	return f

def expandField(x, out, i = 4, sep = ",", header = False):
	i = i - 1
	w = open(out, 'w')
	with open(x, 'r') as lines:
		for line in lines:
			line = line.strip().split()
			if header:
				hdr = line
				unp = line[i]
				header = False
			else:
				J = line[i].split(sep)
				for j in J:
					w.write('\t'.join(line[:i] + [j] + line[i+1:]) + '\n')
	w.close()

def bedTrack(x, out, name = "", description = "", rgb = "204,0,0"):
	w = open(out, 'w')
	w.write('track name="' + name + '" description="' + description +\
	        '" useScore=1 itemRgb="ON" visibility=full\n')
	w.write('#chrom\tstart\tend\tname\tscore\tstrand\t' +\
	        'thickStart\tthickEnd\titemRgb\n')
	
	with open(x, 'r') as lines:
		for line in lines:
			line = line.strip().split()
			score = round(float(line[4]))
			line = line[:4] + [str(score), '.'] + line[1:3] + [rgb]
			w.write('\t'.join(line) + '\n')
	w.close()

F = []
for root, sub, files in os.walk(sys.argv[1]):
	for x in files:
		if x.endswith(sys.argv[2]):
			F += [os.path.join(root, x)]

for x in F:
	
	print('\n# Processing file ' + os.path.basename(x) + ' ...')
	
	tmp = x + '.tmp'
	blist = x.split('.')[0] + '_coverageDrops.bed'
	track = x.split('.')[0] + '_coverageDrops_track.bed'
	track_name = os.path.basename(x).split('.')[0]
	depth = sys.argv[3] + 'x'
	track_desc = track_name + ' blacklisted regions at ' + depth
	drops = x.split('.')[0] + '_geneDrops.bed'
	
	call('sort -k1,1 -k2,2n ' + x + ' > ' + tmp, shell = True)
	os.rename(tmp, x)
	
	call('awk -v OFS=\"\t\" \'{ if ($4 < ' + sys.argv[3] + ') ' +\
	     'print }\' ' + x + ' > ' + blist,
	     shell = True)
	
	call('bedtools merge -i ' + blist + ' > ' + tmp, shell = True)
	os.rename(tmp, blist)
	
	call('bedtools map -a ' + blist + ' -b ' + sys.argv[4] +\
	     ' -c 4 -o distinct > ' + tmp,
	     shell = True)
	expandField(tmp, out = blist)
	os.remove(tmp)
	
	call('awk -v OFS=\"\t\" \'{ if ($4 != \".\") print }\' ' +\
	     blist + ' > ' + tmp,
	     shell = True)
	os.rename(tmp, blist)
	
	call('bedtools map -a ' + blist + ' -b ' + x +\
	     ' -c 4 -o median > ' + tmp,
	     shell = True)
	expandField(tmp, out = blist)
	os.remove(tmp)
	
	bedTrack(blist, out = track,
	         name = track_name,
	         description = track_desc)
	
	call('grep Exon ' + blist + ' > ' + tmp, shell = True)
	
	setHeader(blist, 'chrom dropStart dropEnd region medianCoverage')
	
	call('bedtools map -a ' + sys.argv[5] + ' -b ' + tmp +\
	     ' -c 4 -o count > ' + drops,
	     shell = True)
	os.remove(tmp)
	
	setHeader(drops, 'chrom geneStart geneEnd symbol ' +\
	          'exonCount strand nExonDrops_' +\
	          depth)
	
	print('# Done.')

print()
