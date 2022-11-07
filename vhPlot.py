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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


def freqs(x, vclass = 'VARIANT_CLASS', impact = 'IMPACT',
          log2 = False, melt = False, export = None):
	
	w = pd.crosstab(index = x[vclass], columns = x[impact])
	w = w[set(x[impact])]
	w[vclass] = list(w.index)
	
	if export != None:
		w.to_csv(export, index = False, sep = "\t")
	
	if log2:
		w = w.apply(lambda x: np.log2(x + 1) if np.issubdtype(x.dtype,
		                                        np.number) else x)
	
	if melt:
		w = pd.melt(w, id_vars = vclass)
	
	return w


def fetchCountpaths(path, target = 'geneCounts', nodeFilter = 'byClass',
                    depth = 2):
	f = {}
	for root, sub, files in os.walk(path):
		for x in files:
			newnode = os.path.join(root, x)
			if target in newnode and nodeFilter in newnode:
				key = newnode.split('/')[-1*(depth + 1)]
				if key in f.keys():
					f[key] += [newnode]
				else:
					f[key] = [newnode]
	return f


def dfCombine(x1, x2, ID):
	
	L1 = list(x1)
	L2 = list(x2)
	
	for a in set(L1 + L2):
		if a not in L1 and a in L2:
			x1[a] = 0
		elif a in L1 and a not in L2:
			x2[a] = 0
		if a != ID:
			x1[a] = x1[a].add(x2[a], fill_value = 0)
	
	return x1


def combineCounts(paths, ID = 'SYMBOL', suffix = '_counts.txt',
                  where = None, verbose = False):
	f = {}
	for k in paths.keys():
		x = importTable(paths[k][0])
		if len(paths[k]) > 1:
			for i in range(1, len(paths[k])):
				x = dfCombine(x, importTable(paths[k][i]), ID = ID)
		f[k] = x
		if where != None:
			x.to_csv(os.path.join(where, k + suffix),
			         index = False, sep = "\t")
		if verbose:
			print('[INFO]    Counts collected for class:', k)
	return(f)


def setScore(x, symbols = 'SYMBOL', high_impact = 'HIGH',
             moderate_impact = 'MODERATE', low_impact = 'LOW',
             modifier = 'MODIFIER', w1 = 1.1, w2 = 1, w3 = 0.1,
             w4 = 0.01, exclude = ['NA', '-'], ascending = False):
	
	x = x[~x[symbols].isin(exclude)]
	
	x['rawscore'] = np.log2(x[high_impact] + x[moderate_impact] +\
	                        x[low_impact] + x[modifier] + 1)
	
	x['hiscore'] = np.log2(w1*x[high_impact] + w2*x[moderate_impact] +\
	                       w3*x[low_impact] + 1)
	
	x = x.sort_values(by = 'rawscore', ascending = ascending)
	
	return x


def selectMutants(x, symbols = 'SYMBOL', high_impact = 'HIGH',
                  moderate_impact = 'MODERATE', low_impact = 'LOW',
                  modifier = 'MODIFIER', gscore = True, log2 = True,
                  ascending = False, nonzero = False, top = 30,
                  exclude = ['NA', '-'], w1 = 1.3, w2 = 1, w3 = 0.01,
                  export = None, computeFreqs = True):
	
	if computeFreqs:
		x = freqs(x, vclass = symbols)
	
	g = x[~x[symbols].isin(exclude)]
	
	if export != None:
		g.to_csv(export, index = False, sep = "\t")
	
	if gscore:
		
		if high_impact in list(g):
			h = g[high_impact]
		else:
			h = 0
		if moderate_impact in list(g):
			m = g[moderate_impact]
		else:
			m = 0
		if low_impact in list(g):
			l = g[low_impact]
		else:
			l = 0
		
		if type(h) == int and type(m) == int:
			print('[WARNING] No high/moderate-impact variation found.' +\
			      ' Using low-impact variation.')
			if low_impact in list(g):
				l = g[low_impact]
			else:
				print('[ERROR]   Failed. Execution aborted.')
				exit()
		
		if log2:
			g['gscore'] = np.log2(w1*h + w2*m + w3*l + 1)
		else:
			g['gscore'] = w1*h + w2*m + w3*l
		
		g = g[g['gscore'] > 0]
		g = g.sort_values(by = 'gscore', ascending = ascending)
		
	else:
		g = g.sort_values(by = list(g)[:-1], ascending = ascending)
	
	if nonzero:
		try:
			g = g[g[high_impact] > 0]
		except:
			print('[WARNING] Skipped: high-impact variation not found.')
		try:
			g = g[g[moderate_impact] > 0]
		except:
			print('[WARNING] Skipped: moderate-impact variation not found.')
	
	g = g.head(n = top)
	
	return g


def rawVarcounts(x, symbols = 'SYMBOL', high_impact = 'HIGH',
                 moderate_impact = 'MODERATE', low_impact = 'LOW',
                 modifier = 'MODIFIER', ascending = False, top = 30,
                 exclude = ['NA', '-'], ids = 'RefSeq',
                 computeFreqs = True):
	
	if ids != None:
		x = x[~x[ids].isin(exclude)]
	
	if computeFreqs:
		x = freqs(x, vclass = symbols)
	
	if high_impact in list(x):
		h = x[high_impact]
	else:
		h = 0
	
	if moderate_impact in list(x):
		m = x[moderate_impact]
	else:
		m = 0
	
	if low_impact in list(x):
		l = x[low_impact]
	else:
		l = 0
	
	if modifier in list(x):
		w = x[modifier]
	else:
		w = 0
	
	g = x[~x[symbols].isin(exclude)]
	g['gscore'] = np.log2(h + m + l + w + 1)
	g = g.sort_values(by = 'gscore', ascending = ascending)
	g = g.head(n = top)
	
	return g


def clinsigEncode(x, nullvalues):
	
	csd = {'benign': "B", 'likely_benign': "B",
	       'pathogenic': "P", 'likely_pathogenic': "P",
	       'uncertain_significance': "V",
	       'Uncertain_significance': "V"}
	
	try:
		x = "".join(set([csd[w] for w in x if w not in nullvalues]))
		if x == "BV":
			return "Benign/VUS"
		elif x == "B":
			return "Benign/Likely_benign"
			#return "Benign"
		elif x == "P":
			return "Pathogenic/Likely_pathogenic"
			#return "Pathogenic"
		elif x == "V":
			return "VUS"
		elif x == "PV":
			return "Pathogenic/VUS"
		else:
			return "Benign/Pathogenic"
	except:
		x = "Other"


def clinsigRecode(v, sep = ",", nullvalues = ['not_provided', '-']):
	
	csd = {'benign': "Benign", 'likely_benign': "Likely_benign",
	       'pathogenic': "Pathogenic",
	       'likely_pathogenic': "Likely_pathogenic",
	       'uncertain_significance': "VUS",
	       'Uncertain_significance': "VUS",
	       'Conflicting_interpretations_of_pathogenicity': "Conflicting"}
	
	recoded = []
	for x in v:
		x = x.split(sep)
		if len(x) > 1:
			x = clinsigEncode(x, nullvalues)
		else:
			if x[0] in nullvalues:
				x = "NA"
			else:
				try:
					x = csd[x[0]]
				except:
					if x[0] in csd.values():
						x = x[0]
					else:
						x = "Other"
		recoded += [x]
	return recoded


def setColor(x, colors, labels = ["Benign", "Benign/Likely_benign",
                                  "Likely_benign", "Benign/VUS", "VUS",
                                  "Conflicting", "Benign/Pathogenic",
                                  "Pathogenic/VUS", "Likely_Pathogenic",
                                  "Pathogenic/Likely_Pathogenic",
                                  "Pathogenic", "Other", "NA"]):
	
	keys = list(x['labels'])
	data = []
	
	for i in range(len(labels)):
		if labels[i] in keys:
			data += [x['counts'][keys.index(labels[i])]]
		else:
			data += [0]
	
	return pd.DataFrame({"labels": labels, "data": data, "color": colors})


def naWarning(x, label, warning_title = ""):
	
	if "NA" in list(x[label]):
		
		nas = int(x[x[label] == "NA"].iloc[:, 1])
		nna = x[x[label] != "NA"].iloc[:, 1].sum()
		pct = round(100*nna/(nna + nas), 3)
		
		if warning_title != "":
			warning_title = warning_title + ': '
		
		print('[WARNING] ' + warning_title + str(nna) + '/' +\
		      str(nna + nas) + ' non-NA variants used (' +\
		      str(pct) + '%)')


def checkImpact(x, high = 'HIGH', moderate = 'MODERATE', low = 'LOW',
                modifier = 'MODIFIER'):
	
	pd.options.mode.chained_assignment = None
	L = len(x.axes[0])
	
	if high not in list(x):
		x.loc[:, high] = [0]*L
	
	if moderate not in list(x):
		x.loc[:, moderate] = [0]*L
	
	if low not in list(x):
		x.loc[:, low] = [0]*L
	
	if modifier not in list(x):
		x.loc[:, modifier] = [0]*L
	
	return x


def countbars(x, where = 'mycountbars.pdf', dpi = 450, stacked = False,
              pal = ['#CC0000', '#FF8000', '#FFFF00', '#99CCFF'],
              showgrid = True, vclass = 'VARIANT_CLASS', label = 'labels',
              impact = 'IMPACT', check_impact = False, hue = None,
              order = None, value = 'value', vertical = True,
              title = 'Overall variant frequency', xlab = 'Variant class',
              ylab = 'Frequency', angle = 45, titlesize = 20,
              xlabsize = 16, ylabsize = 16, width = 18, height = 10,
              topgenes = 0, log2 = False, score = 'gscore', gscores = False,
              exclude = ['NA', '-'], silent = False):
	
	if exclude != None:
		if label not in list(x):
			label = vclass
		x[label] = ["NA" if w == "-" else w for w in x[label]]
		if "NA" in list(x[label]) and not silent:
			naWarning(x, label = label)	
		x = x[~x[label].isin(exclude)]
	
	if log2:
		x = x.apply(lambda w: np.log2(w + 1) if np.issubdtype(w.dtype,
		                                        np.number) else w)
	
	if check_impact:
		x = checkImpact(x)
	
	if topgenes > 0:
		if not vertical:
			x = x.sort_values(by = score, ascending = True)
		#x = x.iloc[:, :-1]
		x = x[['HIGH', 'MODERATE', 'LOW', vclass]]
		if len(x.axes[0]) != topgenes:
			topgenes = len(x.axes[0])
		title = title + ' (top ' + str(topgenes) + ' genes)'
		
	elif gscores:
		#pal = ['#00CC00']
		#x = x.iloc[:, :-1]
		if not vertical:
			#x = x.sort_values(by = x.columns[1], ascending = True)
			x = x.sort_values(by = x.columns[-1], ascending = True)
		#x = x.iloc[:, :-1]
		x = x[['HIGH', 'MODERATE', 'LOW', 'MODIFIER', vclass]]
		if len(x.axes[0]) == topgenes:
			title = title + ' (top ' + str(topgenes) + ' genes)'	
	
	if showgrid:
		sns.set(style = 'whitegrid',
		        rc = {"figure.figsize": (width, height)})
	else:
		sns.set(style = 'white',
		        rc = {"figure.figsize": (width, height)})
	
	if stacked:
		if order != None:
			x = x[order + [vclass]]
		if vertical:
			x.set_index(vclass).plot(kind = 'bar', stacked = True,
			                         color = pal)
		else:
			x.set_index(vclass).plot(kind = 'barh', stacked = True,
			                         color = pal)
	else:
		if len(list(x)) != 3:
			x = pd.melt(x, id_vars = vclass)
		if vertical:
			sns.barplot(data = x, x = vclass, y = value, hue = hue,
						hue_order = order, palette = pal)
		else:
			sns.barplot(data = x, y = vclass, x = value, hue = hue,
						hue_order = order, palette = pal)
	
	plt.title(title, fontsize = titlesize)
	plt.xlabel(xlab, fontsize = xlabsize)
	plt.ylabel(ylab, fontsize = ylabsize)
	plt.xticks(rotation = angle)
	#plt.show()
	plt.savefig(where, dpi = dpi)
	plt.clf()
	plt.close()


def cake(counts, where = 'mycake.pdf', exclude = ["NA"], label = 'labels',
         mincount = 1, maxcount = None, dpi = 450, width = 20,
         height = 16, title = 'Overall variant impact', titlesize = 20,
         fontsize = 10, dec = None, shadow = False, explode = None,
         pctdist = 0.6, donut = False, radius = 0.6, detach = 0,
         basicdonut = False, silent = False, warning_title = "",
         pal = ["#0080FF", "#66B2FF", "#99CCFF", "#CCCCFF",
                "#00E100", "#FAE100", "#FFB266", "#FFCCFF",
                "#FF9999", "#FF6666", "#FF0000", "#C0C0C0",
                "#E0E0E0"]):
	
	pie = setColor(counts, colors = pal)
		
	if exclude != None:
		
		if "NA" in list(pie['labels']) and not silent:
			naWarning(pie, label = 'labels', warning_title = warning_title)
		
		pie = pie[~pie['labels'].isin(exclude)]
	
	if exclude == "clean":
		exclude = ["NA", "Benign", "Benign/Likely_benign", "Likely_benign"]
	
	if mincount != None:
		pie = pie[pie['data'] >= mincount]
	
	if maxcount != None:
		pie = pie[pie['data'] <= maxcount]
	
	#sns.set(rc = {"figure.figsize": (width, height)})
	
	if dec != None:
		dec = '%.' + str(dec) + 'f%%'
	
	plt.pie(pie['data'], labels = pie['labels'],
	        colors = pie['color'],
	        autopct = dec,
	        pctdistance = pctdist,
	        shadow = shadow,
	        explode = explode,
	        textprops = {'fontsize': fontsize},
	        wedgeprops = {'linewidth': detach, 'edgecolor': 'white' })
	
	if donut:
		centre_circle = plt.Circle((0, 0), radius, fc = 'white')
		fig = plt.gcf()
		fig.gca().add_artist(centre_circle)
	
	if title != None:
		plt.title(title, fontsize = titlesize)
	
	#plt.show()
	plt.savefig(where, dpi = dpi)
	plt.clf()
	plt.close()


def setdirs(infile, by = "sample", vclass = ""):
	
	overall = os.path.join(os.path.dirname(infile), "overall")
	if not os.path.exists(overall):
		os.mkdir(overall)
	
	base = os.path.join(os.path.dirname(infile),
						os.path.basename(infile).split('.')[0])
	if not os.path.exists(base):
		os.mkdir(base)
	
	if by == "sample":
		
		byx = os.path.join(base, "bySample")
		if not os.path.exists(byx):
			os.mkdir(byx)
	
	elif by == "class":
		
		by = os.path.join(base, "byClass")
		if not os.path.exists(by):
			os.mkdir(by)		
		
		byx = os.path.join(base, "byClass", vclass)
		if not os.path.exists(byx):
			os.mkdir(byx)
	
	else:
		raise ValueError('Analysis can be either by "sample" or by "class"')
	
	counts = os.path.join(byx, "counts")
	if not os.path.exists(counts):
		os.mkdir(counts)
	
	varfreq = os.path.join(counts,
			  os.path.basename(infile).split('.')[0] + '_varFreqs.txt')
	clinsig = os.path.join(counts,
			  os.path.basename(infile).split('.')[0] + '_ClinSig.txt')
	clinvar = os.path.join(counts,
			  os.path.basename(infile).split('.')[0] + '_ClinVar.txt')
	gcounts = os.path.join(counts,
			  os.path.basename(infile).split('.')[0] + '_geneCounts.txt')
	
	stackedbar = os.path.join(byx,
	             os.path.basename(infile).split('.')[0] +\
	             '_overall_variants_stacked.pdf')
	
	freqbar = os.path.join(byx,
	          os.path.basename(infile).split('.')[0] +\
	          '_overall_variants_freq.pdf')
	
	cake1 = os.path.join(byx,
	        os.path.basename(infile).split('.')[0] + '_ClinSig_cake.pdf')
	
	cake2 = os.path.join(byx,
	        os.path.basename(infile).split('.')[0] + '_ClinVar_cake.pdf')
	
	gvariants = os.path.join(byx,
	            os.path.basename(infile).split('.')[0] +\
	            '_genewise_varFreqs.pdf')
	
	rawvar = os.path.join(byx,
	         os.path.basename(infile).split('.')[0] + '_genewise_rawFreqs.pdf')
	
	tree = (varfreq, clinsig, clinvar, gcounts, stackedbar, freqbar,
	        cake1, cake2, gvariants, rawvar)
	
	return tree


def importTable(infile, sep = "\t", lineterm = "\n", lowmem = False,
                      verbose = False):
	
	if verbose:
		print('[INFO]    Reading input file ...')
	
	x = pd.read_csv(infile, sep = sep, lineterminator = lineterm,
	                low_memory = lowmem)
	
	if verbose:
		print('[INFO]    Done.')
		print('[INFO]    Annotation file:', os.path.basename(infile))
	
	return x


def varCharts(x, bysample = True, varfreq = None, clinsig = None,
              clinvar = None, gcounts = None, stackedbar = None,
              freqbar = None, cake1 = None, cake2 = None, gvariants = None,
              rawvar = None, log2 = True, attr_vclass = None,
              vclass = 'VARIANT_CLASS', impact = 'IMPACT', hue = None,
              order = None, exclude = ["NA", "-"], label = 'labels',
              bar_angle = 0, bar_ylab = "log2 frequency",
              attr_clinsig = 'CLIN_SIG', attr_clinvar = 'ClinVar_CLNSIG',
              title_clinsig = None, title_clinvar = None, donut = True,
              pie_detach = 0.6, pie_decimals = 1, pie_pctdist = 0.8,
              pie_fontsize = 16, pie_titlesize = 22, top = 50,
              geneids = 'SYMBOL', mut_title = "Gene variation countbars",
              mut_xlab = "log2 variant frequency", mut_ylab = "Gene",
              raw_xlab = "log2 variant frequency", raw_ylab = "Gene",
              showgrid = True, verbose = True):
	
	if varfreq == None:
		varfreq = "Variant_freqs.txt"
	if clinsig == None:
		clinsig = "ClinSig_freqs.txt"
	if clinvar == None:
		clinvar = "ClinVar_freqs.txt"
	if gcounts == None:
		gcounts = "Gene_counts.txt"
	if stackedbar == None:
		stackedbar = "Overall_variants_stacked.pdf"
	if freqbar == None:
		freqbar = "Overall_variants_freq.pdf"
	if cake1 == None:
		cake1 = "ClinSig_cake.pdf"
	if cake2 == None:
		cake2 = "ClinVar_cake.pdf"
	if gvariants == None:
		gvariants = "Genewise_varFreqs.pdf"
	if rawvar == None:
		rawvar = "Genewise_rawFreqs.pdf"
	
	if verbose:
		print('[INFO]    Generating countbars ...')
	
	if bysample:
		
		f = freqs(x, vclass = vclass, impact = impact, log2 = log2,
		          export = varfreq)
		
		countbars(f, vclass = vclass, impact = impact, where = stackedbar,
		          angle = bar_angle, ylab = bar_ylab, showgrid = showgrid,
		          stacked = True, exclude = exclude, hue = hue,
		          order = order, check_impact = True)
		
		countbars(f, vclass = vclass, impact = impact, where = freqbar,
		          angle = bar_angle, ylab = bar_ylab, showgrid = showgrid,
		          exclude = exclude, hue = hue, order = order,
		          check_impact = True)
	
	if title_clinvar == None:
		title_clinvar = "Overall variant impact\n(" + attr_clinvar + ")"

	if title_clinsig == None:
		title_clinsig = "Overall variant impact\n(" + attr_clinsig + ")"
	
	c1 = clinsigRecode(x[attr_clinvar])
	c1 = pd.Series(c1).value_counts()
	c1 = c1.rename_axis('labels').reset_index(name = 'counts')
	c1.to_csv(clinvar, index = False, sep = "\t")
	
	if verbose:
		print('[INFO]    Done.')
	
	if bysample:
		
		if verbose:
			print('[INFO]    Generating impact pies ...')

		c2 = clinsigRecode(x[attr_clinsig])
		c2 = pd.Series(c2).value_counts()
		c2 = c2.rename_axis('labels').reset_index(name = 'counts')
		c1.to_csv(clinsig, index = False, sep = "\t")
		
		cake(c1, where = cake1, donut = True, detach = pie_detach,
			 dec = pie_decimals, pctdist = pie_pctdist, exclude = exclude,
			 fontsize = pie_fontsize, titlesize = pie_titlesize,
			 title = title_clinsig, warning_title = attr_clinsig)
		
		cake(c2, where = cake2, donut = True, detach = pie_detach,
			 dec = pie_decimals, pctdist = pie_pctdist, exclude = exclude,
			 fontsize = pie_fontsize, titlesize = pie_titlesize,
			 title = title_clinvar, warning_title = attr_clinvar)
	
	if verbose:
		print('[INFO]    Generating gene variation charts ...')
	
	if not bysample:
		
		palette = ["#0080FF", "#66B2FF", "#99CCFF", "#CCCCFF", "#00E100",
		           "#FAE100", "#FFB266", "#FFCCFF", "#FF9999", "#FF6666",
		           "#FF0000", "#C0C0C0", "#E0E0E0"]
		
		c3 = setColor(c1, colors = palette)
		c3 = c3[c3['data'] > 0]
		c3 = c3.sort_values(by = 'data', ascending = False)
		color = c3['color'][~c3[label].isin(exclude)]
		
		countbars(c3, where = freqbar, vclass = label, value = 'data',
				  angle = bar_angle, title = mut_title,
				  xlab = 'log2 frequency', ylab = 'ClinVar',
				  exclude = exclude, vertical = False,
				  showgrid = showgrid, pal = color)
	
	g = selectMutants(x, top = top, export = gcounts)
	countbars(g, where = gvariants, vclass = geneids, check_impact = True,
			  topgenes = top, angle = bar_angle, title = mut_title,
			  xlab = mut_xlab, ylab = mut_ylab, stacked = True,
			  vertical = False, exclude = exclude, showgrid = showgrid)
	
	raw = rawVarcounts(x, top = 50)
	countbars(raw, where = rawvar, vclass = geneids, check_impact = True,
			  angle = 0, xlab = raw_xlab, ylab = raw_ylab, stacked = True,
			  vertical = False, gscores = True, exclude = exclude,
			  showgrid = showgrid)
	
	if verbose:
		print('[INFO]    Done.')


def classCharts(x, infile, attr_vclass = 'VARIANT_CLASS',
                exclude = ["NA", "-"], verbose = True):
	
	for vclass in set(x[attr_vclass]):
		
		v = x[x[attr_vclass] == vclass]
		
		if verbose:
			print('\n[INFO]    Analysing class:', vclass)
		
		tree = setdirs(infile, by = "class", vclass = vclass)
		
		varCharts(v, varfreq = tree[0], clinsig = tree[1],
		          clinvar = tree[2], gcounts = tree[3],
		          stackedbar = tree[4], freqbar = tree[5],
		          cake1 = tree[6], cake2 = tree[7],
		          gvariants = tree[8], rawvar = tree[9],
		          exclude = exclude, bysample = False)

