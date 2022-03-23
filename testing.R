#  VarHound - TSO500 - Coverage diagnostics

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


# R

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(Cairo))
suppressMessages(library(mclust))

#source("~/VarHound/vhLaunch.R")
source("~/VarHound/varhound.R")
source("~/VarHound/vhCoverage.R")


### Unit test 1: Input header.

x <- read.delim("~/varhound_dev/TSO500_COVERAGE/test/16AP_DNA.thresholds.bed", stringsAsFactors = FALSE)
#h <- colnames(x)
#m <- regexpr("[xX]*\\d+[xX]*", h)
#h <- c(h[m == -1], paste0("x", str_replace_all(h[m > 0], "[xX]", "")))
h <- cov.header(x)
h


### Unit Test 2: Coverage data acquisition.

x <- cov.preprocess("~/varhound_dev/TSO500_COVERAGE/test/16AP_DNA.thresholds.bed")
head(x)
dim(x)


### Unit Test 3: Covdata file format.

# Longtype coverage data
# [chrom, start, end, region, ID, time, coverage, length, p]
x.long <- cov.long("~/varhound_dev/TSO500_COVERAGE/test/16AP_DNA.thresholds.bed")
head(x.long)

# Covdata format
# [chrom, start, end, region, ID, time, {SAMPLES_COVERAGE}]
covdata <- cov.data("~/varhound_dev/TSO500_COVERAGE/test/210715_A01423_0008_AH35CWDRXY/coverage")
head(covdata)


### Unit Test 4: Sequencing run-level coverage aggregation.

# Listing available depths
depths <- h$header[seq(1, h$n.depths) + h$offset]
depths

run <- cov.aggregate(covdata, t = depths)

# Extended covdata format (median)
# [chrom, start, end, region, ID, time, {SAMPLES_COVERAGE},
#  median, min, max, Coverage, Depth]
head(run$covdata)

# Run-level coverage stacked barplot
plot(run$barplot)

# Extended covdata format (Q1)
Q1 <- cov.q1(covdata, t = depths)
head(Q1$covdata)
plot(Q1$bar)


### Unit Test 5: Sample-level covdata aggregation.

# Preparing run-level covdata aggregates
run <- vh.panel(covdir = "~/varhound_dev/TSO500_COVERAGE/test/210715_A01423_0008_AH35CWDRXY/coverage",
                runtype = "snv", format = "bed", depths = depths,
				plot.name = "wholePanel_exonCoverage.png",
				y.label = "Percent of covered exons",
				outdir = "./")

R <- vh.reshape(run$covdata, level = "sample")
M <- vh.reorder(R, depth = "x100", level = "sample")
rDataset <- as.factor(M$rDataset[M$Depth == "x100"])
Depth <- as.factor(M$Depth[M$Depth == "x100"])
head(M)


### VarHoundCA: step-by-step coverage analysis in R


# R

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(Cairo))
suppressMessages(library(mclust))

source("~/VarHound/varhound.R")
source("~/VarHound/vhCoverage.R")


# Input settings

runtype <- "rna"

depths <- c("5x", "10x", "50x", "100x", "250x", "500x")

covdir <- "~/varhound_dev/TSO500_COVERAGE/Validation/rna_coverage_test"

outdir <- paste0(covdir, "/CoverageDiagnostics")
if (!file.exists(outdir)) dir.create(outdir)


# Whole-panel diagnostics

run <- vh.panel(covdir = covdir, runtype = runtype,
                format = "bed", depths = depths,
				plot.name = "wholePanel_exonCoverage.png",
				y.label = "Percent of covered exons",
				outdir = outdir)

head(run$covdata)


# Sample-level diagnostics at 5x
x5 <- vh.sample(run$covdata, runtype = runtype, depth = depths[1],
				color = "lavender", line.color = "orchid",
				outdir = outdir)


# Sample-level diagnostics at 10x
x10 <- vh.sample(run$covdata, runtype = runtype, depth = depths[2],
				 color = "azure", line.color = "dodgerblue",
				 outdir = outdir)


# Sample-level diagnostics at 50x
x50 <- vh.sample(run$covdata, runtype = runtype, depth = depths[3],
				 color = "lightblue", line.color = "darkblue",
				 outdir = outdir)


# Sample-level diagnostics at 100x
x100 <- vh.sample(run$covdata, runtype = runtype, depth = depths[4],
				  color = "green3", line.color = "darkgreen",
				  outdir = outdir)


# Sample-level diagnostics at 250x
x250 <- vh.sample(run$covdata, runtype = runtype, depth = depths[5],
				  color = "gold", line.color = "brown",
				  outdir = outdir)


# Sample-level diagnostics at 500x
x500 <- vh.sample(run$covdata, runtype = runtype, depth = depths[6],
				  color = "darkorange", line.color = "darkred",
				  outdir = outdir)


# Sample-level reprot preparation

samples <- vh.covtable(x500$M, level = "sample", depths)


# Diagnostic plots for possible low-coverage genes

blacklist <- vh.blacklist(samples, run$covdata, runtype = runtype,
                          outdir = outdir)


# Diagnostic tables for possible low-coverage samples and genes

bl.genes <- vh.covtable(blacklist$B$ref$M, level = "gene", depths = depths)

vh.covreport(samples, bl.genes, runtype = runtype, outdir = outdir)


# BED files containing blacklisted exons

ebl5 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                     d = depths[1],
	                     outdir = outdir)

ebl10 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                      d = depths[2],
	                      outdir = outdir)

ebl50 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                      d = depths[3],
	                      outdir = outdir)

ebl100 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                       d = depths[4],
	                       outdir = outdir)

ebl250 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                       d = depths[5],
	                       outdir = outdir)

ebl500 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                       d = depths[6],
	                       outdir = outdir)


### VarHoundCA: single-sample gene coverage profiling


# R

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(Cairo))
suppressMessages(library(mclust))

source("~/VarHound/varhound.R")
source("~/VarHound/vhCoverage.R")


# Input settings

runtype <- "rna"

depths <- c("5x", "10x", "50x", "100x", "250x", "500x")

covdir <- "~/varhound_dev/TSO500_COVERAGE/Validation/rna_genecov_test"

outdir <- paste0(covdir, "/CoverageDiagnostics")
if (!file.exists(outdir)) dir.create(outdir)


# Single-sample gene coverage profiling

genecov <- vh.genecov(covdir = covdir, p = 75, d = "500x",
                      depths = depths,
                      runtype = runtype,
                      outdir = outdir)

# BED files containing blacklisted exons

ebl5 <- vh.exonblack(genecov$Q1, runtype = runtype,
                     d = depths[1],
                     outdir = outdir)

ebl10 <- vh.exonblack(genecov$Q1, runtype = runtype,
                      d = depths[2],
                      outdir = outdir)

ebl50 <- vh.exonblack(genecov$Q1, runtype = runtype,
                      d = depths[3],
                      outdir = outdir)

ebl100 <- vh.exonblack(genecov$Q1, runtype = runtype,
                       d = depths[4],
                       outdir = outdir)

ebl250 <- vh.exonblack(genecov$Q1, runtype = runtype,
                       d = depths[5],
                       outdir = outdir)

ebl500 <- vh.exonblack(genecov$Q1, runtype = runtype,
                       d = depths[6],
                       outdir = outdir)







