# VarHound - TSO500 - Coverage diagnostics launcher

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(Cairo))
suppressMessages(library(mclust))

args <- (commandArgs(TRUE))

source(args[1])
source(args[2])

covdir <- args[3]
runtype <- args[4]

if (as.numeric(args[5]) > 1) {
	vh.covrun(covdir, runtype)
} else {
	vh.covgen(covdir, runtype)
}
