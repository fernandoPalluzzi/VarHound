# VarHound - TSO500 - Coverage diagnostics

#library(dplyr)
#library(ggplot2)
#library(reshape2)
#library(Cairo)
#library(mclust)

vh.bar <- function(bardata, plot.name, w = 20, h = 10, u = 'in', r = 450, cairo = TRUE) {
	if (cairo) {
		CairoPNG(plot.name, width = w, height = h, units = u, res = r)
		print(bardata)
		dev.off()
	} else {
		png(plot.name, width = w, height = h, units = u, res = r)
		print(bardata)
		dev.off()

	}
}

vh.reshape <- function(covdata, level,
                       depths = c("50x", "100x", "250x", "500x"),
                       init = 9) {
	
	if (level == "sample") {
		level <- "PatientID"
		a <- ncol(covdata)
		R <- reshape2::melt(covdata[, c(7:(a-5), a)], id = "Depth")
		names(R) <- c("Depth", level, "Coverage")
	
	} else if (level == "gene") {
		covdata <- cov.aggregate(covdata, init = init,
		                         t1 = depths[1], t2 = depths[2],
		                         t3 = depths[3], t4 = depths[4],
		                         barplot = FALSE)$covdata
		R <- covdata[, colnames(covdata) %in% c("Depth", "symbol", "median")]
		names(R) <- c("Gene", "Coverage", "Depth")
		
	} else {
		stop("invalid level choice.")
	}
	return(R)
}

vh.reorder <- function(R, depth, level) {
	if (level == "sample") {
		reordered <- R %>% filter(Depth == depth) %>%
		                   mutate(rDataset = reorder(PatientID,
		                          Coverage,
		                          function(x) 1/(mean(x)+1))) %>%
		                          select(PatientID, rDataset)
		M <- R %>% inner_join(reordered, by = "PatientID")
		return(M)
	} else if (level == "gene") {
		reordered <- R %>% filter(Depth == depth) %>%
		                   mutate(rDataset = reorder(Gene,
		                          Coverage,
		                          function(x) 1/(mean(x)+1))) %>%
		                          select(Gene, rDataset)
		M <- R %>% inner_join(reordered, by = "Gene")
		return(M)
	} else {
		stop("invalid level choice")
	}
}

vh.sample <- function(covdata, depth, level = "sample", runtype = "snv",
                      w = 20, h = 10, u = 'in', r = 450, out.size = -1,
                      color = "gold", line.color = "darkred",
                      x.lab = "Patient ID",
                      y.lab = "Median Percent Coverage (MPC)",
                      outdir = "") {
	
	R <- vh.reshape(covdata, level = level)
	M <- vh.reorder(R, depth = depth, level = level)
	rDataset <- as.factor(M$rDataset[M$Depth == depth])
	Depth <- as.factor(M$Depth[M$Depth == depth])
	
	if (!is.null(runtype)) {
		if (outdir != "") outdir <- paste0(outdir, "/")
		plot.name <- paste0(outdir, "TSO500_", toupper(runtype),
		                    "_coverage_byPatient_",
							depth, ".png")
		png(plot.name, width = w, height = h, units = u, res = r)
		print(ggplot(M[M$Depth == depth,],
			   aes(x = rDataset, y = Coverage, fill = Depth)) +
			   geom_boxplot(outlier.size = out.size) +
			   scale_fill_manual(name = "Depth", values = color) +
			   theme_bw() +
			   theme(panel.border = element_blank(),
			   panel.grid.major = element_line(),
			   panel.grid.minor = element_blank(),
			   axis.line = element_line(colour = "black"),
			   axis.text = element_text(size = 14),
			   axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
			   axis.title = element_text(size = 26, face = "bold"),
			   legend.key.size = unit(1, "cm"),
			   legend.text = element_text(size = 22),
			   legend.title = element_text(size = 24)) +
			   geom_hline(yintercept = 75, linetype = "dashed", color = line.color) +
			   scale_y_continuous(limits = c(0, 100)) +
			   labs(x = x.lab, y = y.lab))
		dev.off()
	}
	
	return(list(M = M, rDataset = rDataset, Depth = Depth))
}

vh.panel <- function(covdir, cleanup = "Exon", runtype = "snv", include = NULL,
                     format = "bed", plot.name = NULL,
                     x.label = "Percent of covered regions",
                     y.label = "Depth of coverage",
                     w = 20, h = 10, u = 'in', r = 450, cairo = TRUE,
                     outdir = NULL) {
	
	W <- cov.data(directory = covdir, include = include, cleanup = cleanup,
	              runtype = runtype, format = format)
	
	W <- cov.aggregate(W, xlab = x.label, ylab = y.label)
	
	if (!is.null(plot.name) & !is.null(outdir)) {
		plot.name = paste0(outdir, "/TSO500_", toupper(runtype), "_", plot.name)
		vh.bar(W$barplot, plot.name, w, h, u, r, cairo = cairo)
	}
	return(W)
}

set.status <- function(Q, ctrl = c(22, 23, 24, 25),
                          depths = c("50x", "100x", "250x", "500x")) {
	status <- vector()
	for (j in 1:nrow(Q)) {
		q <- Q[j, ctrl]
		if (sum(q) == 0) {
			flag <- "VALID"
		} else {
			w <- which(q == 2)
			if (length(w) == 0) {
				m <- "WARNING"
				w <- min(which(q == 1))
			} else {
				m <- "INVALID"
				w <- min(w)
			}
			flag <- paste0(m, "_", depths[w])
		}
		status <- c(status, flag)
	}
	return(status)
}

vh.covtable <- function(M, level, depths = c("50x", "100x", "250x", "500x")) {
	
	if (level == "sample") {
		level <- "PatientID"
	} else if (level == "gene") {
		level <- "Gene"
	} else {
		stop("invalid level choice.")
	}
	
	Qmin <- aggregate(M[M$Depth == depths[1], 3],
	                  list(M[M$Depth == depths[1], 2]), min)
	names(Qmin) <- c(level, "Qmin50")
	Qmin$Qmin100 <- aggregate(M[M$Depth == depths[2], 3],
	                          list(M[M$Depth == depths[2], 2]),
	                          min)[, 2]
	Qmin$Qmin250 <- aggregate(M[M$Depth == depths[3], 3],
	                          list(M[M$Depth == depths[3], 2]),
	                          min)[, 2]
	Qmin$Qmin500 <- aggregate(M[M$Depth == depths[4], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          min)[, 2]
	
	Q1 <- aggregate(M[M$Depth == depths[1], 3],
	                list(M[M$Depth == depths[1], 2]),
	                function(x) quantile(x)[2])
	names(Q1) <- c(level, "Q1x50")
	Q1$Q1x100 <- aggregate(M[M$Depth == depths[2], 3],
	                       list(M[M$Depth == depths[2], 2]),
	                       function(x) quantile(x)[2])[, 2]
	Q1$Q1x250 <- aggregate(M[M$Depth == depths[3], 3],
	                       list(M[M$Depth == depths[3], 2]),
	                       function(x) quantile(x)[2])[, 2]
	Q1$Q1x500 <- aggregate(M[M$Depth == depths[4], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[2])[, 2]
	
	Q2 <- aggregate(M[M$Depth == depths[1], 3],
	                list(M[M$Depth == depths[1], 2]),
	                function(x) quantile(x)[3])
	names(Q2) <- c(level, "Q2x50")
	Q2$Q2x100 <- aggregate(M[M$Depth == depths[2], 3],
	                       list(M[M$Depth == depths[2], 2]),
	                       function(x) quantile(x)[3])[, 2]
	Q2$Q2x250 <- aggregate(M[M$Depth == depths[3], 3],
	                       list(M[M$Depth == depths[3], 2]),
	                       function(x) quantile(x)[3])[, 2]
	Q2$Q2x500 <- aggregate(M[M$Depth == depths[4], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[3])[, 2]
	
	Q3 <- aggregate(M[M$Depth == depths[1], 3],
	                list(M[M$Depth == depths[1], 2]),
	                function(x) quantile(x)[4])
	names(Q3) <- c(level, "Q3x50")
	Q3$Q3x100 <- aggregate(M[M$Depth == depths[2], 3],
	                       list(M[M$Depth == depths[2], 2]),
	                       function(x) quantile(x)[4])[, 2]
	Q3$Q3x250 <- aggregate(M[M$Depth == depths[3], 3],
	                       list(M[M$Depth == depths[3], 2]),
	                       function(x) quantile(x)[4])[, 2]
	Q3$Q3x500 <- aggregate(M[M$Depth == depths[4], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[4])[, 2]
	
	Qmax <- aggregate(M[M$Depth == depths[1], 3],
	                  list(M[M$Depth == depths[1], 2]), max)
	names(Qmax) <- c(level, "Qmax50")
	Qmax$Qmax100 <- aggregate(M[M$Depth == depths[2], 3],
	                          list(M[M$Depth == depths[2], 2]),
	                          max)[, 2]
	Qmax$Qmax250 <- aggregate(M[M$Depth == depths[3], 3],
	                          list(M[M$Depth == depths[3], 2]),
	                          max)[, 2]
	Qmax$Qmax500 <- aggregate(M[M$Depth == depths[4], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          max)[, 2]
	
	Q <- data.frame(Qmin[, 1:2], Q1x50 = Q1$Q1x50, Q2x50 = Q2$Q2x50,
	                Q3x50 = Q3$Q3x50, Qmax50 = Qmax$Qmax50,
	                Qmin100 = Qmin$Qmin100, Q1x100 = Q1$Q1x100,
	                Q2x100 = Q2$Q2x100, Q3x100 = Q3$Q3x100,
	                Qmax100 = Qmax$Qmax100, Qmin250 = Qmin$Qmin250,
	                Q1x250 = Q1$Q1x250, Q2x250 = Q2$Q2x250,
	                Q3x250 = Q3$Q3x250, Qmax250 = Qmax$Qmax250,
	                Qmin500 = Qmin$Qmin500, Q1x500 = Q1$Q1x500,
	                Q2x500 = Q2$Q2x500, Q3x500 = Q3$Q3x500,
	                Qmax500 = Qmax$Qmax500)
	
	Q$FLAGx50 <- ifelse(Q$Q1x50 >= 75, 0, 1)
	Q$FLAGx50[Q$Q3x50 == 0] <- 2
	Q$FLAGx100 <- ifelse(Q$Q1x100 >= 75, 0, 1)
	Q$FLAGx100[Q$Q3x100 == 0] <- 2
	Q$FLAGx250 <- ifelse(Q$Q1x250 >= 75, 0, 1)
	Q$FLAGx250[Q$Q3x250 == 0] <- 2
	Q$FLAGx500 <- ifelse(Q$Q1x500 >= 75, 0, 1)
	Q$FLAGx500[Q$Q3x500 == 0] <- 2
	
	if (level == "PatientID") {
		Q <- Q[order(Q$FLAGx500, Q$FLAGx250, Q$FLAGx100, Q$FLAGx50,
		             Q$PatientID),]
	} else if (level == "Gene") {
		Q <- Q[order(Q$FLAGx500, Q$FLAGx250, Q$FLAGx100, Q$FLAGx50,
		             Q$Gene),]
	}
	Q$status <- set.status(Q)
	Q <- Q[, c(1, 22:26, 2:21)]
	
	return(Q)
}

vh.genecov <- function(covdir, p = 75, d = "500x",
                       depths = c("50x", "100x", "250x", "500x"),
                       runtype = "snv", format = "bed",
                       w = 20, h = 10, u = 'in', r = 450,
                       color = c("lightblue", "green3", "gold", "darkorange"),
                       line.color = c("darkblue", "darkgreen", "brown", "darkred"),
                       x.lab = "Gene",
                       y.lab = "Exon Median Percent Coverage (MPC)",
                       outdir = "", init = 9) {
	
	x <- cov.data(directory = covdir, include = NULL,
	              cleanup = "Exon",
	              runtype = "snv",
	              format = format)
	
	covdata <- cov.aggregate(x, t1 = depths[1], t2 = depths[2],
	                            t3 = depths[3], t4 = depths[4],
	                            init = 7, barplot = FALSE,
	                            percent = TRUE,
	                            xlab = x.lab, ylab = y.lab,
	                            dummy = TRUE)$covdata
	
	info <- data.frame(t(sapply(sapply(covdata$region,
	                   function(j) strsplit(j, "_")), c)))
	names(info) <- c("symbol", "exon", "GeneID")
	
	covdata <- cbind(covdata[, 1:4], info, covdata[, 5:ncol(covdata)])
	R <- R <- covdata[, colnames(covdata) %in% c("Depth", "symbol", "median")]
	names(R) <- c("Gene", "Coverage", "Depth")
	
	M1 <- vh.reorder(R, depth = depths[1], level = "gene")[, c(3, 1, 2, 4)]
	rDataset1 <- as.factor(M1$rDataset[M1$Depth == depths[1]])
	Depth1 <- as.factor(M1$Depth[M1$Depth == depths[1]])
	
	M2 <- vh.reorder(R, depth = depths[2], level = "gene")[, c(3, 1, 2, 4)]
	rDataset2 <- as.factor(M2$rDataset[M2$Depth == depths[2]])
	Depth2 <- as.factor(M2$Depth[M2$Depth == depths[2]])
	
	M3 <- vh.reorder(R, depth = depths[3], level = "gene")[, c(3, 1, 2, 4)]
	rDataset3 <- as.factor(M3$rDataset[M3$Depth == depths[3]])
	Depth3 <- as.factor(M3$Depth[M3$Depth == depths[3]])
	
	M <- vh.reorder(R, depth = depths[4], level = "gene")[, c(3, 1, 2, 4)]
	rDataset <- as.factor(M$rDataset[M$Depth == depths[4]])
	Depth <- as.factor(M$Depth[M$Depth == depths[4]])
	
	B <- list(depth1 = list(M = M1, rDataset = rDataset1, Depth = Depth1),
	          depth2 = list(M = M2, rDataset = rDataset2, Depth = Depth2),
	          depth3 = list(M = M3, rDataset = rDataset3, Depth = Depth3),
	          ref = list(M = M, rDataset = rDataset, Depth = Depth))
	
	if (!is.null(runtype)) {
		if (outdir != "") outdir <- paste0(outdir, "/")
		for (j in 1:4) {
			png(paste0(outdir, "TSO500_", toupper(runtype),
			           "_coverage_byGene_", depths[j], ".png"),
			           width = w, height = h, units = u, res = r)
			
			B[[j]]$M$rDataset <- as.factor(B[[j]]$M$rDataset)
			
			print(ggplot(B[[j]]$M[B[[j]]$M$Depth == depths[j],],
			      aes(x = rDataset, y = Coverage, fill = Depth)) +
			      geom_boxplot(outlier.size = -1) +
			      scale_fill_manual(name = "Depth", values = color[j]) +
			      theme_bw() +
			      theme(panel.border = element_blank(),
			        panel.grid.major = element_line(),
			        panel.grid.minor = element_blank(),
			        axis.line = element_line(colour = "black"),
			        axis.text = element_text(size = 14),
			        axis.text.x = element_text(angle = 50,
			                                   vjust = 1, hjust = 1),
			        axis.title = element_text(size = 26, face = "bold"),
			        legend.key.size = unit(1, "cm"),
			        legend.text = element_text(size = 22),
			        legend.title = element_text(size = 24)) +
			      geom_hline(yintercept = 75, linetype = "dashed",
			                 color = line.color[j]) +
			      scale_y_continuous(limits = c(0, 100)) +
			      labs(x = x.lab, y = y.lab))
			dev.off()
		}
	}
	names(covdata)[11] <- "Q1"
	
	return(list(Q1 = covdata, B = B))
}

vh.blacklist <- function(samples, covdata, p = 75, d = "500x",
                         depths = c("50x", "100x", "250x", "500x"),
                         runtype = "snv", w = 20, h = 10, u = 'in', r = 450,
                         color = c("lightblue", "green3", "gold", "darkorange"),
                         line.color = c("darkblue", "darkgreen", "brown", "darkred"),
                         x.lab = "Gene",
                         y.lab = "Exon Median Percent Coverage (MPC)",
                         outdir = "", init = 9) {
	
	blacklist <- as.character(samples$PatientID[samples$status != "VALID"])
	if (length(samples) == 0) {
		stop("all samples are valid: no blacklisted regions found")
	}
	covdata <- covdata[, colnames(covdata) %in% c("chrom", "start", "end",
	                                              "region", "ID", "time",
	                                              blacklist)]
	info <- data.frame(t(sapply(sapply(covdata$region,
	                   function(j) strsplit(j, "_")), c)))
	names(info) <- c("symbol", "exon", "GeneID")
	
	covdata <- cbind(covdata[, 1:4], info, covdata[, 5:ncol(covdata)])
	R <- vh.reshape(covdata, level = "gene", depths = depths, init = init)
	
	M1 <- vh.reorder(R, depth = depths[1], level = "gene")[, c(3, 1, 2, 4)]
	rDataset1 <- as.factor(M1$rDataset[M1$Depth == depths[1]])
	Depth1 <- as.factor(M1$Depth[M1$Depth == depths[1]])
	
	M2 <- vh.reorder(R, depth = depths[2], level = "gene")[, c(3, 1, 2, 4)]
	rDataset2 <- as.factor(M2$rDataset[M2$Depth == depths[2]])
	Depth2 <- as.factor(M2$Depth[M2$Depth == depths[2]])
	
	M3 <- vh.reorder(R, depth = depths[3], level = "gene")[, c(3, 1, 2, 4)]
	rDataset3 <- as.factor(M3$rDataset[M3$Depth == depths[3]])
	Depth3 <- as.factor(M3$Depth[M3$Depth == depths[3]])
	
	M <- vh.reorder(R, depth = depths[4], level = "gene")[, c(3, 1, 2, 4)]
	rDataset <- as.factor(M$rDataset[M$Depth == depths[4]])
	Depth <- as.factor(M$Depth[M$Depth == depths[4]])
	
	B <- list(depth1 = list(M = M1, rDataset = rDataset1, Depth = Depth1),
	          depth2 = list(M = M2, rDataset = rDataset2, Depth = Depth2),
	          depth3 = list(M = M3, rDataset = rDataset3, Depth = Depth3),
	          ref = list(M = M, rDataset = rDataset, Depth = Depth))
	
	if (!is.null(runtype)) {
		if (outdir != "") outdir <- paste0(outdir, "/")
		for (j in 1:4) {
			png(paste0(outdir, "TSO500_", toupper(runtype),
			           "_coverage_byGene_", depths[j], ".png"),
			           width = w, height = h, units = u, res = r)
			
			B[[j]]$M$rDataset <- as.factor(B[[j]]$M$rDataset)
			
			print(ggplot(B[[j]]$M[B[[j]]$M$Depth == depths[j],],
			      aes(x = rDataset, y = Coverage, fill = Depth)) +
			      geom_boxplot(outlier.size = -1) +
			      scale_fill_manual(name = "Depth", values = color[j]) +
			      theme_bw() +
			      theme(panel.border = element_blank(),
			        panel.grid.major = element_line(),
			        panel.grid.minor = element_blank(),
			        axis.line = element_line(colour = "black"),
			        axis.text = element_text(size = 14),
			        axis.text.x = element_text(angle = 50,
			                                   vjust = 1, hjust = 1),
			        axis.title = element_text(size = 26, face = "bold"),
			        legend.key.size = unit(1, "cm"),
			        legend.text = element_text(size = 22),
			        legend.title = element_text(size = 24)) +
			      geom_hline(yintercept = 75, linetype = "dashed",
			                 color = line.color[j]) +
			      scale_y_continuous(limits = c(0, 100)) +
			      labs(x = x.lab, y = y.lab))
			dev.off()
		}
	}
	
	Q1 <- cov.q1(covdata, init = init, t1 = depths[1], t2 = depths[2],
		         t3 = depths[3], t4 = depths[4],
		         barplot = FALSE)$covdata
	
	return(list(Q1 = Q1, B = B))
}

vh.exonblack <- function(Q1, runtype, d = "500x", p = 75, outdir = "", report = TRUE) {
	exon.blacklist <- Q1[Q1$Depth == d & Q1$Q1 < p,
                         colnames(Q1) %in% c("chrom", "start", "end",
                                             "region", "symbol", "exon",
                                             "GeneID", "ID", "time",
                                             "Q1", "Depth")]
    exon.blacklist <- exon.blacklist[order(exon.blacklist$Q1),]
    
    if (outdir != "") outdir <- paste0(outdir, "/")
    if (report) {
		report <- paste0(outdir, "TSO500_", toupper(runtype),
		                 "_exonBlacklist_", d, ".bed")
		write.table(exon.blacklist, report,
		            quote = FALSE, row.names = FALSE, sep = "\t")
	}
	return(exon.blacklist)
}

vh.covreport <- function(samples, genes, runtype, outdir = "") {
	if (outdir != "") outdir <- paste0(outdir, "/")
	
	report.s <- paste0(outdir, "TSO500_", toupper(runtype),"_sampleCoverage.txt")
	write.table(samples, report.s,
	            quote = FALSE, row.names = FALSE, sep = "\t")
	
	report.g <- paste0(outdir, "TSO500_", toupper(runtype),"_geneCoverage.txt")
	write.table(genes, report.g,
	            quote = FALSE, row.names = FALSE, sep = "\t")
}

vh.covrun <- function(covdir, runtype = "snv",
                      depths = c("50x", "100x", "250x", "500x"),
                      suffix = "bed") {
	
	message("\n# Creating coverage diagnostics paths ...")
	if (runtype == "auto") {
		runtype <- strsplit(covdir, "/")[[1]]
		runtype <- tolower(runtype[length(runtype)])
	}
	if (!(runtype %in% c("snv", "cnv", "rna"))) stop("invalid run type.")
	
	outdir <- paste0(covdir, "/CoverageDiagnostics")
	if (!file.exists(outdir)) {
		dir.create(outdir)
	}
	message("# Done.\n")
	
	message("# Whole-run diagnostics ...")
	x <- vh.panel(covdir, runtype = runtype,
				  plot.name = "wholePanel_exonCoverage.png",
				  y.label = "Percent of covered exons",
				  outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[1], " ..."))
	x50 <- vh.sample(x$covdata, runtype = runtype, depth = depths[1],
					 color = "lightblue", line.color = "darkblue",
					 outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[2], " ..."))
	x100 <- vh.sample(x$covdata, runtype = runtype, depth = depths[2],
					  color = "green3", line.color = "darkgreen",
					  outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[3], " ..."))
	x250 <- vh.sample(x$covdata, runtype = runtype, depth = depths[3],
					  color = "gold", line.color = "brown",
					  outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[4], " ..."))
	x500 <- vh.sample(x$covdata, runtype = runtype, depth = depths[4],
					  color = "darkorange", line.color = "darkred",
					  outdir = outdir)
	message("# Done.\n")
	
	message("# Sample-level report preparation ...")
	Q <- vh.covtable(x500$M, level = "sample")
	message("# Done.\n")
	
	message("# Searching for blacklisted regions ...")
	blacklist <- vh.blacklist(Q, x$covdata, runtype = runtype, outdir = outdir)
	B <- vh.covtable(blacklist$B$ref$M, level = "gene")
	message("# Done.\n")
	
	message("# Gene-level reports preparation ...")
	
	vh.covreport(Q, B, runtype = runtype, outdir = outdir)
	
	ebl50 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                      d = depths[1],
	                      outdir = outdir)
	ebl100 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                       d = depths[2],
	                       outdir = outdir)
	ebl250 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                       d = depths[3],
	                       outdir = outdir)
	ebl500 <- vh.exonblack(blacklist$Q1, runtype = runtype,
	                       d = depths[4],
	                       outdir = outdir)
	message("# Done.\n")
}

vh.covgen <- function(covdir, runtype = "snv",
                      depths = c("50x", "100x", "250x", "500x"),
                      suffix = "bed") {
	
	message("\n# Creating coverage diagnostics paths ...")
	if (runtype == "auto") {
		runtype <- strsplit(covdir, "/")[[1]]
		runtype <- tolower(runtype[length(runtype)])
	}
	if (!(runtype %in% c("snv", "cnv", "rna"))) stop("invalid run type.")
	outdir <- paste0(covdir, "/CoverageDiagnostics")
	dir.create(outdir)
	message("# Done.\n")
	
	message("# Exon coverage diagnostics ...")
	G <- vh.genecov(covdir, outdir = outdir, format = suffix)
	message("# Done.\n")
	
	message("# Extracting blacklisted regions ...")
	ebl50 <- vh.exonblack(G$Q1, runtype = runtype,
	                      d = depths[1],
	                      outdir = outdir)
	ebl100 <- vh.exonblack(G$Q1, runtype = runtype,
	                       d = depths[2],
	                       outdir = outdir)
	ebl250 <- vh.exonblack(G$Q1, runtype = runtype,
	                       d = depths[3],
	                       outdir = outdir)
	ebl500 <- vh.exonblack(G$Q1, runtype = runtype,
	                       d = depths[4],
	                       outdir = outdir)
	message("# Done.\n")
}
