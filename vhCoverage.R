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

vh.reshape <- function(covdata, level, depths, init = 9) {
	
	if (level == "sample") {
		level <- "PatientID"
		a <- ncol(covdata)
		R <- reshape2::melt(covdata[, c(7:(a-5), a)], id = "Depth")
		names(R) <- c("Depth", level, "Coverage")
	
	} else if (level == "gene") {
		covdata <- cov.aggregate(covdata, init = init, t = depths,
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
                     depths = c("5x", "10x", "50x", "100x", "250x", "500x"),
                     format = "bed", plot.name = NULL,
                     x.label = "Percent of covered regions",
                     y.label = "Depth of coverage",
                     w = 20, h = 10, u = 'in', r = 450, cairo = TRUE,
                     outdir = NULL) {
	
	W <- cov.data(directory = covdir, include = include, cleanup = cleanup,
	              runtype = runtype, format = format)
		
	W <- cov.aggregate(W, t = depths, xlab = x.label, ylab = y.label)
	
	if (!is.null(plot.name) & !is.null(outdir)) {
		plot.name = paste0(outdir, "/TSO500_", toupper(runtype), "_", plot.name)
		vh.bar(W$barplot, plot.name, w, h, u, r, cairo = cairo)
	}
	return(W)
}

set.status <- function(Q, depths = c("5x", "10x", "50x", "100x",
                                     "250x", "500x")) {
	d <- length(depths)
	ctrl <- 1:d + (5*d + 1)
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

vh.qmin <- function(M, level, depths) {
	
	Qmin <- aggregate(M[M$Depth == depths[1], 3],
		              list(M[M$Depth == depths[1], 2]), min)
	names(Qmin) <- c(level, "Qmin50")
	
	for (j in 1:length(depths)) {
		Qmin$Qmin100 <- aggregate(M[M$Depth == depths[2], 3],
	                              list(M[M$Depth == depths[2], 2]),
	                              min)[, 2]
		Qmin$Qmin250 <- aggregate(M[M$Depth == depths[3], 3],
	                              list(M[M$Depth == depths[3], 2]),
	                              min)[, 2]
		Qmin$Qmin500 <- aggregate(M[M$Depth == depths[4], 3],
	                              list(M[M$Depth == depths[4], 2]),
	                              min)[, 2]
	}
}

vh.covtable <- function(M, level, depths) {
	
	if (level == "sample") {
		level <- "PatientID"
	} else if (level == "gene") {
		level <- "Gene"
	} else {
		stop("invalid level choice.")
	}
	
	Qmin <- aggregate(M[M$Depth == depths[1], 3],
	                  list(M[M$Depth == depths[1], 2]), min)
	names(Qmin) <- c(level, "Qmin5")
	Qmin$Qmin10 <- aggregate(M[M$Depth == depths[2], 3],
	                          list(M[M$Depth == depths[2], 2]),
	                          min)[, 2]
	Qmin$Qmin50 <- aggregate(M[M$Depth == depths[3], 3],
	                         list(M[M$Depth == depths[3], 2]),
	                         min)[, 2]
	Qmin$Qmin100 <- aggregate(M[M$Depth == depths[4], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          min)[, 2]
	Qmin$Qmin250 <- aggregate(M[M$Depth == depths[5], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          min)[, 2]
	Qmin$Qmin500 <- aggregate(M[M$Depth == depths[6], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          min)[, 2]
	
	Q1 <- aggregate(M[M$Depth == depths[1], 3],
	                list(M[M$Depth == depths[1], 2]),
	                function(x) quantile(x)[2])
	names(Q1) <- c(level, "Q1x5")
	Q1$Q1x10 <- aggregate(M[M$Depth == depths[2], 3],
	                       list(M[M$Depth == depths[2], 2]),
	                       function(x) quantile(x)[2])[, 2]
	Q1$Q1x50 <- aggregate(M[M$Depth == depths[3], 3],
	                       list(M[M$Depth == depths[3], 2]),
	                       function(x) quantile(x)[2])[, 2]
	Q1$Q1x100 <- aggregate(M[M$Depth == depths[4], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[2])[, 2]
	Q1$Q1x250 <- aggregate(M[M$Depth == depths[5], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[2])[, 2]
	Q1$Q1x500 <- aggregate(M[M$Depth == depths[6], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[2])[, 2]
	
	Q2 <- aggregate(M[M$Depth == depths[1], 3],
	                list(M[M$Depth == depths[1], 2]),
	                function(x) quantile(x)[3])
	names(Q2) <- c(level, "Q2x5")
	Q2$Q2x10 <- aggregate(M[M$Depth == depths[2], 3],
	                       list(M[M$Depth == depths[2], 2]),
	                       function(x) quantile(x)[3])[, 2]
	Q2$Q2x50 <- aggregate(M[M$Depth == depths[3], 3],
	                       list(M[M$Depth == depths[3], 2]),
	                       function(x) quantile(x)[3])[, 2]
	Q2$Q2x100 <- aggregate(M[M$Depth == depths[4], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[3])[, 2]
	Q2$Q2x250 <- aggregate(M[M$Depth == depths[5], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[3])[, 2]
	Q2$Q2x500 <- aggregate(M[M$Depth == depths[6], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[3])[, 2]
	
	Q3 <- aggregate(M[M$Depth == depths[1], 3],
	                list(M[M$Depth == depths[1], 2]),
	                function(x) quantile(x)[4])
	names(Q3) <- c(level, "Q3x5")
	Q3$Q3x10 <- aggregate(M[M$Depth == depths[2], 3],
	                       list(M[M$Depth == depths[2], 2]),
	                       function(x) quantile(x)[4])[, 2]
	Q3$Q3x50 <- aggregate(M[M$Depth == depths[3], 3],
	                       list(M[M$Depth == depths[3], 2]),
	                       function(x) quantile(x)[4])[, 2]
	Q3$Q3x100 <- aggregate(M[M$Depth == depths[4], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[4])[, 2]
	Q3$Q3x250 <- aggregate(M[M$Depth == depths[5], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[4])[, 2]
	Q3$Q3x500 <- aggregate(M[M$Depth == depths[6], 3],
	                       list(M[M$Depth == depths[4], 2]),
	                       function(x) quantile(x)[4])[, 2]
	
	Qmax <- aggregate(M[M$Depth == depths[1], 3],
	                  list(M[M$Depth == depths[1], 2]), max)
	names(Qmax) <- c(level, "Qmax5")
	Qmax$Qmax10 <- aggregate(M[M$Depth == depths[2], 3],
	                          list(M[M$Depth == depths[2], 2]),
	                          max)[, 2]
	Qmax$Qmax50 <- aggregate(M[M$Depth == depths[3], 3],
	                          list(M[M$Depth == depths[3], 2]),
	                          max)[, 2]
	Qmax$Qmax100 <- aggregate(M[M$Depth == depths[4], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          max)[, 2]
	Qmax$Qmax250 <- aggregate(M[M$Depth == depths[5], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          max)[, 2]
	Qmax$Qmax500 <- aggregate(M[M$Depth == depths[6], 3],
	                          list(M[M$Depth == depths[4], 2]),
	                          max)[, 2]
	
	Q <- data.frame(Qmin[1:2],
	                Q1x5 = Q1$Q1x5, Q2x5 = Q2$Q2x5, Q3x5 = Q3$Q3x5,
	                Qmax5 = Qmax$Qmax5,
	                Qmin10 = Qmin$Qmin10, Q1x10 = Q1$Q1x10,
	                Q2x10 = Q2$Q2x10, Q3x10 = Q3$Q3x10,
	                Qmax10 = Qmax$Qmax10,
	                Qmin50 = Qmin$Qmin50, Q1x50 = Q1$Q1x50,
	                Q2x50 = Q2$Q2x50, Q3x50 = Q3$Q3x50,
	                Qmax50 = Qmax$Qmax50,
	                Qmin100 = Qmin$Qmin100, Q1x100 = Q1$Q1x100,
	                Q2x100 = Q2$Q2x100, Q3x100 = Q3$Q3x100,
	                Qmax100 = Qmax$Qmax100,
	                Qmin250 = Qmin$Qmin250, Q1x250 = Q1$Q1x250,
	                Q2x250 = Q2$Q2x250, Q3x250 = Q3$Q3x250,
	                Qmax250 = Qmax$Qmax250,
	                Qmin500 = Qmin$Qmin500,
	                Q1x500 = Q1$Q1x500, Q2x500 = Q2$Q2x500,
	                Q3x500 = Q3$Q3x500, Qmax500 = Qmax$Qmax500)
	
	Q$FLAGx5 <- ifelse(Q$Q1x5 >= 75, 0, 1)
	Q$FLAGx5[Q$Q3x5 == 0] <- 2
	Q$FLAGx10 <- ifelse(Q$Q1x10 >= 75, 0, 1)
	Q$FLAGx10[Q$Q3x10 == 0] <- 2
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
		             Q$FLAGx10, Q$FLAGx5, Q$PatientID),]
	} else if (level == "Gene") {
		Q <- Q[order(Q$FLAGx500, Q$FLAGx250, Q$FLAGx100, Q$FLAGx50,
		             Q$FLAGx10, Q$FLAGx5, Q$Gene),]
	}
	
	Q$status <- set.status(Q)
	d <- length(depths)
	flags <- 1:d + (5*d + 1)
	flags <- c(flags, flags[length(flags)] + 1)
	Q <- Q[, c(1, flags, 2:(flags[1] - 1))]
	
	return(Q)
}

vh.genecov <- function(covdir, p = 75, d = "500x",
                       depths = c("5x", "10x", "50x", "100x", "250x", "500x"),
                       runtype = "snv", format = "bed",
                       w = 20, h = 10, u = 'in', r = 450,
                       color = c("lavender", "azure", "lightblue",
                                 "green3", "gold", "darkorange"),
                       line.color = c("orchid", "dodgerblue", "darkblue",
                                      "darkgreen", "brown", "darkred"),
                       x.lab = "Gene",
                       y.lab = "Exon Median Percent Coverage (MPC)",
                       outdir = "", init = 9) {
	
	x <- cov.data(directory = covdir, include = NULL,
	              cleanup = "Exon",
	              runtype = "snv",
	              format = format)
	
	covdata <- cov.aggregate(x, t = depths, init = 7, barplot = FALSE,
	                         percent = TRUE, xlab = x.lab, ylab = y.lab,
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
	
	M4 <- vh.reorder(R, depth = depths[4], level = "gene")[, c(3, 1, 2, 4)]
	rDataset4 <- as.factor(M4$rDataset[M4$Depth == depths[4]])
	Depth4 <- as.factor(M4$Depth[M4$Depth == depths[4]])

	M5 <- vh.reorder(R, depth = depths[5], level = "gene")[, c(3, 1, 2, 4)]
	rDataset5 <- as.factor(M5$rDataset[M5$Depth == depths[5]])
	Depth5 <- as.factor(M5$Depth[M5$Depth == depths[5]])

	M <- vh.reorder(R, depth = depths[6], level = "gene")[, c(3, 1, 2, 4)]
	rDataset <- as.factor(M$rDataset[M$Depth == depths[6]])
	Depth <- as.factor(M$Depth[M$Depth == depths[6]])
	
	B <- list(depth1 = list(M = M1, rDataset = rDataset1, Depth = Depth1),
	          depth2 = list(M = M2, rDataset = rDataset2, Depth = Depth2),
	          depth3 = list(M = M3, rDataset = rDataset3, Depth = Depth3),
	          depth4 = list(M = M4, rDataset = rDataset4, Depth = Depth4),
	          depth5 = list(M = M5, rDataset = rDataset5, Depth = Depth5),
	          ref = list(M = M, rDataset = rDataset, Depth = Depth))
	
	if (!is.null(runtype)) {
		if (outdir != "") outdir <- paste0(outdir, "/")
		for (j in 1:6) {
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
                         depths = c("5x", "10x", "50x", "100x", "250x", "500x"),
                         runtype = "snv", w = 20, h = 10, u = 'in', r = 450,
                         color = c("lavender", "azure", "lightblue",
                                   "green3", "gold", "darkorange"),
                         line.color = c("orchid", "dodgerblue", "darkblue",
                                        "darkgreen", "brown", "darkred"),
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
	
	M4 <- vh.reorder(R, depth = depths[4], level = "gene")[, c(3, 1, 2, 4)]
	rDataset4 <- as.factor(M4$rDataset[M4$Depth == depths[4]])
	Depth4 <- as.factor(M4$Depth[M4$Depth == depths[4]])

	M5 <- vh.reorder(R, depth = depths[5], level = "gene")[, c(3, 1, 2, 4)]
	rDataset5 <- as.factor(M5$rDataset[M5$Depth == depths[5]])
	Depth5 <- as.factor(M5$Depth[M5$Depth == depths[5]])
	
	M <- vh.reorder(R, depth = depths[6], level = "gene")[, c(3, 1, 2, 4)]
	rDataset <- as.factor(M$rDataset[M$Depth == depths[6]])
	Depth <- as.factor(M$Depth[M$Depth == depths[6]])
	
	B <- list(depth1 = list(M = M1, rDataset = rDataset1, Depth = Depth1),
	          depth2 = list(M = M2, rDataset = rDataset2, Depth = Depth2),
	          depth3 = list(M = M3, rDataset = rDataset3, Depth = Depth3),
	          depth4 = list(M = M4, rDataset = rDataset4, Depth = Depth4),
	          depth5 = list(M = M5, rDataset = rDataset5, Depth = Depth5),
	          ref = list(M = M, rDataset = rDataset, Depth = Depth))
	
	if (!is.null(runtype)) {
		if (outdir != "") outdir <- paste0(outdir, "/")
		for (j in 1:6) {
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
	
	Q1 <- cov.q1(covdata, init = init, t = depths, barplot = FALSE)$covdata
	
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
                      depths = c("5x", "10x", "50x", "100x", "250x", "500x"),
                      suffix = "bed") {
	
	message("\n# Creating coverage diagnostics paths ...")
	if (runtype == "auto") {
		runtype <- strsplit(covdir, "/")[[1]]
		runtype <- tolower(runtype[length(runtype)])
	}
	if (!(runtype %in% c("snv", "cnv", "rna", "fusion"))) {
		stop("invalid run type.")
	}
	
	outdir <- paste0(covdir, "/CoverageDiagnostics")
	if (!file.exists(outdir)) {
		dir.create(outdir)
	}
	message("# Done.\n")
	
	message("# Whole-run diagnostics ...")
	x <- vh.panel(covdir, runtype = runtype, format = suffix,
				  plot.name = "wholePanel_exonCoverage.png",
				  y.label = "Percent of covered exons",
				  outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[1], " ..."))
	x5 <- vh.sample(x$covdata, runtype = runtype, depth = depths[1],
					color = "lavender", line.color = "orchid",
					outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[2], " ..."))
	x10 <- vh.sample(x$covdata, runtype = runtype, depth = depths[2],
					 color = "azure", line.color = "dodgerblue",
					 outdir = outdir)
	message("# Done.\n")

	message(paste0("# Sample-level diagnostics at ", depths[3], " ..."))
	x50 <- vh.sample(x$covdata, runtype = runtype, depth = depths[3],
					 color = "lightblue", line.color = "darkblue",
					 outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[4], " ..."))
	x100 <- vh.sample(x$covdata, runtype = runtype, depth = depths[4],
					  color = "green3", line.color = "darkgreen",
					  outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[5], " ..."))
	x250 <- vh.sample(x$covdata, runtype = runtype, depth = depths[5],
					  color = "gold", line.color = "brown",
					  outdir = outdir)
	message("# Done.\n")
	
	message(paste0("# Sample-level diagnostics at ", depths[6], " ..."))
	x500 <- vh.sample(x$covdata, runtype = runtype, depth = depths[6],
					  color = "darkorange", line.color = "darkred",
					  outdir = outdir)
	message("# Done.\n")
	
	message("# Sample-level report preparation ...")
	Q <- vh.covtable(x500$M, level = "sample", depths = depths)
	message("# Done.\n")
	
	message("# Searching for blacklisted regions ...")
	blacklist <- vh.blacklist(Q, x$covdata, runtype = runtype, outdir = outdir)
	B <- vh.covtable(blacklist$B$ref$M, level = "gene", depths = depths)
	message("# Done.\n")
	
	message("# Gene-level reports preparation ...")
	
	vh.covreport(Q, B, runtype = runtype, outdir = outdir)
	
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
	message("# Done.\n")
}

vh.covgen <- function(covdir, runtype = "snv",
                      depths = c("5x", "10x", "50x", "100x", "250x", "500x"),
                      suffix = "bed") {
	
	message("\n# Creating coverage diagnostics paths ...")
	if (runtype == "auto") {
		runtype <- strsplit(covdir, "/")[[1]]
		runtype <- tolower(runtype[length(runtype)])
	}
	if (!(runtype %in% c("snv", "cnv", "rna", "fusion"))) {
		stop("invalid run type.")
	}
	outdir <- paste0(covdir, "/CoverageDiagnostics")
	if (!file.exists(outdir)) {
		dir.create(outdir)
	}
	message("# Done.\n")
	
	message("# Exon coverage diagnostics ...")
	G <- vh.genecov(covdir, outdir = outdir, format = suffix)
	message("# Done.\n")
	
	message("# Extracting blacklisted regions ...")
	ebl5 <- vh.exonblack(G$Q1, runtype = runtype,
	                      d = depths[1],
	                      outdir = outdir)
	ebl10 <- vh.exonblack(G$Q1, runtype = runtype,
	                      d = depths[2],
	                      outdir = outdir)
	ebl50 <- vh.exonblack(G$Q1, runtype = runtype,
	                      d = depths[3],
	                      outdir = outdir)
	ebl100 <- vh.exonblack(G$Q1, runtype = runtype,
	                       d = depths[4],
	                       outdir = outdir)
	ebl250 <- vh.exonblack(G$Q1, runtype = runtype,
	                       d = depths[5],
	                       outdir = outdir)
	ebl500 <- vh.exonblack(G$Q1, runtype = runtype,
	                       d = depths[6],
	                       outdir = outdir)
	message("# Done.\n")
}
