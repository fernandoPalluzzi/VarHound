#  VarHound - TSO500 - Diagnostics core functions

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


cov.header <- function(x, re1 = "[xX]*\\d+[xX]*", re2 = "[xX]", prefix = "x") {
	h <- colnames(x)
	h[1] <- str_replace(h[1], "X.", "")
	m <- regexpr(re1, h)
	h <- c(h[m == -1], paste0(prefix, str_replace_all(h[m > 0], re2, "")))
	return(list(header = h, offset = sum(m == -1), n.depths = sum(m > 0)))
}

cov.preprocess <- function(file, include = NULL, cleanup = NULL, runtype = NULL) {
	x <- read.delim(file, stringsAsFactors = FALSE)
	h <- cov.header(x)
	names(x) <- h$header
	x$ID <- 1:nrow(x)
	if (is.null(include)) {
		include <- ""
	} else {
		include <- paste(include, collapse = "|")
	}
	if (!is.null(runtype)) {
		coll <- ifelse(include == "", include, "|")
		if (runtype == "snv") {
			g <- paste0("PIK3CA|EGFR|ALK|ROS1|BRAF|RET|",
						"BRCA1|BRCA2|KIT|PDGFRA|KRAS|NRAS|",
						"NTRK1|NTRK2|NTRK3")
		} else if (runtype == "cnv") {
			g <- paste0("ALK|BRAF|BRCA1|BRCA2|EGFR|KIT|KRAS|NRAS|",
			            "PDGFRA|PIK3CA|RET")
		} else if (runtype == "rna") {
			g <- paste0("ALK|BRAF|EGFR|KRAS|NRAS|NTRK1|NTRK2|NTRK3|",
			            "PIK3CA|RET")
		} else if (runtype == "fusion") {
			g <- paste0("AKT2|ALK|AR|ATM|BRAF|BRCA1|BRCA2|CCND1|CCND3|",
			            "CCNE1|CDK4|CDK6|CHEK1|CHEK2|EGFR|ERBB2|ERBB3|",
			            "ERCC1|ERCC2|ESR1|FGF1|FGF10|FGF14|FGF19|FGF2|",
			            "FGF23|FGF3|FGF4|FGF5|FGF6|FGF7|FGF8|FGF9|",
			            "FGFR1|FGFR2|FGFR3|FGFR4|JAK2|KIT|KRAS|LAMP1|",
			            "MDM2|MDM4|MET|MYC|MYCL1|MYCN|NRAS|NRG1|",
			            "PDGFRA|PDGFRB|PIK3CA|PIK3CB|PTEN|RAF1|RET|",
			            "RICTOR|RPS6KB1|TFRC")
		}
		include <- paste(c(include, g), collapse = coll)
	}
	x <- x[grep(include, x$region),]
	if (!is.null(cleanup)) {
		cleanup <- paste(cleanup, collapse = "|")
		x <- x[grep(cleanup, x$region),]
	}
	return(x)
}

cov.yield <- function(file, yield = "sequential", g = 2, runtype = "snv") {
	
	x <- cov.preprocess(file)
	if (yield == "plainCoverage") {
		W <- x[, 5:8]
		j <- c(1, 1, 1, 2, 3, 4)
		if (runtype %in% c("dna", "snv", "cnv")) {
			labs <- c("x50", "x100", "x250", "x500")
		} else if (runtype %in% c("rna", "fusion")) {
			labs <- c("x5", "x10", "x50", "x100")
		} else {
			stop("invalid runtype choice.")
		}
	} else if (yield == "sequential") {
		if (runtype %in% c("dna", "snv", "cnv")) {
			W <- data.frame(CY100 = -100*(x$x100 - x$x50)/(x$x50 + 1),
							CY250 = -100*(x$x250 - x$x100)/(x$x100 + 1),
							CY500 = -100*(x$x500 - x$x250)/(x$x250 + 1))
			labs <- c("CY100", "CY250", "CY500")
		} else if (runtype %in% c("rna", "fusion")) {
			W <- data.frame(CY10 = -100*(x$x10 - x$x5)/(x$x5 + 1),
							CY50 = -100*(x$x50 - x$x10)/(x$x10 + 1),
							CY100 = -100*(x$x100 - x$x50)/(x$x50 + 1))
			labs <- c("CY10", "CY50", "CY100")
		} else {
			stop("invalid runtype choice.")
		}
		j <- c(1, 2, 1, 2, 3, 3)
	} else if (yield == "reference") {
		if (runtype %in% c("dna", "snv", "cnv")) {
			W <- data.frame(CY100 = -100*(x$x100 - x$x50)/(x$x50 + 1),
							CY250 = -100*(x$x250 - x$x50)/(x$x50 + 1),
							CY500 = -100*(x$x500 - x$x50)/(x$x50 + 1))
			labs <- c("CY100", "CY250", "CY500")
		} else if (runtype %in% c("rna", "fusion")) {
			W <- data.frame(CY10 = -100*(x$x10 - x$x5)/(x$x5 + 1),
							CY50 = -100*(x$x50 - x$x10)/(x$x5 + 1),
							CY100 = -100*(x$x100 - x$x50)/(x$x5 + 1))
			labs <- c("CY10", "CY50", "CY100")
		} else {
			stop("invalid runtype choice.")
		}
		j <- c(1, 2, 1, 2, 3, 3)
	}
	
	if (is.null(g)) {
		Xb <- mclustBIC(W)
	} else {
		Xb <- mclustBIC(W, G = 1:g)
	}
	Xs <- summary(Xb, data = W)
	x$mclust <- Xs$classification
	message("Done.")
	
	message("Plots creation ...")
	png("Coverage_yield_plot1.png", width = 20, height = 10, units = 'in', res = 450)
	coordProj(data = W, dimens = c(labs[j[1]], labs[j[4]]),
	          what = "classification",
	          parameters = Xs$parameters, z = Xs$z)
	if (yield != "plainCoverage") {
		abline(a = 0, b = 1, lty = 2)
		abline(a = 0, b = 2, lty = 2)
		abline(a = 0, b = 4, lty = 2)
	}
	dev.off()
	
	png("Coverage_yield_plot2.png", width = 20, height = 10, units = 'in', res = 450)
	coordProj(data = W, dimens = c(labs[j[2]], labs[j[5]]),
	          what = "classification",
	          parameters = Xs$parameters, z = Xs$z)
	if (yield != "plainCoverage") {
		abline(a = 0, b = 1, lty = 2)
		abline(a = 0, b = 2, lty = 2)
		abline(a = 0, b = 4, lty = 2)
	}
	dev.off()
	
	png("Coverage_yield_plot3.png", width = 20, height = 10, units = 'in', res = 450)
	coordProj(data = W, dimens = c(labs[j[3]], labs[j[6]]),
	          what = "classification",
	          parameters = Xs$parameters, z = Xs$z)
	if (yield != "plainCoverage") {
		abline(a = 0, b = 1, lty = 2)
		abline(a = 0, b = 2, lty = 2)
		abline(a = 0, b = 4, lty = 2)
	}
	dev.off()
	message("Done.")
	
	return(W)
}

cov.long <- function(file, include = NULL, cleanup = NULL, runtype = NULL,
                     varying = NULL) {
	x <- cov.preprocess(file, include, cleanup, runtype)
	if (is.null(varying)) {
		h <- cov.header(x)
		varying <- list(seq(1, h$n.depths) + h$offset - 1)
	}
	L <- reshape(x, idvar = "ID", varying = varying,
	             v.names = "coverage",
	             direction = "long")
	L$length <- L$end - L$start
	L$p <- 100*L$coverage/L$length
	return(L)
}

cov.data <- function(directory, include = NULL, cleanup = NULL,
                     runtype = NULL, format = "bed") {
	check <- TRUE
	for (file in list.files(directory, recursive = TRUE)) {
		if (endsWith(file, format)) {
			sample.name <- strsplit(file, "\\.")[[1]][1]
			file.path <- paste0(directory, "/", file)
			if (check) {
				R <- cov.long(file.path, include, cleanup,
				              runtype)[, c(1:6, 9)]
				check <- FALSE
			} else {
				R <- data.frame(R, p = cov.long(file.path, include,
				                                cleanup, runtype)[, 9])
			}
			sample.name <- basename(sample.name)
			names(R)[ncol(R)] <- sample.name
		}
	}
	return(R)
}

cov.aggregate <- function(x, t = NULL, init = 7, barplot = TRUE,
                          percent = TRUE, xlab = "", ylab = "",
                          dummy = FALSE) {
	
	if (dummy) {
		x$median <- x[, init]
		x$min <- x[, init]
		x$max <- x[, init]
	} else {
		x$median <- apply(x[, init:ncol(x)], 1, median)
		x$min <- apply(x[, init:ncol(x)], 1, min)
		x$max <- apply(x[, init:ncol(x)], 1, max)
	}
	
	x$Coverage <- "0"
	x$Coverage[x$median > 0 & x$median <= 25] <- "(0, 25]"
	x$Coverage[x$median > 25 & x$median <= 50] <- "(25, 50]"
	x$Coverage[x$median > 50 & x$median <= 75] <- "(50, 75]"
	x$Coverage[x$median > 75 & x$median < 100] <- "(75, 100)"
	x$Coverage[x$median == 100] <- "100"
	x$Coverage <- factor(x$Coverage, levels = c("0", "(0, 25]",
	                                            "(25, 50]", "(50, 75]",
	                                            "(75, 100)", "100"))
	
	x$Depth <- "None"
	for (k in 1:length(t)) {
		x$Depth[x$time == k] <- t[k]
	}
	x$Depth <- factor(x$Depth, levels = t)
	
	if (barplot) {
		bar <- ggplot(x, aes(x = Depth, fill = Coverage))
		if (percent) {
			bar <- bar + geom_bar(position = "fill")
		} else {
			bar <- bar + geom_bar(stat = "identity")
		}
		bar <- bar + scale_fill_manual(name = "MPC",
		                               values = c("red2", "orange",
		                                          "gold", "green3",
		                                          "deepskyblue2",
		                                          "royalblue")) +
		theme_bw() +
		theme(panel.border = element_blank(),
		      panel.grid.major = element_line(),
		      panel.grid.minor = element_blank(),
		      axis.line = element_line(colour = "black"),
		      axis.text = element_text(size = 24),
		      axis.title = element_text(size = 26, face = "bold"),
		      legend.key.size = unit(1, "cm"),
		      legend.text = element_text(size = 22),
		      legend.title = element_text(size = 24)) +
		labs(x = xlab, y = ylab) +
		scale_y_continuous(labels = scales::percent)
	} else {
		bar <- NULL
	}
	
	return(list(covdata = x, barplot = bar))
}

cov.q1 <- function(x, t = NULL, init = 7, barplot = TRUE, percent = TRUE,
                   xlab = "", ylab = "") {
	
	x$Q1 <- apply(x[, init:ncol(x)], 1, function(v) quantile(v, probs = 0.25))
	x$min <- apply(x[, init:ncol(x)], 1, min)
	x$max <- apply(x[, init:ncol(x)], 1, max)
	
	x$Coverage <- "0"
	x$Coverage[x$Q1 > 0 & x$Q1 <= 25] <- "(0, 25]"
	x$Coverage[x$Q1 > 25 & x$Q1 <= 50] <- "(25, 50]"
	x$Coverage[x$Q1 > 50 & x$Q1 <= 75] <- "(50, 75]"
	x$Coverage[x$Q1 > 75 & x$Q1 < 100] <- "(75, 100)"
	x$Coverage[x$Q1 == 100] <- "100"
	x$Coverage <- factor(x$Coverage, levels = c("0", "(0, 25]",
	                                            "(25, 50]", "(50, 75]",
	                                            "(75, 100)", "100"))
	
	x$Depth <- "None"
	for (k in 1:length(t)) {
		x$Depth[x$time == k] <- t[k]
	}
	x$Depth <- factor(x$Depth, levels = t)
	
	if (barplot) {
		bar <- ggplot(x, aes(x = Depth, fill = Coverage))
		if (percent) {
			bar <- bar + geom_bar(position = "fill")
		} else {
			bar <- bar + geom_bar(stat = "identity")
		}
		bar <- bar + scale_fill_manual(name = "MPC",
		                               values = c("red2", "orange",
		                                          "gold", "green3",
		                                          "deepskyblue2",
		                                          "royalblue")) +
		theme_bw() +
		theme(panel.border = element_blank(),
		      panel.grid.major = element_line(),
		      panel.grid.minor = element_blank(),
		      axis.line = element_line(colour = "black"),
		      axis.text = element_text(size = 24),
		      axis.title = element_text(size = 26, face = "bold"),
		      legend.key.size = unit(1, "cm"),
		      legend.text = element_text(size = 22),
		      legend.title = element_text(size = 24)) +
		labs(x = xlab, y = ylab) +
		scale_y_continuous(labels = scales::percent)
	} else {
		bar <- NULL
	}
	
	return(list(covdata = x, barplot = bar))
}
