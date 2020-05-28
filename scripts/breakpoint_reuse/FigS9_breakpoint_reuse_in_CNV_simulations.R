##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Plots the main figures 2G and 2I, and the supplementary figures S9A and S9B
#
# Usage:

library(data.table)
library(openxlsx)
library(dftdLowCov)
library(stringr)

#' Generates breakpoint reuse tables for each simulation replicate
process_sims <- function(sims, origtable) {
    POS.SENTINEL <- .Machine$integer.max
    NEG.SENTINEL <- -POS.SENTINEL
    
    # Set up a data.table to map New.CNV_ID_ext to clone type
    dt <- origtable[, .(CNV_ID_ext, Cancer_type, Sample_group)]
    dt[Cancer_type == "DFT1", Category := "DFT1"]
    dt[Cancer_type == "DFT2", Category := "DFT2"]
    dt[Cancer_type == "Non-DFTD", Category := dftdLowCov::scan_col("(\\d+T\\d?) ", Sample_group, c('s'), c('sample'))]
    dt[Cancer_type == "Non-DFTD", Category := paste("NonDFTD", Category)]
    dt[, Sample_group := NULL]
    
    # Reconstitute original New.CNV_ID_ext from Unique.ID (which may have an appendix 'a' or 'b'
    # if the sim CNV has been wrapped)
    sims[, CNV_ID_ext := sub("*\\.[ab]", "\\1", Unique.ID)]
    
    # Add in the category according to the dt
    sims[dt, Category := Category, on = c("CNV_ID_ext")]
    
    # Round starts and ends to nearest 100kb
    sims[, c("start", "end") := .(as.integer(round(as.integer(start), -5)), as.integer(round(as.integer(end), -5)))]
    sims[, width := end - start]
    
    # Breakpoint reuse
    # Do some tidying...
    setnames(sims, old = c("seqnames", "State", "Unique.ID"), new = c("chr", "Type", "CNV_ID_ext"))
    sims[, Order := 1:.N, by = rep]
    sims <- sims[width > 0]
    
    # Deal with wrapped CNVs. I only want to count the first start and last end points among
    # wrapped subCNVs, not the breakpoints that are induced by chromosome starts and ends.
    # I will set these to deliberately unreasonable values, and filter them out later.
    POS.SENTINEL <- .Machine$integer.max
    NEG.SENTINEL <- -POS.SENTINEL
    
    # Easy to spot starts...
    sims[(wrapped) & start == 0, start := NEG.SENTINEL]
    # Harder to spot ends...
    chrends <- sims[, .(chrend = max(end)), by = chr][order(chr)]
    for (chrom in chrends[, chr]) {
        chrend <- chrends[chr == chrom, chrend]
        sims[(wrapped) & chr == chrom & end >= (chrend - 100000), end := POS.SENTINEL]
    }
    
    sims
}


#' Make segment match table for each simulation replicate
make_segmatches <- function(sims) {
    POS.SENTINEL <- .Machine$integer.max
    NEG.SENTINEL <- -POS.SENTINEL
    
    segmatches <- list()
    replicate_indices <- sort(unique(sims[, rep]))
    for (replicate in replicate_indices) {
        sm <- make_segmatch(sims[rep==replicate])
        sm[, keep := TRUE]
        sm[shared_breakpoint == "A_START_ONLY" & start == NEG.SENTINEL, keep := FALSE]
        sm[shared_breakpoint == "B_END_ONLY" & end == POS.SENTINEL, keep := FALSE]
        sm[shared_breakpoint == "C_SEGMENT" & (start == NEG.SENTINEL | end == POS.SENTINEL), keep := FALSE]
        sm <- sm[(keep)]
        sm[, rep := replicate]
        segmatches[[replicate]] <- sm
        print(replicate)
    }
    segmatches <- rbindlist(segmatches)
    segmatches <- segmatches[dft1Count + dft2Count + nonDftdCount > 1]
    
    segmatches[, multipleDFT1 := dft1Count > 1]
    segmatches[, multipleDFT2 := dft2Count > 1]
    segmatches[, multipleNonDFT := nonDftdCount > 1]
    segmatches[, anyDFT1 := grepl("DFT1", tumourTypes)]
    segmatches[, anyDFT2 := grepl("DFT2", tumourTypes)]
    segmatches[, anyNonDFTD := grepl("NonDFTD", tumourTypes)]
    segmatches
}


# Load original data
cnvtable <- dftdLowCov::load_cnv_table()

# Set some filepaths (TODO - get from command line)
dft1_simfile <- "/Users/kg8/Documents/projects/DFTD_low_coverage/dftd_breakpoint_reuse/simulation_results/DFT1_breakpoint_reuse_set/samples.tsv"
allclones_simfile <- "/Users/kg8/Documents/projects/DFTD_low_coverage/dftd_breakpoint_reuse/simulation_results/All_clones_breakpoint_reuse_set/samples.tsv"

#####################
# Process simulations
#####################
# Load previous simulations
sims_dft1 <- process_sims(
    fread(dft1_simfile),
    cnvtable)

sims_all <- process_sims(
    fread(allclones_simfile),
    cnvtable)

segmatches_dft1 <- make_segmatches(sims_dft1)
fwrite(segmatches_dft1, file.path(outdir, "FigS9_breakpoint_reuse_in_CNV_simulations_DFT1.csv"))

segmatches_all <- make_segmatches(sims_all)
fwrite(segmatches_all, file.path(outdir, "FigS9_breakpoint_reuse_in_CNV_simulations_all.csv"))

segmatches_all[(multipleDFT1) & !(multipleDFT2) & !(multipleNonDFT), rowlab := "DFT1"]
segmatches_all[(multipleDFT1) & !(multipleDFT2) & (multipleNonDFT), rowlab := "D1ND"]
segmatches_all[(multipleDFT1) & (multipleDFT2) & !(multipleNonDFT), rowlab := "D1D2"]
segmatches_all[(multipleDFT1) & (multipleDFT2) & (multipleNonDFT), rowlab := "D12N"]
segmatches_all[!(multipleDFT1) & !(multipleDFT2) & !(multipleNonDFT), rowlab := "SNGL"]
segmatches_all[!(multipleDFT1) & !(multipleDFT2) & (multipleNonDFT), rowlab := "NDFT"]
segmatches_all[!(multipleDFT1) & (multipleDFT2) & !(multipleNonDFT), rowlab := "DFT2"]
segmatches_all[!(multipleDFT1) & (multipleDFT2) & (multipleNonDFT), rowlab := "D2ND"]

segmatches_all[, collab := NA_character_]
segmatches_all[(anyDFT1) & (anyDFT2) & !(anyNonDFTD), collab := "dft1_2"]
segmatches_all[(anyDFT1) & !(anyDFT2) & (anyNonDFTD), collab := "dft1_N"]
segmatches_all[!(anyDFT1) & (anyDFT2) & (anyNonDFTD), collab := "dft2_N"]
segmatches_all[str_count(tumourTypes, "NonDFTD") > 1 & !(anyDFT1) & !(anyDFT2), collab := "dftN_N"]
segmatches_all[(anyDFT1) & (anyDFT2) & (anyNonDFTD), collab := "dft12N"]

get_totals <- function(data) {
    total <- data[, .N]
    multiple_dft1_only <- data[(multipleDFT1) & !(multipleDFT2) & !(multipleNonDFT), .N]
    multiple_dft1_and_dft2 <- data[(multipleDFT1) & (multipleDFT2) & !(multipleNonDFT), .N]
    multiple_dft2_only <- data[!(multipleDFT1) & (multipleDFT2) & !(multipleNonDFT), .N]
    multiple_dft2_and_nondft <- data[!(multipleDFT1) & (multipleDFT2) & (multipleNonDFT), .N]
    multiple_nondft_only <- data[!(multipleDFT1) & !(multipleDFT2) & (multipleNonDFT), .N]
    multiple_dft1_and_nondft <- data[(multipleDFT1) & !(multipleDFT2) & (multipleNonDFT), .N]
    multiple_everything <- data[(multipleDFT1) & (multipleDFT2) & (multipleNonDFT), .N]
    single_everything <- data[!(multipleDFT1) & !(multipleDFT2) & !(multipleNonDFT), .N]
    output <- 
        c(`DFT1` = multiple_dft1_only,
          `DFT2` = multiple_dft2_only,
          `NDFT` = multiple_nondft_only,
          `D1D2` = multiple_dft1_and_dft2,
          `D1ND` = multiple_dft1_and_nondft,
          `D2ND` = multiple_dft2_and_nondft,
          `D12N` = multiple_everything,
          `SNGL` = single_everything)
    stopifnot(total == sum(output))
    output
}

collect_totals <- function(segmatches) {
    # DFT1 and DFT2
    dft1_2 <- get_totals(segmatches[nTumourTypes>1][grepl("DFT1", tumourTypes) & grepl("DFT2", tumourTypes) & !(grepl("NonDFTD", tumourTypes))])
    
    # DFT1-NonDFTD
    dft1_N <- get_totals(segmatches[nTumourTypes>1][grepl("DFT1", tumourTypes) & !(grepl("DFT2", tumourTypes)) & grepl("NonDFTD", tumourTypes)])
    
    # DFT2-NonDFTD
    dft2_N <- get_totals(segmatches[nTumourTypes>1][!(grepl("DFT1", tumourTypes)) & grepl("DFT2", tumourTypes) & grepl("NonDFTD", tumourTypes)])
    
    # Non-DFTDs
    dftN_N <- get_totals(segmatches[nTumourTypes>1][!(grepl("DFT1", tumourTypes)) & !grepl("DFT2", tumourTypes) & grepl("NonDFTD", tumourTypes)])
    
    # All clones
    dft12N <- get_totals(segmatches[nTumourTypes>1][(grepl("DFT1", tumourTypes)) & grepl("DFT2", tumourTypes) & grepl("NonDFTD", tumourTypes)])
    
    data <- cbind(dft1_2, dft1_N, dft2_N, dftN_N, dft12N)
    data
}

mycolours <- list(
    DFT1 = "#1887BE",
    DFT2 = "#D1403B",
    NDFT = "#20A130", #"#00EB70",
    SNGL = "grey80",
    D1D2 = "#94698C",
    D1ND = "#4FA8A9", # "#10BF9B",
    D2ND = "#E96C29", #"#93AC59",
    D12N = "#BE7672" #"#79A083"
)

NSIMS <- sims_all[, max(rep)]

##############################################################
# Plot main figures 2G and 2I, and supplementary figures S9A-B
##############################################################

setwd("outdir")
pdf(paste0("FigS9B.pdf"), width = 9, height = 8)
par(mar=c(2.5,3,2,2.5), lwd=0.25)
catorder <- c("DFT1", "D1ND", "NDFT", "D1D2", "D2ND", "DFT2", "D12N", "SNGL")
data <- collect_totals(segmatches_all) / NSIMS
barplot(data[catorder,], names.arg = c("DFT1/DFT2", "DFT1/Non-DFTD", "DFT2/Non-DFTD", "Non-DFTD/Non-DFTD", "DFT1/DFT2/Non-DFTD"),
        border = "white", col = simplify2array(mycolours[catorder]),
        space = 1, lwd = 1,  main = NA, cex.names = 0.7, beside=F, ylim = c(0, 20))
dev.off()

# Digression - make segmatches for the real data, too
coordinates <- dftdLowCov::make_coordinates(cnvtable)
segmatch_real <- make_segmatch(coordinates, ignore_direction = TRUE)
segmatch_real[, multipleDFT1 := dft1Count > 1]
segmatch_real[, multipleDFT2 := dft2Count > 1]
segmatch_real[, multipleNonDFT := nonDftdCount > 1]
segmatch_real[, anyDFT1 := grepl("DFT1", tumourTypes)]
segmatch_real[, anyDFT2 := grepl("DFT2", tumourTypes)]
segmatch_real[, anyNonDFTD := grepl("NonDFTD", tumourTypes)]
realdata <- collect_totals(segmatch_real)

pdf("Fig2I.pdf", width = 2.65, height = 2.75)
par(mar=c(4,3,2,2.5))
plot(c(0, 5), c(0, 90), type = "n",
     xaxt = "n", bty = "l",
     xlab = NA, ylab = NA, las = 2, cex.axis = 0.5, tck = -0.02, mgp=c(3, .3, 0))
axis(1, at = 1:5 - 0.5,
     labels = c("DFT1\nDFT2", "DFT1\nNon-DFTD", "DFT2\nNon-DFTD", "Non-DFTD\nNon-DFTD", "DFT1\nDFT2\nNon-DFTD"),
     las = 2, cex.axis = 0.5, tick = F, line=-0.7)
catorder <- c("DFT1", "D1ND", "NDFT", "D1D2", "D2ND", "DFT2", "D12N", "SNGL")
rd <- apply(realdata[catorder,], 2, cumsum)
sd <- apply(simdata[catorder,], 2, cumsum)
for (c in 1:ncol(rd)) {
    rbottom <- c(0, rd[1:(nrow(rd)-1), c])
    rtop <- rd[1:nrow(rd), c]
    rleft <- c - 0.48
    rright <- c - 0.15
    i <- (rtop - rbottom) > 0
    rect(rleft, rbottom[i], rright, rtop[i],
         col = unlist(mycolours[rownames(rd)][i]),
         border = unlist(mycolours[rownames(rd)][i]), lwd = 0.6)
    
    sbottom <- c(0, sd[1:(nrow(sd)-1), c])
    stop <- sd[1:nrow(sd), c]
    sleft <- c - 1 + 0.15
    sright <- c - 0.52
    j <- (stop - sbottom) > 0
    rect(sleft, sbottom[j], sright, stop[j],
         col = unlist(mycolours[rownames(sd)][j]),
         border = unlist(mycolours[rownames(sd)][j]),
         density = 30, angle = 30, lwd = 0.6)
}
dev.off()

# Within - simulations
# Data prep
dt <- segmatches_dft1[dft1Count > 1]
counts <- dt[, .N, by = .(rep, dft1Count)]
tmp <- list()
for (count in counts[, unique(dft1Count)]) {
    present <- counts[dft1Count==count]
    missing <- setdiff(1:NSIMS, present$rep)
    if (length(missing) > 0) {
        absent <- data.table(rep=setdiff(1:NSIMS, present$rep), dft1Count = count, N = 0)
    } else {
        absent <- present[0]
    }
    tmp[[length(tmp) + 1]] <- present
    tmp[[length(tmp) + 1]] <- absent
}
tmp <- rbindlist(tmp)[order(rep, dft1Count)]
stats <- tmp[, .(mean = mean(N), sd = sd(N), stderr = sd(N)/sqrt(.N)), by = dft1Count]
stats[, errlower := mean - sd]#1.96*stderr]
stats[, errupper := mean + sd]#1.96*stderr]

# Grab real data
dt <- segmatch_real[dft1Count > 1]
counts <- dt[, .N, by = dft1Count][order(dft1Count)]
missing_countvals <- counts[, setdiff(seq(2,max(dft1Count)), dft1Count)]

# Custom barplot for Figure 2G
{
    pdf("Figure2G.pdf", width=2.8, height=2.6)
    edge_leeway <- 0.15
    mid_leeway <- 0.025
    data <- stats[counts, , on = "dft1Count"]
    plot(c(1.2, 10), c(0, 80), type = "n",
         xaxt = "n", bty = "l",
         ylab = "Count", xlab = NA,
         las = 2, cex.axis = 0.5, tck = -0.02, mgp=c(0.7, .3, 0), cex.lab = .5)
    axis(1, at = 2:10 - 0.5,
         labels = 2:10,
         las = 1, cex.axis = 0.5, tick = F, line=-1.1)
    title(xlab = "Breakpoint reuse within DFT1", cex.lab = 0.5, line = 0.3)
    for (c in data$dft1Count) {
        left <- c - 1
        midpoint <- c - 0.5
        right <- c
        simheight <- data[dft1Count == c, mean]
        realheight <- data[dft1Count == c, N]
        # Simulation
        if (!is.na(simheight)) {
            rect(left + edge_leeway, 0, midpoint - mid_leeway, simheight,
                 col = scales::alpha(mycolours$DFT1, 0.6), border = mycolours$DFT1,
                 density = 30, angle = 30, lwd = 0.8)
            arrows((midpoint - mid_leeway + left + edge_leeway) / 2,
                   stats[dft1Count == c, errlower],
                   (midpoint - mid_leeway + left + edge_leeway) / 2,
                   stats[dft1Count == c, errupper], angle = 90, code = 3, length = .015, lwd = 0.5)
        } else {
            rect(left + edge_leeway, 0, midpoint - mid_leeway, 0,
                 col = scales::alpha(mycolours$DFT1, 0.6), border = mycolours$DFT1,
                 density = 30, angle = 30, lwd = 0.5)
        }
        
        # Real data
        rect(midpoint + mid_leeway, 0, right - edge_leeway, realheight,
             col = scales::alpha(mycolours$DFT1, 0.6), border = mycolours$DFT1, lwd = 0.8)
    }
    # Legend
    rect(6, 75, 6.325, 81, col = scales::alpha(mycolours$DFT1, 0.6), border = mycolours$DFT1,
         density = 30, angle = 30, lwd = 0.8)
    text(6.6, 78, "Simulated", cex = 0.3, adj=0)
    
    arrows(6.1625, 62, 6.1625, 68,
           angle = 90, code = 3, length = .015, lwd = 0.5)
    text(6.6, 65, "Standard deviation", cex = 0.3, adj=0)
    
    rect(6, 49, 6.325, 55, col = scales::alpha(mycolours$DFT1, 0.6), border = mycolours$DFT1,
         lwd = 0.8)
    text(6.6, 52, "Observed", cex = 0.3, adj=0)
    dev.off()
}


{
    pdf(paste0("FigS9A", timestamp, ".pdf"), width = 210/25.4, height = 297/25.4)
    par(mar=c(5.4, 6.4, 4.4, 3.4))
    dt <- segmatches_dft1[dft1Count > 1]
    counts <- dt[, .N, by = dft1Count][order(dft1Count)]
    missing_countvals <- counts[, setdiff(seq(2,max(dft1Count)), dft1Count)]
    if (length(missing_countvals) > 0) {
        missing <- data.table(dft1Count = counts[, setdiff(seq(2,max(dft1Count)), dft1Count)], N = 0)
        counts <- rbindlist(list(counts, missing))[order(dft1Count)]
    }
    bp <- barplot(counts[, N/NSIMS], names.arg = counts[, dft1Count], col = scales::alpha(mycolours$DFT1, 0.6),
                  border = mycolours$DFT1, main = NA, ylim = c(0, 20), xlab = "Number of CNVs sharing a breakpoint",
                  ylab = "Average frequency per simulation replicate", width = 1, space = 0.21, xaxt = "n", las = 2, cex.axis = 1.5, cex.lab = 2)
    axis(1, at = bp, labels = 2:4, cex.axis = 1.5)
    #text(x = bp, y = -1, labels = 2:max(counts[, 1]), xpd = TRUE, cex = 1.0)
    box(which = "plot", bty = "l")
    dev.off()
}
