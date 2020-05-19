##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Generates supplementary figures S10A-C and accompanying tables. Writes output to a directory 
# specified on the command line, which is created if it doesn't exist.
#
# Usage:
# Rscript FigS10_intersecting_CNV_maps.R <output_directory>


library(data.table)
library(logging); basicConfig()
library(readxl)
library(dftdLowCov)
library(colorspace)


#####################
# Handle command line
#####################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("No output directory provided. Usage: Rscript FigS10_intersecting_CNV_maps <output_dir>")
}

if (length(args) > 1) {
    stop("Too many arguments. Usage: Rscript FigS10_intersecting_CNV_maps <output_dir>")
}

outdir <- args[1]
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

######################
# Function definitions
######################
make_coordinates <- function(cnvs) {
    cnvs <- copy(cnvs)
    cnvs[Cancer_type == "DFT1", Category := "DFT1"]
    cnvs[Cancer_type == "DFT2", Category := "DFT2"]
    cnvs[Cancer_type == "Non-DFTD", Category := dftdLowCov::scan_col("(\\d+T\\d?) ", Sample_group, c('s'), c('sample'))]
    cnvs[Cancer_type == "Non-DFTD", Category := paste("NonDFTD", Category)]
    
    coordinates <- cnvs[, .(Order=.I, chr=Chr, start=Start, end=End, Type, CNV_ID_ext, Category)]
    return(coordinates)
}

get_cnvs_from_subset <- function(subset) {
    all_cnvs <- subset[, strsplit(paste(commonCnvIds, collapse=','), ',')][[1]]
    all_cnvs
}

select_cnvs_from_table <- function(cnvs, table) {
    table[CNV_ID_ext %in% cnvs]
}

prepare_for_chromosome_map <- function(selected_cnvs) {
    dt <- copy(selected_cnvs)
    setnames(dt, old = c("chr", "Type", "CNV_ID_ext"),
             new = c("seqnames", "State", "Unique.ID"))
    dt[, width := end - start + 1]
    dt
}

simple_chromosome_map <- function(chr_lengths, data, category_colours, ...) {
    darker_colours <- as.list(colorspace::darken(category_colours, 0.3, method = "absolute"))
    names(darker_colours) <- names(category_colours)
    print(darker_colours)
    
    # (Half) plotting height of chromosome
    ch = 1
    # Layer height
    lh = 1
    # Layer separation
    sep = 0.1
    
    # plot y range
    maxy <- data[State=="gain", max(layer)] * lh + ch
    miny <- -(data[State=="loss", max(layer)] * lh + ch)
    
    # plot x range
    minx <- chr_lengths[, min(start)]
    maxx <- chr_lengths[, max(end)]
    
    # Empty plot
    plot(c(minx, maxx), c(miny, maxy), type = "n", xaxt = "n", yaxt = "n",
         bty = "n", xlab = NA, ylab = NA, ...)
    
    # Draw chromosomes
    rect(chr_lengths$start, -ch+sep, chr_lengths$end, ch-sep, col = "grey")
    text(chr_lengths$start+2e6, 0, labels = paste0("Chr", chr_lengths$seqnames), adj=0)
    
    # Draw gain segments
    colours <- data[State == "gain", unlist(category_colours[Category])]
    borders <- data[State == "gain", unlist(darker_colours[Category])]
    data[State == "gain",
         rect(start, (layer-1)*lh+ch+sep, end, layer*lh+ch-sep,
              col = colours, border = borders, lwd = 0.5)]
    
    # Draw loss segments
    colours <- data[State == "loss", unlist(category_colours[Category])]
    borders <- data[State == "loss", unlist(darker_colours[Category])]
    data[State == "loss",
         rect(start, -((layer-1)*lh+ch+sep), end, -(layer*lh+ch-sep),
              col = colours, border = borders, lwd = 0.5)]
}

#################
# Data processing
#################
# Prep
cnvtable <- dftdLowCov::load_cnv_table()

# Chr lengths
chr_lengths <- fread("~/Documents/projects/DFTD_low_coverage/dftd_chromosome_maps/devil_chrom_lengths.tsv")
chr_lengths[, chr := toupper(sub("Chr", "", CHROM))]

coordinates <- make_coordinates(cnvtable)

# Collapse 106T1 106T2 and 106T3 into a single category (maybe unnecessary as they always co-occur)
coordinates[Category %like% "106", Category := sub("106T\\d+", "106Tx", Category)]

segmatch <- make_segmatch(coordinates, ignore_direction = TRUE)

shared_segments <- segmatch[shared_breakpoint == "C_SEGMENT"][order(chr, start, end)]

#################
# Collect subsets
#################
subset_dft1_dft2 <- shared_segments[dft1Count>0 & dft2Count>0 & nonDftdCount==0]
subset_dft1_nondftd <- shared_segments[dft1Count>0 & dft2Count==0 & nonDftdCount>0]
subset_dft2_nondftd <- shared_segments[dft1Count==0 & dft2Count>0 & nonDftdCount>0]
subset_dft1_dft2_nondftd <- shared_segments[dft1Count>0 & dft2Count>0 & nonDftdCount>0]

######################
# Draw chromosome maps
######################
# Prepare tables for plotting
chrlengths <- load_chromosome_lengths()$truncated
chrlengths[, c("start", "end") := .(1, LENGTH)]
for (i in 2:nrow(chrlengths)) {
    chrlengths[i, c("start", "end") := .(start + chrlengths[i-1, end], end + chrlengths[i-1, end])]
}
chrlengths[, seqnames := toupper(sub("Chr", "", CHROM))]

cc <- list(
    "DFT1"="#6597E9",
    "DFT2"="#FA0D11",
    "NonDFTD 340T"="#95DB01",
    "NonDFTD 349T1"="#79AE0B",
    "NonDFTD 352T1"="#628D08",
    "NonDFTD 359T1"="#436300",
    "NonDFTD 435T"="#EBC102",
    "NonDFTD 997T1"="#BA9B0B",
    "NonDFTD 106Tx"="#977D08",
    "NonDFTD 131T"="#695700")
    
filenames <- c("FigS10A_DFT1_and_DFT2",
               "FigS10B_DFT1_and_NonDFTD",
               "FigS10C_DFT2_and_NonDFTD")

subsets <- list(subset_dft1_dft2,
             subset_dft1_nondftd,
             subset_dft2_nondftd)

for (i in 1:3) {
    cnvs <- get_cnvs_from_subset(subsets[[i]])
    selected_rows <- select_cnvs_from_table(cnvs, coordinates)
    plot_data <- prepare_for_chromosome_map(selected_rows)
    plot_data[chrlengths, c("start", "end", "start_within_chr", "end_within_chr") := .(start + i.start - 1, end + i.start - 2, start, end), on = "seqnames"]
    plot_data <- rbindlist(list(rectangle_packing4(plot_data[State == "gain"]),
                                rectangle_packing4(plot_data[State == "loss"])))
    fwrite(plot_data, file.path(outdir, paste0(filenames[i], ".csv")))
    pdf(file.path(outdir, paste0(filenames[i], ".pdf")), width = 30, height = 6)
    simple_chromosome_map(chrlengths, plot_data, cc, ylim = c(-29, 11))
    lnames <- sort(intersect(names(cc), plot_data$Category))
    lfill <- cc[lnames]
    legend("bottomright", legend = lnames, fill = unlist(lfill))
    
    dev.off()
}
