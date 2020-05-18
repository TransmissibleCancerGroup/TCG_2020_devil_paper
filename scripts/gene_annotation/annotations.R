##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Generates supplementary tables S6A-D. Requires several input tables, which are included in the
# repository under 'scripts/gene_annotation/data'. The input data directory is specified on the 
# command line. Writes output to a directory specified on the command line, which is created if
# it doesn't exist.
#
# Usage:
# Rscript annotations.R <data_directory> <output_directory>

library(readxl)
library(openxlsx)
library(data.table)
library(logging); basicConfig();# addHandler(writeToFile, file = "check_totals_log.txt")
library(dftdLowCov)
library(argparse)

#############################################
# Get the data directory from the commandline
#############################################
parser <- ArgumentParser()
parser$add_argument("indir", help = "Input directory containing lookup tables and COSMIC data")
parser$add_argument("outdir", help = "Output directory to write output supplementary tables")

args <- parser$parse_args()
data.dir <- args$indir
outdir <- args$outdir

if (!dir.exists(data.dir)) {
    stop("Provided directory does not exist.")
}

if (!dir.exists(outdir)) {
    dir.create(outdir)
}

##################################################################################
# Read in some lookup tables to match gene names between human and Devil gene sets
##################################################################################

# Define some translation tables
# 'position_lookup' translates between contig naming schemes, and maps contig positions to chromosome positions
# 'opossum_to_devil' maps opossum ensembl IDs (provided in Deakin et al.) to Devil ensembl IDs and contig positions
position_lookup <- fread(file.path(data.dir, "Devil_position_lookup_table.csv"))
opossum_to_devil <- fread(file.path(data.dir, "opossum_to_devil_genes.txt"))
opossum_to_human_and_devil <- fread(file.path(data.dir, "opossum_to_human_and_devil.txt"))
human_to_devil <- fread(file.path(data.dir, "human_to_devil.txt"))
human_to_opossum_and_devil <- fread(file.path(data.dir, "human_to_opossum_and_devil.txt"))
devil_gene_list <- fread(file.path(data.dir, "devil_gene_list.txt"))

# Read in the COSMIC gene list
cosmic <- fread(file.path(data.dir, "Census_allMon_Nov_18_00_12_15_2019.csv"))

# Read in auxiliary Devil gene name translation file
gene_translation <- fread(file.path(data.dir, "Devil_gene_translation.txt"))

#############################################
# Load the CNV data and do some preprocessing
#############################################
cnv_table <- dftdLowCov::load_cnv_table()

# Change column names of CNV table to match the other tables
setnames(cnv_table, old = c("Chr", "Start", "End"), new = c("chr", "start", "end"))

##############################
# Preprocess the lookup tables
##############################
gene_list <- position_lookup[devil_gene_list, on = "CONTIG.ALT", nomatch=0L]
gene_list[, contig_startpos := start]
gene_list[, contig_endpos := end]
gene_list[, start := OFFSET + contig_startpos]
gene_list[, end := OFFSET + contig_endpos]
gene_list[, chr := toupper(sub("Chr", "", CHROM))]

setkey(gene_list, chr, start, end)
setkey(cnv_table, chr, start, end)

########################
# Annotations Dec 2019 #
########################

make_annotation_table <- function(data, filename, required_depth = 3) {
    chrlengths <- data.table(chr=c(1:6, "X"), LENGTH=c(582483498, 598991806, 519046246, 412394343, 249841267, 215187208, 54800000))
    
    # List intervals that have gain coverage of ≥3 or loss coverage of ≥3 (we can make separate tables for gains and for losses)
    # List the interval coverage
    gain_intervals <- dftdLowCov::get_depth(data[Type == "gain"], chrlengths)[depth >= required_depth]
    gain_intervals[, Type := "gain"]
    loss_intervals <- dftdLowCov::get_depth(data[Type == "loss"], chrlengths)[depth >= required_depth]
    loss_intervals[, Type := "loss"]
    intervals <- rbindlist(list(gain_intervals, loss_intervals))
    intervals[, intervalID := .I]
    rm(gain_intervals)
    rm(loss_intervals)
    
    # List the CNVs covering the interval
    data[, end := end - 1]
    keycols <- c("chr", "Type", "start", "end")
    setkeyv(data, keycols)
    setkeyv(intervals, keycols)
    ol <- foverlaps(data, intervals, nomatch=0L)
    
    # List the samples and clades associated with the CNV IDs (you could just add the data in the CNV table samples_CNV_group column, making it clear which CNV ID each of the samples relate to)
    intervals <- ol[, .(chr, Type, start, end, depth, width, intervalID, CNV_ID, CNV_ID_ext, Sample_group, `Linkage_group_gains`, `Linkage_group_losses`)]
    rm(ol)
    
    # List the genes annotated fully or partially in the interval, defaulting to ENSSHAG ID if HGNC symbol not available
    gene_list <- gene_translation[gene_list, on = c(ID="DevilID")]
    gene_list[HGNC %in% c("NA", ""), HGNC := NA_character_]
    gene_list[`External Name` %in% c("NA", ""), `External Name` := NA_character_]
    gene_list[DevilGeneName %in% c("NA", ""), DevilGeneName := NA_character_]
    
    gene_list[!is.na(HGNC), gene_name := HGNC]
    gene_list[is.na(HGNC) & !is.na(DevilGeneName), gene_name := DevilGeneName]
    gene_list[is.na(HGNC) & is.na(DevilGeneName) & !is.na(`External Name`), gene_name := `External Name`]
    gene_list[is.na(HGNC) & is.na(DevilGeneName) & is.na(`External Name`), gene_name := ID]
    setcolorder(gene_list, c("chr", "start", "end", "gene_name"))
    
    # If any genes on the list are present in the COSMIC Cancer Gene Census, then add additional columns with the COSMIC Tier (1 or 2), and the data from the following columns ‘Molecular Genetics’, ‘Role in Cancer’, ‘Mutation Types’
    cosmic[, EnsemblID := scan_col("(ENSG\\d+)", Synonyms, c('c'), c("EnsemblID"))]
    gene_list[human_to_devil, HumanID := HumanID, on = c(ID = "DevilID")]
    gene_list[cosmic[!is.na(EnsemblID)], c("COSMIC_GeneSymbol", "COSMIC_Tier", "COSMIC_MolecularGenetics", "COSMIC_RoleInCancer", "COSMIC_MutationTypes") := .(`Gene Symbol`, Tier, `Molecular Genetics`, `Role in Cancer`, `Mutation Types`),
              on = c("HumanID" = "EnsemblID")]
    
    setkey(gene_list, chr, start, end)
    setkey(intervals, chr, start, end)
    ol <- foverlaps(gene_list, intervals, nomatch=0L)
    
    # Break extra info into separate tables
    cnv_info <- unique(intervals[, .(intervalID, Type, CNV_ID, CNV_ID_ext, Sample_group, `Linkage_group_gains`, `Linkage_group_losses`)])[order(intervalID)]
    gene_info <- unique(ol[, .(intervalID, Type, gene_name, chr, start = i.start, end = i.end, COSMIC_GeneSymbol, COSMIC_Tier, COSMIC_MolecularGenetics, COSMIC_RoleInCancer, COSMIC_MutationTypes)])[order(intervalID)]
    
    # Make summary string of cosmic info
    gene_info[!is.na(COSMIC_GeneSymbol), COSMIC_Summary := sprintf("%s (%s|%s|%s|%s)",
                                                                   COSMIC_GeneSymbol,
                                                                   COSMIC_Tier,
                                                                   COSMIC_MolecularGenetics,
                                                                   COSMIC_RoleInCancer,
                                                                   COSMIC_MutationTypes)]
    
    interval_gene_lists <- ol[, genes := paste(unique(gene_name), collapse = ","), by = intervalID]
    interval_cnv_lists <- intervals[, cnvs := paste(unique(CNV_ID_ext), collapse = ","), by = intervalID]
    contains_cosmic <- unique(ol[!is.na(COSMIC_GeneSymbol), .(intervalID)])
    
    intervals[interval_gene_lists, genes := genes, on = "intervalID"]
    intervals[interval_cnv_lists, cnvs := cnvs, on = "intervalID"]
    intervals <- unique(intervals[, .(chr, Type, start, end, depth, width, intervalID, cnvs, genes)])[order(intervalID)]
    intervals[, cosmic := FALSE]
    intervals[contains_cosmic, cosmic := TRUE, on = "intervalID"]
    intervals[gene_info[!is.na(COSMIC_Summary), .(pasted_summary = paste(COSMIC_Summary, collapse = ", ")), by = intervalID],
              COSMIC_Summary := pasted_summary,
              on = "intervalID"]
    
    wb=openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "intervals")
    openxlsx::addWorksheet(wb, "cnvs")
    openxlsx::addWorksheet(wb, "genes")
    openxlsx::writeData(wb, "intervals", intervals)
    openxlsx::writeData(wb, "cnvs", cnv_info)
    openxlsx::writeData(wb, "genes", gene_info)
    openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
}

dft1_cnvs <- isDFT1(cnv_table)
cellline_cnvs <- isCellLine(cnv_table)
dft2_cnvs <- isDFT2(cnv_table)
nondftd_cnvs <- isNonDFTD(cnv_table)

# DFT1-only biopsy data
data <- copy(cnv_table[CNV_ID_ext %in% setdiff(dft1_cnvs, cellline_cnvs)])
filename <- file.path(outdir, "TableS6A_annotated_intervals_DFT1_biopsies.xlsx")
make_annotation_table(data, filename)

# DFT1-only cell-line data
data <- copy(cnv_table[CNV_ID_ext %in% intersect(dft1_cnvs, cellline_cnvs)])
filename <- file.path(outdir, "TableS6B_annotated_intervals_DFT1_celllines.xlsx")
make_annotation_table(data, filename)

# DFT2-only data
data <- copy(cnv_table[CNV_ID_ext %in% dft2_cnvs])
filename <- file.path(outdir, "TableS6C_annotated_intervals_DFT2.xlsx")
make_annotation_table(data, filename, required_depth = 2)

# Non-DFTD data
data <- copy(cnv_table[CNV_ID_ext %in% nondftd_cnvs])
filename <- file.path(outdir, "TableS6D_annotated_intervals_NonDFTD.xlsx")
make_annotation_table(data, filename, required_depth = 2)
