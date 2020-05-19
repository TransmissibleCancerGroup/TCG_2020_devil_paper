##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 01 of the pipeline for analysing recurrent alleles.
# Loads the data and does some preprocessing.


library(ape)
library(data.table)
library(openxlsx)
library(dftdLowCov)
library(logging); basicConfig()
logger <- getLogger("PIPELINE_01")

########################################################
# Define the file paths used to load data              #
# There should be no other hard coded paths than these #
########################################################
DATA.DIR <- "/Users/kg8/Documents/projects/DFTD_low_coverage/data"
TUMOUR_DATA_FILE <- "/Users/kg8/Documents/projects/DFTD_low_coverage/data/tumour_data/haplotype_database_LATEST.xlsx"
# TUMOUR_DATA_FILE <- "/Users/kg8/Documents/projects/DFTD_low_coverage/data/tumour_data/2019-12-09_1232_haplotype_database.xlsx"
MIPS_TREE_DATA_FILE <- "/Users/kg8/Documents/projects/DFTD_low_coverage/dftd_combined_data_trees/trees/mitochondrial_partition/part.nex.best_scheme.nex.treefile" # file.path(DATA.DIR, "treeMP_nimputed_duplicated_corrected_amb_50_threshold_50_maxlikelihood_bootstrapped_4_12_18.nwk")
CNV_TABLE_FILE <- file.path(DATA.DIR, "cnv_table_finalised/CNV_summary_table_LATEST.xlsx")
MIPS_VAF_TABLE_FILE <- file.path(DATA.DIR, "mips_data/latest_from_gdrive/VAF_table_for_Kevin_easier_format_for_matching.csv")
MIPS_PRES_ABS_TABLE_FILE <- file.path(DATA.DIR, "mips_data/latest_from_gdrive/presence_absence_mips_table.csv")

split_counts <- function(dt) {
    dt[, vaf_counts := sub("'", "", vaf_counts)]
    counts <- unique(dt[, .(vaf_counts)])
    counts[, c("alt", "total") := scan_col("(.*)/(.*)", vaf_counts, c('i', 'i'))]
    counts[, vaf := alt / total]
    dt[counts, c("alt", "total", "vaf") := .(alt, total, vaf), on = .(vaf_counts)]
    dt[1:.N]
}

load_new_vaf_table <- function(filename) {
    raw = fread(filename)
    C <- ncol(raw)
    samplenames   <- as.character(unlist(raw[V1=="Sample.TCG_ID", 2:C]))
    mode_vaf      <- as.numeric(unlist(raw[V1=="Mode_VAF", 2:C]))
    median_vaf    <- as.numeric(unlist(raw[V1=="Median_VAF", 2:C]))
    mt_haplogroup <- as.character(unlist(raw[V1=="Mt Haplogroup", 2:C]))
    clade         <- as.character(unlist(raw[V1=="Clade", 2:C]))
    metadata <- data.table(sample = samplenames,
                           mode_vaf = mode_vaf,
                           median_vaf = median_vaf,
                           mt_haplogroup = mt_haplogroup,
                           clade = clade)
    
    raw <- raw[startsWith(V1, "Chr")]
    colnames = c("mip_name", samplenames)
    setnames(raw, colnames)
    vaftable <- melt.data.table(raw, id.vars = "mip_name", value.name = "vaf_counts", variable.name = "sample", variable.factor = FALSE, value.factor = FALSE)
    vaftable <- split_counts(vaftable)
    setkey(vaftable, mip_name, sample)
    setcolorder(vaftable)
    return(list(vaftable = vaftable[1:.N], metadata = metadata[1:.N]))
}

load_new_vaf_table_for_me <- function(filename) {
    raw <- fread(filename)
    setnames(raw, old = c("Sample.TCG_ID", "VAF"), new = c("sample", "vaf_counts"))
    vaftable <- split_counts(raw)
    setnames(vaftable, old=c("chrom_chr", "chrom_pos"), new=c("seqnames", "mip_position"))
    vaftable[, CHROM := NULL]
    vaftable[, POS := NULL]
    setkey(vaftable, mip_name, sample)
    setcolorder(vaftable)
    vaftable[1:.N]
}

load_new_presence_absence_table <- function(filename) {
    raw <- fread(filename)
    raw <- raw[startsWith(mip_name, "Chr")]
    suppressWarnings(
        presabs <- melt.data.table(raw, id.vars = c("mip_name", "REF", "ALT"),
                                   value.name = "indicator",
                                   variable.name = "sample",
                                   variable.factor = FALSE,
                                   value.factor = FALSE)
    )
    presabs[, mip_allele := ifelse(indicator == "0", "REF", ifelse(indicator == "1", "ALT", "N"))]
    #presabs[, indicator := NULL]
    setkey(presabs, mip_name, sample)
    setcolorder(presabs)
    presabs[1:.N]
}

load_tumour_data <- function(excel_filename) {
    require(data.table)
    # Peek
    tumour_data <- as.data.table(read.xlsx(TUMOUR_DATA_FILE, detectDates = TRUE))
    tumour_data
}

load_cnv_data <- function(filepath) {
    if (endsWith(filepath, "xlsx")) {
        cnv_data <- as.data.table(read.xlsx(filepath))
    } else {
        cnv_data <- fread(cnv_path)
    }
    setnames(cnv_data, make.unique(colnames(cnv_data)))
    #setnames(cnv_data, old = c("204T", "204T.1"), new = c("Clade_204T", "204T"))
    if (!("seqnames" %in% colnames(cnv_data))) {
        cnv_data[, seqnames := chr]
    }
    setnames(cnv_data,
             old = c("New.CNV_ID_ext", "Type"),
             new = c("cnv_id", "State"))
    samplecols = grep("^\\d+T[1-9.a]*", colnames(cnv_data), value = T)
    cnv_data[, (samplecols) := lapply(.SD, as.integer), .SDcols = samplecols]
    cnv_data <- melt(cnv_data, id.vars = c("seqnames", "start", "end", "State", "cnv_id"),
         measure.vars = patterns("^\\d+T[1-9.ab]*"), value.name = "cnv_present_in_sample", variable.name = "sample",
         value.factor = FALSE, variable.factor = FALSE)
    setkey(cnv_data, cnv_id, sample)
    setcolorder(cnv_data)
    cnv_data[cnv_id != ""]
}

# Getting the copy number in each aberrant sample, for each unique CNV...
# Data: cnv data table
# Relevant columns: New.CNV_ID, Sample_CNV_Group
get_cnv_copynumber <- function(cnv_path) {
    parse_element <- function(s) {
        p <- strsplit(s, ",")[[1]]
        pattern <- "([0-9Tab.]*) \\(.*\\) [>=]*([+-]?[0-9]*[.]?[0-9]+)"
        as.data.table(scan_col(pattern, p, c('c', 'n'), names = c("sample", "copy_number")))
    }
    
    if (endsWith(cnv_path, "xlsx")) {
        cnvs <- as.data.table(read.xlsx(cnv_path))
    } else {
        cnvs <- fread(cnv_path)
    }
    elements <- lapply(cnvs$Sample_CNV_Group, parse_element)
    for (i in 1:length(elements)) {
        elements[[i]][, cnv_id := cnvs[i, New.CNV_ID_ext]]
    }
    cn <- rbindlist(elements)
    setkey(cn, cnv_id, sample)
    setcolorder(cn)
    cn[1:.N]
}

load_tree_data <- function(filename) {
    tree <- read.tree(filename)
    tree$tip.label <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
    ladderize(tree)
}



newvaftableforme <- load_new_vaf_table_for_me(MIPS_VAF_TABLE_FILE)
newpresabstable <- load_new_presence_absence_table(MIPS_PRES_ABS_TABLE_FILE)

newpresabstable <- newpresabstable[sample != "reference"]
mips_data <- newvaftableforme[newpresabstable, on = .(mip_name, sample, REF, ALT)]
setkey(mips_data, mip_name, sample)

rm(newvaftableforme, newpresabstable)

haplotype_data <- load_tumour_data(TUMOUR_DATA_FILE)
suppressWarnings(haplotype_data[, purity := as.numeric(Tumor.Purity.MIP)])
suppressWarnings(haplotype_data[is.na(purity), purity := as.numeric(`Tumour.purity..mean.CNVs..MIPs.`)])
mips_data[haplotype_data, c("purity", "clade") := .(purity, Group), on = c(sample = "Sample.TCG_ID")]
logger$info("Loaded MIP VAF data as 'mips_data'")

cnv_data <- load_cnv_data(CNV_TABLE_FILE)
logger$info("Loaded CNV data as 'cnv_data'")


cn_data <- get_cnv_copynumber(CNV_TABLE_FILE)
cn_data <- cn_data[cnv_id != ""]

# 11/09/19 Liz's by-hand correction to copy number
cn_data[cnv_id == "10R10BM14" & sample == "619T1", copy_number := 2.5]
cn_data[cnv_id == "10R10BM15" & sample == "620T1", copy_number := 2.5]
cn_data[cnv_id == "10R10BM4" & sample == "218T",   copy_number := 2.5]
cn_data[cnv_id == "10R10BM7" & sample == "221T",   copy_number := 2.5]
cn_data[cnv_id == "288" & sample == "993T1.1",     copy_number := 1.5]
cn_data[cnv_id == "3R1BM23" & sample == "135T6a",  copy_number := 2.5]
# 24/09/19
cn_data[cnv_id %in% c("3R1", "9", "10"), copy_number := 3.0] # In all samples, these CNVs have a copy number of 3

logger$info("Loaded copy number data as 'cn_data'")

tree <- load_tree_data(MIPS_TREE_DATA_FILE)
logger$info("Loaded tree data as 'tree'")