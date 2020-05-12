##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# This file prepares the three sequence alignments used for making the phylogenetic tree.
#
# Details:
# We want to make three sequence alignments, from the following data sources:
#  - mitochondrial DNA variants
#  - CNV presence-absence data
#  - Molecular inversion probe sequencing (MIPS) DNA variants
# Each alignment will contain only variant sites (the Lewis ascertainment bias correction model
# is used in the phylogenetic analysis to account for this).
#
# This script extracts the Mitochondrial and CNV data from the supplementary tables and converts
# them into fasta alignment files. The MIPS data was provided to me already in fasta format, and
# is included as external data in the dftdLowCov package.
#
# After generating the alignments, this script does a filtering step to exclude certain low
# information sequences.
#
# Finally, the filtered alignments are written to an output folder, along with a constraint tree
# for the phylogenetic analysis.

library(dftdLowCov)
library(data.table)
library(readxl)
library(seqinr)
library(logging); basicConfig()

# Output files
mito_fasta_filename <- "mitochondrial_alignment.fas"
cnv_fasta_file <- "cnv_alignment.fas"
output_folder <- "filtered_alignments"  # Filtered alignments written to this folder

############################
# Begin Function Definitions
############################

# Extract mitochondrial variant info from the Mitochondrial.variants column
parse_variants <- function(s) {
    if (!grepl(">", s)) {
        return (NA_character_)
    }
    if (grepl(",", s) & !grepl(";", s)) {
        comma_split <- strsplit(s, ",")
        return (comma_split[[1]])
    } else {
        if (grepl(";", s) & !grepl(",", s)) {
            semicolon_split <- strsplit(s, ";")
            return (semicolon_split[[1]])
        } else {
            return (s)
        }
    }
    stop(sprintf("Shouldn't happen (caused by %s)", s))
}

# Reads a single row of the table
read_row <- function(i) {
    id <- hapdata[i, Sample_ID]
    variants <- parse_variants(hapdata[i, Mitochondrial.variants])
    return (data.table(Sample_ID = id, variant_string = variants))
}

write_fasta_from_matrix <- function(m, filename) {
    sequences <- as.matrix(apply(m, 1, paste0, collapse=""))
    write(paste0(">", rownames(sequences),"\n",sequences,"\n"), file = filename)
}

seq_to_matrix <- function(table) {
    matrix(unlist(sapply(apply(table[, .(sequence)], 1, strsplit, ""), `[`, 1)), byrow = TRUE, ncol = nchar(table[1, sequence]))
}

matrix_to_seq <- function(mat) {
    apply(mat, 1, paste0, collapse="")
}

output_all_fasta_files <- function(seqtable, outpath = ".") {
    # MIPs
    m <- seq_to_matrix(seqtable[origin == "MIP"])
    mip_variant <- apply(m, 2, is_variant_site, alphabet = "DNA")
    seqtable[origin == "MIP", sequence := matrix_to_seq(m[, mip_variant])]
    
    # Mitochondria
    m <- seq_to_matrix(seqtable[origin == "MITO"])
    mito_variant <- apply(m, 2, is_variant_site, alphabet = "DNA")
    seqtable[origin == "MITO", sequence := matrix_to_seq(m[, mito_variant])]
    
    # CNVs
    m <- seq_to_matrix(seqtable[origin == "CNV"])
    cnv_variant <- apply(m, 2, is_variant_site, alphabet = "BINARY")
    seqtable[origin == "CNV", sequence := matrix_to_seq(m[, cnv_variant])]
    
    # Write filtered FASTA files
    timenow <- format(Sys.time(), "%Y-%m-%d")
    if (!dir.exists(outpath)) dir.create(outpath)
    mip_filename <- paste0("mips_filtered_dna_alignment_dft1_", timenow, ".fa")
    write_fasta(file.path(outpath, mip_filename), seqtable[origin == "MIP"])
    loginfo("Wrote MIP alignment to %s", file.path(outpath, mip_filename))
    
    mito_filename <- paste0("mito_filtered_dna_alignment_dft1_", timenow, ".fa")
    write_fasta(file.path(outpath, mito_filename), seqtable[origin == "MITO"])
    loginfo("Wrote mitochondrial alignment to %s", file.path(outpath, mito_filename))
    
    cnv_filename <- paste0("cnv_filtered_alignment_dft1_", timenow, ".fa")
    write_fasta(file.path(outpath, cnv_filename), seqtable[origin == "CNV"])
    loginfo("Wrote CNV alignment to %s", file.path(outpath, cnv_filename))
}

write_constraint_tree <- function(seqtable, outfile) {
    all_names <- unique(seqtable[, label])
    
    clade_a1 <- grep("Clade_A1", all_names, value=TRUE)
    clade_a2 <- grep("Clade_A2", all_names, value=TRUE)
    clade_b <- grep("Clade_B", all_names, value=TRUE)
    clade_b2 <- grep("Clade_B-2_", all_names, value=TRUE)
    clade_b_excl_b2 <- setdiff(clade_b, clade_b2)
    clade_c <- grep("Clade_C", all_names, value=TRUE)
    clade_d <- grep("Clade_D", all_names, value=TRUE)
    clade_e <- grep("Clade_E", all_names, value=TRUE)
    reference <- grep("reference", all_names, value=TRUE)
    
    constr <- sprintf("(%s,((%s),(%s)),(%s,(%s)),(%s),(%s),%s);",
                      reference,
                      paste(clade_a1, collapse = ","),
                      paste(clade_a2, collapse = ","),
                      paste(clade_b_excl_b2, collapse = ","),
                      paste(clade_b2, collapse = ","),
                      paste(clade_c, collapse = ","),
                      paste(clade_d, collapse = ","),
                      paste(clade_e, collapse = ","))
    
    cat(constr, file = outfile)
    loginfo("Wrote constraint tree to %s", outfile)
}
# End function definitions

#########################################################################
# Step 1: Produce the mitochondrial alignment using data in the haplotype 
# table (Table S1A)
#########################################################################

haplotype_db <- dftdLowCov::load_sample_table()
hapdata <- haplotype_db[Cancer.type == "DFT1" & Mitochondrial.variants != "Unknown", .(Sample_ID, Mitochondrial.variants)]

variants <- rbindlist(lapply(1:nrow(hapdata), read_row))
variants[, c("position", "ref", "alt") := scan_col("^\\s?(\\d+)(\\w)>(\\w)", variant_string, c('i','c','c'))]

Sample_ID = variants$Sample_ID
variant_string = variants[!is.na(variant_string), variant_string]
setkey(variants, Sample_ID, variant_string)
all_variants <- variants[CJ(Sample_ID, variant_string, unique = TRUE)][!is.na(variant_string)]

all_variants[, variant_is_present := !is.na(position)]
all_variants[, c("position", "ref", "alt") := scan_col("^\\s?(\\d+)(\\w)>(\\w)", variant_string, c('i','c','c'))]
all_variants[, value := ifelse(variant_is_present, alt, ref)]
setorder(all_variants, position)
sequences <- all_variants[, .(sequence = paste(value, collapse = "")), by = Sample_ID]
refseq <- all_variants[Sample_ID == "9T2", paste(ref, collapse = "")]
sequences <- rbindlist(list(sequences, data.table(Sample_ID = "reference", sequence = refseq)))

# Check for invariant sites
m <- matrix(unlist(sapply(apply(sequences[, .(sequence)], 1, strsplit, ""), `[`, 1)), byrow = TRUE, ncol = nchar(sequences[1, sequence]))
check <- apply(m, 2, is_variant_site, alphabet = "DNA")
stopifnot(all(check))

# Annotate labels
seqnames <-
    haplotype_db[sequences,
                 on = "Sample_ID",
                 paste(Sample_ID, Clade.group, gsub(" ", "_", Location), sep = "_")]
sequences[, label := seqnames]
sequences[, label := sub("\\(captive\\)", "captive", label)]

# Write fasta
dftdLowCov::write_fasta(mito_fasta_filename,
                        sequences)
loginfo("Wrote mitochondrial fasta file to %s", mito_fasta_filename)


#########################################################################
# Step 2: Produce the CNV alignment using data in the copy number variant
# table (Table S3)
#########################################################################


cnv_table <- dftdLowCov::load_cnv_table()

# Collect the DFT1 samples
dft1_samples <- haplotype_db[Cancer.type == "DFT1", Sample_ID]

# Restrict to the set for which we have data in the table
dft1_samples <- intersect(dft1_samples, colnames(cnv_table))

# Remove the 'Do not use to build tree' variants
cnv_table <- cnv_table[Cancer_type == "DFT1" & is.na(Do_not_use_to_build_tree)]
m <- as.matrix(cnv_table[, .SD, .SDcols = dft1_samples])
row.names(m) <- cnv_table[, CNV_ID_ext]

# Filter out rows (CNVs) which are phylogenetically uninformative
variant_sites <- apply(m, 1, is_variant_site, alphabet = "BINARY")
m <- m[variant_sites, ]

# Convert to sequence strings
sequences <- apply(m, 2, paste0, collapse = "")
sequences <- data.table(Sample_ID = names(sequences), sequence = sequences)

# Annotate with clade and location data
haplotype_db[, annotation := paste(Sample_ID, Clade.group, Location, sep = "_")]
sequences[haplotype_db, label := annotation, on = "Sample_ID"]
sequences[, label := gsub(" ", "_", label)]
sequences[, label := sub("\\(captive\\)", "captive", label)]

# Write in FASTA format
dftdLowCov::write_fasta(cnv_fasta_file, sequences)
loginfo("Wrote CNV fasta file to %s", cnv_fasta_file)

#####################################
# Step 3: Filter the three alignments
#####################################

# Read in the FASTA alignments
mip_aln <- seqinr::read.alignment(system.file("extdata", "mips_dna_alignment_dft1.fa", package = "dftdLowCov"),
                                  "fasta")
cnv_aln <- seqinr::read.alignment(cnv_fasta_file, "fasta")
mito_aln <- seqinr::read.alignment(mito_fasta_filename, "fasta")

# Do some filtering and checking of sample names between files
# !!Remove 356T1b from MIP alignment!!
index <- !grepl("356T1b", mip_aln$nam) # Bool vector TRUE for everything that is not 356T1b
mip_aln <- seqinr::as.alignment(
    nb = sum(index),
    nam = mip_aln$nam[index],
    seq = mip_aln$seq[index],
    com = mip_aln$com
)
rm(index)

# Sanitise MIP alignment names
mip_aln$nam <- sub("\\(captive\\)", "captive", mip_aln$nam)

# Do a check to make sure the long form label for each sample is consistent
# between the three alignments (got burned by this before!)
all_names <- Reduce(union, list(mip_aln$nam, mito_aln$nam, cnv_aln$nam))
checker <- data.table(long_name = all_names)
checker[, short_name := tstrsplit(long_name, "_")[1]]
# Should only ever be 1 unique long name per unique short name
stopifnot(checker[, .(check=uniqueN(long_name)), by = short_name][, all(check==1)])
stopifnot(checker[, .(check=uniqueN(short_name)), by = long_name][, all(check==1)])
rm(all_names, checker)

# Put all information into data table for processing
mip_dt <- data.table(label = mip_aln$nam, sequence = toupper(mip_aln$seq), origin = "MIP")
cnv_dt <- data.table(label = cnv_aln$nam, sequence = toupper(cnv_aln$seq), origin = "CNV")
mito_dt <- data.table(label = mito_aln$nam, sequence = toupper(mito_aln$seq), origin = "MITO")
seq_dt <- rbindlist(list(mip_dt, cnv_dt, mito_dt))
seq_dt[, sample_ID := sapply(strsplit(label, "_"), `[`, 1)]

# Samples to exclude
EXCLUDE <- c("8T2", "9T", "120T", "166T4", "390T1", "460T1", "478T2", "1070T1", "1072T1")
DUPLICATES <- haplotype_db[Cancer.type == "DFT1" & !is.na(Duplicate.or.time.course.biopsy)][Duplicate.or.time.course.biopsy == "Duplicate DNA extraction from same tumour", Sample_ID]

# Flag samples that are in our exclusion list, UNINFORMATIVE_SAMPLES + DUPLICATES
seq_dt[, is_duplicate := sample_ID %in% DUPLICATES]
seq_dt[, is_exclude := sample_ID %in% EXCLUDE]

# Write filtered FASTA files and constraint trees
if (!dir.exists(output_folder)) dir.create(output_folder)
output_all_fasta_files(seq_dt[!(is_duplicate) & !(is_exclude)], output_folder)
write_constraint_tree(seq_dt[!(is_duplicate) & !(is_exclude)], file.path(output_folder, "constraint.tree"))
