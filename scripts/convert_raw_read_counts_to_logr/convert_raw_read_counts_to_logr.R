##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: September 2020
#
# WHAT THIS FILE DOES:
# Converts Table S4 from binned read counts to binned logR estimates for host and tumour samples
#
# Usage:
# Rscript [this_file.R] <INPUT_FILE> <OUTPUT_FILE_HOSTS> <OUTPUT_FILE_TUMOURS>
#
# Inputs:
# A CSV file containing the raw read counts for host and tumour samples.
# This is provided as supplementary table S4. This script will automatically
# unpack the file if it is compressed, as long as the filename ends '.zip',
# '.gz' or '.bz2'.
#
# Outputs:
# Two CSV files, one containing the host samples' logR values, the other
# the tumour samples' logR values.

library(argparse)
library(data.table)

###########################
#       Young Mi's        #
# Normalization functions #
###########################
# Internal Function
log2R_samples<-function(x){
    b<-log2(x/median(x, na.rm = TRUE))
    b[which(is.infinite(b))]<-0
    return(b)
}

# Note: this is not the GATK version:
# DataFrames Format: Chr, Start, End, Sample name 1, Sample name 2, ...)
PCA_normalization<-function(host, tumor){

    # Step 1: LogR the Raw Read Counts
    cat("  - Taking logs of conversion data\n")
    host_log<-apply(host[,4:ncol(host)], 2, log2R_samples)
    cat("  - Taking logs of normalisation data\n")
    tumor_log<-apply(tumor[,c(4:ncol(tumor))], 2, log2R_samples)


    # Step 2: Calculate PCA
    cat("  - Doing PCA\n")
    res<-prcomp(t(host_log), center = TRUE, scale = FALSE)

    # Selected the First 4 PC:
    test1<-predict(res, newdata=t(tumor_log))[,1:4] %*% t(res$rotation[,1:4])

    # Created your noise signal for each tumor and centered it on your PON
    cat("  - Identifying noise signal\n")
    trunc <- scale(test1, center = -1 * res$center, scale=FALSE)

    # Remove the noise signal from your tumor
    cat("  - Removing noise from converted logRs\n")
    tumors_pca<-sapply(c(1:ncol(tumor_log)), function(x){
        df<-tumor_log[,x]-trunc[x,]
    })
    colnames(tumors_pca)<-colnames(tumor_log)

    return(tumors_pca)}


#############################
# Handle command line input #
#############################
args <- commandArgs(trailingOnly = TRUE)

parser <- argparse::ArgumentParser(description = "Converts Table S4's raw read counts to logR")
parser$add_argument("infile", help = "Read counts file, i.e. S4_Data.zip.")
parser$add_argument("output_filename_host", help = "File to write host sample logRs in .csv format")
parser$add_argument("output_filename_tumour", help = "File to write tumour sample logRs in .csv format")

args <- parser$parse_args()

filename <- args$infile
output_filename_host <- args$output_filename_host
output_filename_tumour <- args$output_filename_tumour

if (!file.exists(filename)) {
    stop(sprintf("Input file %s does not exist.", filename))
}

if (!dir.exists(dirname(output_filename_host))) {
    stop(sprintf("Output file %s is in a directory that does not exist.", output_filename_host))
}

if (!dir.exists(dirname(output_filename_tumour))) {
    stop(sprintf("Output file %s is in a directory that does not exist.", output_filename_tumour))
}

if (endsWith(filename, ".zip")) {
    suppressWarnings(read_depths <- fread(cmd = paste0("unzip -p ", filename)))
} else {
    read_depths <- fread(filename)
}

hosts <- read_depths[, .SD, .SDcols = c(1,2,3, grep("H", colnames(read_depths)))]
tumours <- read_depths[, .SD, .SDcols = c(1,2,3, grep("T", colnames(read_depths)))]

setDF(hosts)
setDF(tumours)

cat("Converting host samples to logR\n")
hosts_logr <- PCA_normalization(hosts, hosts)
cat("Converting tumour samples to logR\n")
tumours_logr <- PCA_normalization(hosts, tumours)

cat(sprintf("Writing Host logRs to %s\n", output_filename_host))
hosts_logr <- cbind(as.data.table(hosts[, 1:3]), as.data.table(hosts_logr))
fwrite(hosts_logr, output_filename_host)

cat(sprintf("Writing Tumour logRs to %s\n", output_filename_tumour))
tumours_logr <- cbind(as.data.table(tumours[, 1:3]), as.data.table(tumours_logr))
fwrite(tumours_logr, output_filename_tumour)
