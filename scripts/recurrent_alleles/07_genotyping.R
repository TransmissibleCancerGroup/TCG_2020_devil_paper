##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 07 of the pipeline for analysing recurrent alleles.
# Allele-specific genotyping of copy number changes.


logger <- getLogger("PIPELINE_07")

#' Calculate what the VAF would be for a given genotype.
#' The input row contains information on the copy number and purity of a sample
#' at its CNV locus. This function generates a sequence of genotypes in 0.5 copy number intervals
#' and computes what the VAF is expected to be, at the given purity and copy number.
idealisedVAFfn <- function(row) {
    
    idealised_vaf <- function(ta, tb, ha, hb, purity) {
        (hb * (1 - purity) + tb * purity) / ((ha+hb) * (1-purity) + (ta+tb) * purity)
    }
    
    # Read local variables from row
    copy_number_ <- as.numeric(row[["copy_number"]])
    cnv_id_ <- row[["cnv_id"]]
    mip_name_ <- row[["mip_name"]]
    sample_ <- row[["sample"]]
    State_ <- row[["State"]]
    purity_ <- as.numeric(row[["purity"]])
    
    result_on_failure <- data.table(cnv_id = cnv_id_,
                                    mip_name = mip_name_,
                                    sample = sample_,
                                    State = State_,
                                    genotype = NA,
                                    idealised_VAF = NA)
    if (is.na(copy_number_)) {
        logwarn("Failure due to missing copynumber")
        return(result_on_failure)
    }
    
    subclonal = (copy_number_ %% 1 > 0)
    if (subclonal & copy_number_ %% 1 != 0.5) {
        logwarn("Failure due to subclonality != 0.5")
        return (result_on_failure)
    }
    
    ta <- seq(0, copy_number_, 0.5)
    
    vaf <- sapply(ta, function(tref) {
        idealised_vaf(tref, copy_number_ - tref, 2, 0, purity_)
    })
    
    data.table(cnv_id = cnv_id_,
               mip_name = mip_name_,
               sample = sample_,
               State = State_,
               best_fit_genotype = paste(ta, copy_number_ - ta, sep = "|"),
               idealised_vaf = vaf)
}

#' Infer the amount of B-allele present in the tumour, using the
#' estimated VAF, purity and copy number values.
tb <- function(vaf, pur, totcn) {
    vaf * (2 * (1 - pur) + totcn * pur) / pur
}

logger$info("Computing best fit genotype")
# Of the series of generated genotype VAFs, choose the one that is closest to the observed VAF
data <- mips_data[!is.na(copy_number) & !is.na(vaf)]#[(CVG.PASS) & (INFO.PASS) & cnv_present_in_sample == 1]
tmp <- rbindlist(apply(data, 1, idealisedVAFfn))
merged <- merge(data, tmp, by = c("cnv_id", "mip_name", "sample", "State"), all = TRUE)
merged[, score := abs(vaf - idealised_vaf)]
merged[, best := score == min(score), by = .(cnv_id, mip_name, sample, State)]
genotyped <- unique(merged[(best == TRUE) | is.na(best)])

logger$info("Computing estimated genotype")
# Use tb function to infer the genotype directly, with the simplification that if the
# amount of B-allele is greater than the total copynumber, then there is 0 A-allele, not a negative
# amount.
genotyped[, est_alt := tb(vaf, purity, copy_number)]
genotyped[, est_ref := pmax(copy_number - est_alt, 0)]
genotyped[, estimated_genotype := paste(round(est_ref, 2), round(est_alt, 2), sep = "|")]
genotyped[, diff := est_ref - est_alt]


logger$info("Annotating 'mips_data' with genotype")
mips_data[genotyped,
          c("best_fit_genotype", "estimated_genotype", "diff", "est_ref", "est_alt") := .(best_fit_genotype, estimated_genotype, diff, est_ref, est_alt),
          on = intersect(colnames(mips_data), colnames(genotyped))]

mips_data[, FILTERS.PASS := (CVG.PASS) & (INFO.PASS) & cnv_present_in_sample == 1]



logger$info("Attempting to automatically assign implicated allele")
mips_data[cnv_present_in_sample == 1 & diff > 0.5, implicated_allele := ifelse(State == "gain", "REF", "ALT")]
mips_data[cnv_present_in_sample == 1 & diff > 0.4 & diff <= 0.5, implicated_allele := ifelse(State == "gain", "BORDERLINE_REF", "BORDERLINE_ALT")]
mips_data[cnv_present_in_sample == 1 & diff < -0.5, implicated_allele := ifelse(State == "gain", "ALT", "REF")]
mips_data[cnv_present_in_sample == 1 & diff < -0.4 & diff >= -0.5, implicated_allele := ifelse(State == "gain", "BORDERLINE_ALT", "BORDERLINE_REF")]
mips_data[cnv_present_in_sample == 1 & diff >= -0.4 & diff <= 0.4, implicated_allele := "AMBIGUOUS"]
mips_data[cnv_present_in_sample == 1 & State == "gain" & est_ref < 0.25, implicated_allele := "ANCESTOR_HOMZYG_ALT"]
mips_data[cnv_present_in_sample == 1 & State == "gain" & (est_alt < 0.1 | INFO.PASS == FALSE), implicated_allele := "ANCESTOR_HOMZYG_REF"]
mips_data[cnv_present_in_sample == 1 & State == "loss" & cnv_id %like% "BM", implicated_allele := "BACKMUT"]

logger$info("Determining pre-back mutated state")
# If we select the back mutating losses, we can intersect their positions with known gain cnvs in the same samples
# and try to infer the genotype before the back mutation occurred
back_mutations <- mips_data[cnv_present_in_sample == 1 & State == "loss" & cnv_id %like% "BM"]
setkey(back_mutations, seqnames, mip_name, sample, start, end)
gains <- mips_data[cnv_present_in_sample == 1 & State == "gain"]
setkeyv(gains, key(back_mutations))
ol <- foverlaps(gains, back_mutations, nomatch = 0L)
ol[copy_number != i.copy_number]
