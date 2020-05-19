##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 06 of the pipeline for analysing recurrent alleles.
# Coverage filter. Remove MIPS SNVs that are not associated with at least 2 CNVs.


logger <- getLogger("PIPELINE_06")

logger$info("Filtering by CNV coverage - each MIP must be associated with at least 2 informative CNVs")
logger$info("either both in the same direction (2 gains / 2 losses), or involving a back mutation.")

# Currently I have a table that contains every combination of cnv_id, mip_name and sample, among the
# full set of cnvs, mips, and samples for which there is cnv or mip information
# I need to filter this down to 
# * IDs (mip,cnv,sample,State) in which the cnv is present in the sample
# * MIPs for which there are at least 2 combinations of mip+cnv among either state=loss or state=gain
filter_coverage_independent_gains_losses <- function(mips_data_table) {
    logger$info("Filtering for mips that are associated with at least two independent gains or two independent losses")
    tmp <- mips_data_table[cnv_present_in_sample == 1 & (INFO.PASS), .(n_unique_cnvs = length(unique(cnv_id))), by = .(mip_name, State)]
    mips_data_table[tmp, CVG.PASS := n_unique_cnvs > 1, on = .(mip_name, State)]
    return (mips_data_table)
}

filter_coverage_account_for_back_mutations <- function(mips_data_table) {
    logger$info("Filtering to include MIPs that fail the two independent gains/losses filter, but instead have back mutations")
    tmp = mips_data_table[mip_name %in% mips_data_table[cnv_id %like% "BM" & cnv_present_in_sample == 1 & (INFO.PASS), unique(mip_name)] & cnv_present_in_sample == 1,
                    .(n_unique_cnvs = length(unique(cnv_id))), by = .(mip_name)][n_unique_cnvs > 1]
    mips_data_table[mip_name %in% tmp[, mip_name] & cnv_present_in_sample==1 & (INFO.PASS) & CVG.PASS == FALSE, CVG.PASS := TRUE]
    return (mips_data_table)
}

mips_data <- filter_coverage_independent_gains_losses(mips_data)
mips_data <- filter_coverage_account_for_back_mutations(mips_data)

# Mips that occur within both a present cnv gain and a present cnv loss
# bm_mips <- mips_data[cnv_present_in_sample == 1 & (INFO.PASS), length(unique(State)), by = .(mip_name, sample)][V1>1]
# bm_mips[, .N, by = mip_name]

# gain_mips <- mips_data[State == "gain" & (INFO.PASS) & cnv_present_in_sample == 1][, .(cnv_id, mip_name, sample, State, seqnames, start, end, best_fit_genotype)]
# loss_mips <- mips_data[State == "loss" & (INFO.PASS) & cnv_present_in_sample == 1][, .(cnv_id, mip_name, sample, State, seqnames, start, end, best_fit_genotype)]
# View(gain_mips[loss_mips, on = .(mip_name, sample, seqnames), nomatch=0L])
