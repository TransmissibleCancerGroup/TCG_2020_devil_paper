##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 02 of the pipeline for analysing recurrent alleles.
# Finds MIPS SNVs that overlap with CNVs


stopifnot(length(mips_data[, unique(sample)]) == 449) # 449 samples with MIPs data

logger <- getLogger("PIPELINE_02")

mips_data[, c("start", "end") := mip_position]
setkey(cnv_data, seqnames, sample, start, end)
setkeyv(mips_data, key(cnv_data))
ol <- foverlaps(cnv_data, mips_data, nomatch=0L)
ol[, c("start", "end") := NULL]
setnames(ol, old = c("i.start", "i.end"), new = c("start", "end"))
setkey(ol, cnv_id, mip_name, sample, seqnames, start, end)
setcolorder(ol)
mips_data <- ol
stopifnot(length(mips_data[, unique(sample)]) == 418) # Only 418 samples with both MIP and CNV data

rm(ol)
logger$info("Overlapped MIPs with CNVs - data stored as `mips_data`")
