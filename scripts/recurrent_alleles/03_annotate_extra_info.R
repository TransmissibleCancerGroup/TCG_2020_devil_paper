##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 03 of the pipeline for analysing recurrent alleles.
# Adds an annotation: 
#   if a clade carries an ALT allele, this is recorded in the 'clades_alt' column


logger <- getLogger("PIPELINE_03")

tmp <- mips_data[cnv_present_in_sample==1 & mip_allele == "ALT", .(clades_alt = unique(clade)), by = .(mip_name, cnv_id)]
mips_data[, clades_alt := FALSE]
mips_data[tmp, clades_alt := TRUE, on = c("mip_name", "cnv_id", clade = "clades_alt")]
rm(tmp)
logger$info("Added clades_alt annotation - all clades that carry an ALT allele are marked clades_alt='TRUE'")

mips_data[cn_data, copy_number := copy_number, on = .(cnv_id, sample)]
logger$info("Added copy number annotation to `mips_data` - cnv_id-sample pairs are annotated with their estimated copy number")
