##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 05 of the pipeline for analysing recurrent alleles.
# Adds phylogenetic information - requires manual inspection of the tree.


logger <- getLogger("PIPELINE_05")

table_path <- file.path(PIPELINE.DIR, "intermediate_data", "tree_informative_losses.xlsx")
logger$warn(paste("Ensure manual information in", table_path, "is correct"))
logger$info("Reading precomputed table of informative mip+cnv_id combinations for loss cnvs into 'informative'")
manual_table <- as.data.table(read.xlsx(table_path, "final"))

manual_table[, State := "loss"]
manual_table[, informative := TRUE]

logger$info("Adding 'tree_informative' column based on manual annotation")
mips_data[, tree_informative := FALSE]
mips_data[manual_table, tree_informative := i.tree_informative, on = .(cnv_id, mip_name, State, informative)]


mips_data[, INFO.PASS := FALSE]
mips_data[State == "gain", INFO.PASS := informative]
mips_data[State == "loss", INFO.PASS := tree_informative]
mips_data[, informative := NULL]
mips_data[, tree_informative := NULL]
