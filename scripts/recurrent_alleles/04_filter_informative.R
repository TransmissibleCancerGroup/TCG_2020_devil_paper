##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 04 of the pipeline for analysing recurrent alleles.
# Filters for allele-specific copy number-informative MIPS SNVs


logger <- getLogger("PIPELINE_04")

# Currently I have a table that contains every combination of cnv_id, mip_name and sample, among the
# full set of cnvs, mips, and samples for which there is cnv or mip information
# I need to filter this down to 
# * IDs (mip,cnv,sample,State) in which the cnv is present in the sample
# * MIPs that are informative of allele-specific copy number
# * MIPs for which there are at least 2 combinations of mip+cnv among either state=loss or state=gain


logger$info("Within gains, marking allele-specific copy number-informative MIPs")
logger$info("(MIP is informative if it carries the ALT allele)")
mips_data[State == "gain", informative := mip_allele == "ALT"]


logger$info("Within losses, the decision as to whether the mip site is ASCN informative is made by inspecting the tree")
logger$info("The first step is to identify potentially informative MIPs by looking at the high-level clade designation")


# Annotate clade prefix
tmp <- mips_data[State == "loss", .(cnv_id, mip_name, sample, clade, mip_allele, State)]
tmp2 <- unique(tmp[, .(clade)])
tmp2[, c("clade_prefix", "clade_suffix") :=
        scan_col("(Clade_\\w+)[-]*(\\w+)*", clade, c('c', 'c'))]
tmp[tmp2, clade_prefix := clade_prefix, on = .(clade)]
rm(tmp2)
tmp[, clade_prefix_contains_alt := any(mip_allele == "ALT"), by = .(mip_name, clade_prefix)]
mips_data[tmp, informative := clade_prefix_contains_alt, on = .(cnv_id, mip_name, sample, State, clade, mip_allele)]
rm(tmp)
logger$info("Marked potentially informative losses (informative=TRUE) based on the cnv+mip+sample combination's mip site carrying the ALT allele at the clade level")

treeplotfile <- file.path(PIPELINE.DIR, "intermediate_data", "informative_losses_trees_new_tree.pdf")
tree <- ladderize(root(keep.tip(tree, c("reference", unique(mips_data$sample))), "reference"))

logger$info(sprintf("Plotting tree figure to %s", treeplotfile))

pdf(treeplotfile, width = 10, height = 30)
for (my_mip in mips_data[State == "loss" & (informative), sort(unique(mip_name))]) {
    data <- mips_data[State == "loss" & mip_name == my_mip & (informative) & cnv_present_in_sample == 1,
                      .(cnvs = paste(ifelse(cnv_id == "", "NoLabel", cnv_id), collapse = ",")),
                      by = sample]
    if (data[,.N] == 0) next
    
    plot(tree, label.offset = .0015, cex = .45, main = my_mip, cex.main = 1.0)
    
    # Red spot wherever a sample has an ALT allele
    # Find all ALT samples in table
    table_samples_with_alt <- mips_data[mip_name == my_mip & mip_allele == "ALT" & State == "loss", unique(sample)]
    # Limit these to only samples that are also present in the tree
    tree_samples <- intersect(tree$tip.label, table_samples_with_alt)
    # Look up the tip numbers of these samples
    tipnums <- which(tree$tip.label %in% tree_samples)
    # Draw a red spot at each of these samples
    tiplabels(pch=19, tip = tipnums, col = "red")
    
    # Blue ring wherever a mip-cnv intersection occurs
    tipnums <- which(tree$tip.label %in% intersect(data[, sample], tree$tip.label))
    treetipnums <- data.table(sample = tree$tip.label, id = seq_along(tree$tip.label))
    data[treetipnums, id := id, on = .(sample)]
    tiplabels(pch=21, tip = data[!is.na(id), id], col = "blue", bg = NA, lwd = 2)
    if (data[!is.na(id), .N] > 0) {
        tiplabels(data[!is.na(id), cnvs], tip = data[!is.na(id), id], col = "blue", offset = 0.007, frame = "n", bg = NA, cex = 0.25)
    }
}

dev.off()
