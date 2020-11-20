##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Plots the three tree figures - supplementary figure S1, figure 1A, figure 3A
#
# Usage:
# Rscript [this_file.R] iqtree_output_file.treefile
#
# Output:
# Three PDF files, with filenames generated from the input file name, i.e.
# if the input_file is tree.nwk, the outputs will be tree.nwk.supplementary.pdf
# tree.nwk.fig1A.pdf, tree.nwk.fig3A.pdf.


library(ggtree)
library(ggplot2)
library(treeio)
library(readxl)
library(data.table)
library(ape)
library(dftdLowCov)

haplotype_db <- load_sample_table()



# Samples to exclude
EXCLUDE <- c("8T2", "9T", "120T", "166T4", "390T1", "460T1", "478T2", "1070T1", "1072T1")
DUPLICATES <- haplotype_db[Cancer.type == "DFT1" & !is.na(Duplicate.or.time.course.biopsy)][Duplicate.or.time.course.biopsy == "Duplicate DNA extraction from same tumour", Sample_ID]

# Labels to italicise - these are samples low in phylogenetic information, which we indicate in
# the supplementary figure by italicising the label.
# N.B. ggtree can italicise labels, but they came out as a mixture of italicised and regular
# characters, so this script wraps the label in three underscores ('label' -> '___label___'),
# and I italicised them afterwards in Inkscape
#ITALICISE <- c("10T2", "55T3", "160T", "162T", "187T1", "188T1", "196T1", "196T2", "204T", "272T1", "272T2", "272T3", "463T1", "463T3", "469T1", "470T1", "475T2", "480T1", "626T1", "636T1", "826T2", "1009T1", "1078T4", "1079T1", "1079T2")
ITALICISE <- haplotype_db[Cancer.type == "DFT1" & `Included.in.tree.(Figure.1A,.Figure.S1)` %like% "italic", Sample_ID]

# Loads the tree and does a small amount of preprocessing to standardise the look of the tree
load_tree <- function(treefile) {
    tree <- read.newick(treefile)
    tree <- ape::drop.tip(tree, "reference_NA_NA")
    tree <- ape::ladderize(tree)
    tree
}

# Produces a table of tip labels from the tree
make_labels <- function(tree) {
    labels <- data.table(long_name = tree$tip.label)
    labels[, shortname := sapply(strsplit(long_name, "_"), `[[`, 1)]
    labels[, clade := ""]
    labels[long_name %like% "Clade_A1", clade := "A1"]
    labels[long_name %like% "Clade_A2", clade := "A2"]
    labels[long_name %like% "Clade_B", clade := "B"]
    labels[long_name %like% "Clade_C", clade := "C"]
    labels[long_name %like% "Clade_D", clade := "D"]
    labels[long_name %like% "Clade_E", clade := "E"]
    labels
}

# Get the tree file from the commandline
treefile = commandArgs(trailingOnly = TRUE)[1]
tree <- load_tree(treefile)
labels <- make_labels(tree)

labels[haplotype_db, long_name := paste(long_name, Year.of.sampling, sep = "_"), on = c("shortname" = "Sample_ID")]

# Reduce the tip label to just the sample ID (removes location and date info)
tree$tip.label <- sapply(strsplit(tree$tip.label, "_"), `[[`, 1)

# Attempt to define the base node of each clade via a pair of divergent
# samples within the clade/subclade. The MRCA of the pair should be a
# good proxy for the basal node of the clade.
pairs <- list(
    A1 = c("52T2", "136T"),
    A2 = c("381T1", "180T1"),
    B = c("826T1", "1059T2"),
    C = c("464T1", "46T1"),
    D = c("101T2", "102T2")
)

clade_e <- length(tree$tip.label) + 1
clade_d <- MRCA(tree, pairs[["D"]])
clade_c <- MRCA(tree, pairs[["C"]])
clade_b <- MRCA(tree, pairs[["B"]])
clade_a1 <- MRCA(tree, pairs[["A1"]])
clade_a2 <- MRCA(tree, pairs[["A2"]])



colours <- dftdLowCov::clade_colours()

# Put long tip labels back in (add back location and date info)
tree$tip.label <- labels[data.table(shortname=tree$tip.label), long_name, on = "shortname"]

# Enable colouring of tree by clade
tree_by_clades <- groupClade(tree, c(E=clade_e, D=clade_d, C=clade_c, B=clade_b, A1=clade_a1, A2=clade_a2))
ggtree(tree_by_clades, aes(color=group)) + scale_colour_manual(values = colours)

# Support values are encoded as "aLRT/Bootstrap". Here I drop the aLRT half and
# keep the bootstrap.
if (any(grepl("/", tree_by_clades$node.label))) {
    tree_by_clades$node.label <- tstrsplit(tree_by_clades$node.label, "/")[[2]]
}

####################
# Supplementary plot
####################
offset <- 0.0065
text_offset <- 0.001
supplementary_plot <-
    ggtree(tree_by_clades, aes(colour=group), size = 0.5) +
    scale_colour_manual(values = colours) +
    #geom_tiplab(size = 1.0, colour="black", align=FALSE, offset = 0.0001, linetype = NA) +
    geom_nodelab(aes(x=branch), vjust=-1, size=0.5, colour = "gray40", fontface=3) +
    labs(colour = "Clade") +
    geom_cladelabel(clade_a1, label = "A1", offset = offset, offset.text = text_offset, colour = colours["A1"]) +
    geom_cladelabel(clade_a2, label = "A2", offset = offset, offset.text = text_offset, colour = colours["A2"]) +
    geom_cladelabel(clade_b, label = "B", offset = offset, offset.text = text_offset, color = colours["B"]) +
    geom_cladelabel(clade_c, label = "C", offset = offset, offset.text = text_offset, color = colours["C"]) +
    geom_cladelabel(clade_d, label = "D", offset.text = -0.0031, barsize = NA, color = colours["D"]) +
    geom_cladelabel(clade_e-1, label = "E", offset.text = -0.0016, color = colours["E"]) +
    ylim(1, length(tree_by_clades$tip.label)) +
    theme(plot.margin = unit(-c(15,-5,15,0), units="mm"))

mydf <- labels[, .(tip.label=long_name, clade=clade, shortname)]
mydf[, is_wukalina := grepl("[Ww]ukalina", tip.label)]
mydf[, sz := as.numeric(ifelse(is_wukalina, 1, NA))]
# mydf[, tiplab := paste0('`', tip.label, '`')]
mydf[, tiplab := tip.label]
# mydf[shortname %in% ITALICISE, tiplab := paste0('italic(', tiplab, ')')]
mydf[shortname %in% ITALICISE, tiplab := paste0('___', tiplab, '___')]
mydf[tip.label %like% "Clade_A1", tiplabcolour := colours['A1']]
mydf[tip.label %like% "Clade_A2", tiplabcolour := colours['A2']]
mydf[tip.label %like% "Clade_B", tiplabcolour := colours['B']]
mydf[tip.label %like% "Clade_C", tiplabcolour := colours['C']]
mydf[tip.label %like% "Clade_D", tiplabcolour := colours['D']]
mydf[tip.label %like% "Clade_E", tiplabcolour := colours['E']]

supplementary_plot <- supplementary_plot %<+% mydf +
    geom_tippoint(pch = '*',
                  size = mydf$sz * 2,
                  fill = "black",
                  colour = "black",
                  show.legend = FALSE) +
    geom_tiplab(aes(label = tiplab), parse = FALSE,
                size = 1.0, colour="black", align=FALSE, offset = 0.0001, linetype = NA)


ggsave(plot=supplementary_plot,
       file=paste0(treefile, ".supplementary.pdf"),
       width = 210, height = 594, units = "mm")

################
# Main figure 1A
################
fig1A <-
    ggtree(tree_by_clades, aes(colour=group), size = 0.35) +
    scale_colour_manual(values = colours) +
    geom_tippoint(aes(colour=group), size = 0.75) +
    ylim(1, length(tree_by_clades$tip.label)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), units="mm"),
          legend.position = "none",
          plot.background = element_blank(),
          panel.background = element_blank())

ggsave(plot=fig1A %<+% mydf + geom_tippoint(pch = '*',
                                            size = mydf$sz * 2,
                                            fill = "black",
                                            colour = "black",
                                            show.legend = FALSE),
       file=paste0(treefile, ".fig1A.pdf"),
       width = 64, height = 111, units = "mm")

################
# Main figure 3A
################
labels[haplotype_db, is_cell_line := (Tissue.or.cell.line != "Tissue"), on = c("shortname" = "Sample_ID")]
mydf <- labels[, .(tip.label=shortname, clade=clade, is_cell_line=is_cell_line)]
mydf[, sz := as.numeric(ifelse(is_cell_line, 1.1, NA))]


t <- ggtree(tree_by_clades, size = 0.25) + scale_colour_manual(values = colours)
fig3A <- t %<+% mydf +
    geom_tippoint(pch = 21,
                  size = mydf$sz,
                  fill = "#f10f0e",
                  colour = "#af0f0e",
                  show.legend = FALSE) +
    ylim(0, 670) +
    labs(colour = "Clade") +
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 6))

ggsave(plot=fig3A,
       file=paste0(treefile, ".fig3A.pdf"),
       width = 1.8, height = 2.6)
