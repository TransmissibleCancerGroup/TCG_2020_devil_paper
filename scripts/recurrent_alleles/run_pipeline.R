##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Runs the analysis pipeline.


PIPELINE.DIR = "/Users/kg8/Documents/projects/DFTD_low_coverage/dftd_recurrent_alleles/recurrency_pipeline"
suppressWarnings(source(file.path(PIPELINE.DIR, "01_load_data.R")))
suppressWarnings(source(file.path(PIPELINE.DIR, "02_overlap_mips_with_cnvs.R")))
suppressWarnings(source(file.path(PIPELINE.DIR, "03_annotate_extra_info.R")))
suppressWarnings(source(file.path(PIPELINE.DIR, "04_filter_informative.R")))
suppressWarnings(source(file.path(PIPELINE.DIR, "05_add_manual_annotation_from_tree.R")))
suppressWarnings(source(file.path(PIPELINE.DIR, "06_coverage_filter.R")))
suppressWarnings(source(file.path(PIPELINE.DIR, "07_genotyping.R")))
fwrite(mips_data, file.path(PIPELINE.DIR, "intermediate_data", "mips_data_07.csv"))
suppressWarnings(source(file.path(PIPELINE.DIR, "08_manual_genotyping.R")))
fwrite(mips_data[(FILTERS.PASS) & !is.na(implicated_allele)], file.path(PIPELINE.DIR, "intermediate_data", "mips_data_08.csv"))
suppressWarnings(source(file.path(PIPELINE.DIR, "09_plot_chromosome_map.R")))
