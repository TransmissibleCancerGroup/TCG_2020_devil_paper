##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 08 of the pipeline for analysing recurrent alleles.
# Refine genotyping after manual inspection of the data.


logger <- getLogger("PIPELINE_08")

logger$warn("The following genotyping step (assigning the implicated allele among back mutations)")
logger$warn("is done with manual intervention. If the data has changed, this step will need to be")
logger$warn("revised. Ensure the assignments are correct.") # up-to-date Nov-12-2019

do_exploration = FALSE
if (do_exploration) {
    # Find mips that intersect with back mutations, and pass filters (and are losses?)
    mips <- mips_data[(cnv_id %like% "BM") & seqnames == "4" & State == "loss" & (INFO.PASS) & (cnv_present_in_sample == 1), sort(unique(mip_name))]
    
    # Select rows from mips_data where the forward counterpart CNV to the back mutations is present in the sample (CNVs 9, 10, 3R1).
    # AND the back mutation is NOT present. This selects estimates of the genotype when the back mutation is not present, and we can take
    # this as a proxy for the ancestral state, i.e. the genotype that existed prior to the back mutation.
    fwd_selecter <- mips_data[cnv_present_in_sample == 1, .(has_fwd = any(cnv_id %in% c("9", "10", "3R1")), has_bkwd = any(cnv_id %like% "BM")), by = sample][(has_fwd == TRUE) & (has_bkwd == FALSE)]
    fwd <- mips_data[(FILTERS.PASS) & mip_name %in% mips][fwd_selecter, on = .(sample), nomatch=0L]
    
    # Split apart the estimated genotype into separate estimates of REF and ALT copy numbers
    fwd[, c("est_ref", "est_alt") := scan_col("([0-9.]+)\\|([0-9.]+)", estimated_genotype, c('n','n'), c("est_ref", "est_alt"))]
    
    # Count how frequently each best fit genotype occurs, per MIP. The frequency distribution will point towards the ancestral state.
    counts <- fwd[, .(.N, mean_ref = mean(est_ref), mean_alt = mean(est_alt)), by = .(mip_name, best_fit_genotype)][order(mip_name, -N)]
    summary_counts <- counts[, .(N, maxN = max(N), best_fit_genotype), by = .(mip_name)][N==maxN]
    
    # MANUAL INTERVENTION - these MIPs point to ancestral states of 2|1 (REF|ALT).
    # not included are Chr4_supercontig_000000126:67219_(Clade_A2_small_groups)_0463 (3|0, not CN informative)
    # and Chr4_supercontig_000000259:2267971_(Clade_A1_Clade_A2)_0485 (1.5|1.5, ambiguous/mis-called CNV)
    informative_mips <- c(
        "Chr4_supercontig_000000080:939540_(Clade_A2_Tarraleah_exclusive)_0456",
        "Chr4_supercontig_000000093:2874221_(Clade_A2_DFT1_0_WPP-DFT1_1)_0457",
        "Chr4_supercontig_000000104:302458_(Clade_A1_Channel-Forestier-LB)_0460",
        "Chr4_supercontig_000000113:4263166_(Clade_A1)_0461",
        # "Chr4_supercontig_000000126:67219_(Clade_A2_small_groups)_0463",
        "Chr4_supercontig_000000147:190559_(Clade_A2_small_groups)_0467",
        # "Chr4_supercontig_000000259:2267971_(Clade_A1_Clade_A2)_0485",
        "Chr4_supercontig_000000268:70223_(Clade_A1_Channel-Forestier)_0486",
        "Chr4_supercontig_000000276:2159131_(Clade_B_small_groups)_0489",
        "Chr4_supercontig_000000299:1103698_(Clade_C)_0491"
    )
    
    # 459 and 462 deleted because majority of estimated genotypes are ambiguous, meaning any REF or ALT calls
    # are likely to be noise. 463 and 485 deleted for reasons above.
    mips_to_delete <- c("Chr4_supercontig_000000100:672798_(Clade_B_DFT1_2_Bronte-DFT1_2_WPP)_0459",
                        "Chr4_supercontig_000000114:868144_(Clade_A2_DFT1_0_WPP_exclusive)_0462",
                        "Chr4_supercontig_000000126:67219_(Clade_A2_small_groups)_0463",
                        "Chr4_supercontig_000000259:2267971_(Clade_A1_Clade_A2)_0485")
    
    # Now look at each informative MIP in turn, in the samples that DO carry a back mutation, to see
    # how the genotype has changed from the inferred 2|1 state
    bkwd <- mips_data[(CVG.PASS) & (INFO.PASS) & cnv_present_in_sample == 1 & cnv_id %like% "BM" & mip_name %in% informative_mips]
    
    
    f<-fwd[mip_name %in% informative_mips & cnv_present_in_sample == 1, .(mean_ref = mean(est_ref), mean_alt = mean(est_alt), mean_diff = mean(diff)), by = mip_name]
    b<-bkwd[mip_name %in% informative_mips & cnv_present_in_sample == 1, .(mean_ref = mean(est_ref), mean_alt = mean(est_alt), mean_diff = mean(diff)), by = mip_name]
    f[b, on = .(mip_name)][order(mip_name)]
} else {
    # 459 and 462 deleted because majority of estimated genotypes are ambiguous, meaning any REF or ALT calls
    # are likely to be noise. 463 and 485 deleted for reasons above.
    mips_to_delete <- c("Chr4_supercontig_000000100:672798_(Clade_B_DFT1_2_Bronte-DFT1_2_WPP)_0459",
                        "Chr4_supercontig_000000114:868144_(Clade_A2_DFT1_0_WPP_exclusive)_0462",
                        "Chr4_supercontig_000000126:67219_(Clade_A2_small_groups)_0463",
                        "Chr4_supercontig_000000259:2267971_(Clade_A1_Clade_A2)_0485")
}

do_plot = FALSE
if (do_plot & do_exploration & require(cnpipe)) {
    # Exploratory plotting
    library(cnpipe)
    pdf(file.path(PIPELINE.DIR, "intermediate_data", "backmut_corrected_vafplots_2019_11_12.pdf"), height = 10, width = 20)
    par(mfrow = c(2,4))
    for (mip in informative_mips) {
        
        fwd_raw_vaf <- fwd[mip_name == mip & cnv_present_in_sample == 1, vaf]
        fwd_vaf <- apply(fwd[mip_name == mip & cnv_present_in_sample == 1,
                             .(total - alt, alt, 100, 1, idealised_logr(2, 1, 2, 0, clamp(purity, 0.01, 0.99), 2.0, 2.0), clamp(purity, 0.01, 0.99))],
                         1, function(row) v.estimate_tumour_vaf(row[[1]], row[[2]], row[[3]], row[[4]], row[[5]], row[[6]]))
        
        bkwd_raw_vaf <- bkwd[mip_name == mip & cnv_present_in_sample == 1, vaf]
        bkwd_vaf <- apply(bkwd[mip_name == mip & cnv_present_in_sample == 1,
                               .(total - alt, alt, 100, 1, idealised_logr(2, 1, 2, 0, clamp(purity, 0.01, 0.99), 2.0, 2.0), clamp(purity, 0.01, 0.99))],
                          1, function(row) v.estimate_tumour_vaf(row[[1]], row[[2]], row[[3]], row[[4]], row[[5]], row[[6]]))
        
        bkwd_labels <- bkwd[mip_name == mip & cnv_present_in_sample == 1, sample]
        
        FWD_Y <- fwd_vaf
        BKWD_Y <- bkwd_vaf
        
        plot(rep(0.5, length(FWD_Y)), FWD_Y, pch = 20, col = rgb(0,0,1,0.5), ylab = "corrected VAF", xaxt = "n", xlab = NA, ylim = c(0, 1), xlim = c(0, 1), main = mip, cex.main = 0.85)
        
        abline(h = FWD_Y, col = rgb(0,0,1,0.1))
        
        points(rep(0.5, length(BKWD_Y)), BKWD_Y, pch = 1, cex = 2, col = "red", lwd = 3)
        text(x = 0.5 + 0.1 * (-1)^seq_len(length(BKWD_Y)), y = BKWD_Y, labels = bkwd_labels, col = "red", cex = 0.7)
    }
    dev.off()
}


# Based on the exploration and the plots, I conclude the following changes:

logger$info("Removing inconsistent MIPs")
mips_data[(FILTERS.PASS) & mip_name %in% mips_to_delete, FILTERS.PASS := FALSE]

logger$info("Applying manually selected back mutation allele designations")
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000080:939540_(Clade_A2_Tarraleah_exclusive)_0456" & cnv_id %like% "BM", implicated_allele := "REF"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000093:2874221_(Clade_A2_DFT1_0_WPP-DFT1_1)_0457" & cnv_id %like% "BM", implicated_allele := "REF"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000104:302458_(Clade_A1_Channel-Forestier-LB)_0460" & cnv_id %like% "BM", implicated_allele := "REF"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000113:4263166_(Clade_A1)_0461" & cnv_id %like% "BM" & sample == "996T1", implicated_allele := "AMBIGUOUS"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000113:4263166_(Clade_A1)_0461" & cnv_id %like% "BM" & sample != "996T1", implicated_allele := "REF"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000147:190559_(Clade_A2_small_groups)_0467" & cnv_id %like% "BM", implicated_allele := "ALT"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000268:70223_(Clade_A1_Channel-Forestier)_0486" & cnv_id %like% "BM" & sample != "993T1.1", implicated_allele := "REF"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000268:70223_(Clade_A1_Channel-Forestier)_0486" & cnv_id %like% "BM" & sample == "993T1.1", implicated_allele := "AMBIGUOUS"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000276:2159131_(Clade_B_small_groups)_0489" & cnv_id %like% "BM" & sample %in% c("615T1", "134T6", "700T2", "874T1"), implicated_allele := "REF"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000276:2159131_(Clade_B_small_groups)_0489" & cnv_id %like% "BM" & sample %in% c("614T1", "85T", "427T1", "331T"), implicated_allele := "ALT"]
mips_data[(FILTERS.PASS) & mip_name == "Chr4_supercontig_000000299:1103698_(Clade_C)_0491" & cnv_id %like% "BM", implicated_allele := "REF"]
