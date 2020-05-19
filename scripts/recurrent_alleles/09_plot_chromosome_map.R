##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# Step 09 of the pipeline for analysing recurrent alleles.
# Plot the figure.


logger <- getLogger("PIPELINE_09")

#' Calculates the depth at
coverage_histogram <- function(dt, binsize = 100000) {
    bins <- data.table(start = seq(0, dt[, max(end)] - binsize, binsize))
    bins[, end := start + binsize - 1]
    setkey(bins, start, end)
    setkey(dt, NULL)
    foverlaps(dt, bins, nomatch = 0L)[, .(N = length(unique(cnv_id))), by = .(start, end)]
}

#' A faster algorithm to pack aligned reads into layers. Doesn't need IRanges.
rectangle_packing4 <- function(dt, sortby = "width", stack_per_sample = FALSE) {
    if (nrow(dt) == 0) {
        dt[, layer := 0]
        return (dt)
    }
    
    # Function to assign the lowest available layer, avoiding
    # any of the layers in `disallowed`
    assignLayer <- function(disallowed) {
        if (length(disallowed) == 0 || (length(disallowed) == 1 && disallowed == 0)) return (1)
        M <- max(disallowed)
        candidate <- setdiff(1:M, disallowed)
        ifelse(length(candidate) == 0, M+1, candidate)
    }
    
    dt[, origOrder := .I]
    setkey(dt, seqnames, start, end)
    if (stack_per_sample) {
        dt[, group := .GRP, by = .(cnv_id, sample)]
    } else {
        dt[, group := .GRP, by = .(cnv_id)]
    }
    
    hits <- foverlaps(dt, dt)[(group != i.group)][, .(queryHits = group, subjectHits = i.group)]
    
    if (length(sortby) == 1 && sortby == "width") {
        setorder(dt, -width)
    } else {
        if (length(intersect(sortby, colnames(dt))) == length(sortby)) {
            cat(sprintf("Sorting by %s\n", paste(sortby, collapse = ", then ")))
            setorderv(dt, sortby)
        } else {
            setorder(dt, start)
        }
    }
    dt[, layer := as.integer(0)]
    
    # No conflicts possible for first segment, so place in bottom layer
    dt[1, layer := 1]
    
    # processed <- 1
    N <- dt[, .N]
    for (currID in dt[2:.N, group]) {
        # Search for any segments that overlap the current segment. If any layers are assigned
        # for the hits, then they are disallowed for our current segment.
        disallowed_layers <- dt[group %in% hits[queryHits == currID, subjectHits], unique(layer)]
        assigned_layer <- assignLayer(disallowed_layers)
        dt[group == currID, layer := assigned_layer]
        
        # Just for counting and showing progress
        # processed <- processed + 1
        # print(paste("(", processed, "/", N, ")"))
    }
    
    # Put the dt back how we found it
    setorder(dt, origOrder)
    dt[, origOrder := NULL]
    dt[1:.N]
}


#' Scale a vector to min(v)=0 and max(v)=1
scale01 <- function(v) {
    v <- v - min(v)
    v/max(v)
}

make_curve <- function(x1, y1, x2, y2, ...){
    ys <- seq(y1, y2, length.out = 41)
    unscaled_xs <-
        if (y1 < y2) {
            plogis(ys, scale = 0.1 * abs(y2 - y1), loc = (y1 + y2) /2)
        } else {
            plogis(rev(ys), scale = 0.1 * abs(y2 - y1), loc = (y1 + y2) /2)
        }
    xs <- scale01(unscaled_xs) * (x2-x1) + x1
    list(x=xs, y=ys)
}

draw_curve <- function(curve, ...) {
    lines(curve$x, curve$y, ...)
}

x_transform <- function(xvalue, xscale, chromwidth, chrombuffer, state) {
    chrom_spacing <- (chromwidth / 2 + chrombuffer - 0.5) * if (state == "loss") {-1} else {1}
    xscale * (xvalue + chrom_spacing)
}


.plot_handle_options <- function(plot_options = NULL) {
    default_plot_options <- list(
        GAINS.SEGMENT.COLOUR  = "#EA6A58",
        LOSSES.SEGMENT.COLOUR = "#94A6F1",
        GAINS.HIST.COLOUR     = "pink",
        LOSSES.HIST.COLOUR    = "skyblue",
        CHROM.FILL.COLOUR     = "#E1E1E1",
        CHROM.BORDER.COLOUR   = "black",
        CHROM.STRIPE.COLOUR   = "black",
        CHROM.WIDTH           = 2,
        CHROM.BUFFER          = 1,
        SEGMENT.LWD           = 0.5,
        MIN.X                 = -40,
        MAX.X                 = 40,
        X.SCALE               = 1.0,
        CNV.TEXT.CEX          = 0.2,
        MIP.POINT.SHAPE       = 22,
        MIP.POINT.SIZE        = .5,
        MIP.POINT.LWD         = 0.3,
        CONNECTOR.LWD         = 0.5,
        CONNECTOR.COLOUR      = "grey40"
    )
    
    if (is.null(plot_options)) {
        return(default_plot_options)
    }
    
    for (key in names(default_plot_options)) {
        if (!(key %in% names(plot_options))) {
            plot_options[[key]] <- default_plot_options[[key]]
        }
    }
    return(plot_options)
}

.plot_make_depth_histogram <- function(dt) {
    # Depth histograms will be plotted underneath the segments
    if (nrow(dt) > 0) {
        histo <- coverage_histogram(dt)
    } else {
        histo <- data.table(start=1, end=1, N=1)[0]
    }
    histo
}

.plot_open_new_plot_with_dimensions <- function(min_x, max_x, min_y, max_y) {
    plot(c(min_x, max_x),
         c(min_y, max_y),
         type = "n", ylab = "Position", xlab = "Coverage",
         xaxt = "n", yaxt = "n",
         bty = "n",
         xlim = c(min_x, max_x),
         ylim = c(min_y, max_y))
}


#' Plot an individual chromosome's map. Assumes dt.gains and dt.losses
#' have already been produced, including setting the plotting layers
#' for each segment with rectangle packing
plot_chromosome_map <- function(dt_gains, dt_losses, chr_lengths, plot_options = NULL, get_plot_options = FALSE, chr_override = NULL) {
    plot_options <- .plot_handle_options(plot_options)
    if (get_plot_options) {
        return (plot_options)
    }
    
    CHR <- unique(dt_gains$seqnames)
    if(length(CHR) == 0) {
        CHR <- chr_override
    } else {
        if(length(CHR) > 1) stop("More than one chromosome in input data")
        if(nrow(dt_losses) > 0 && unique(dt_losses$seqnames) != CHR) stop("Chromosomes do not match between tables")
    }
    
    # Depth histograms will be plotted underneath the segments
    histo_gains <- .plot_make_depth_histogram(dt_gains)
    histo_losses <- .plot_make_depth_histogram(dt_losses)
    
    # Make sure losses are measured as negative, so they go on left side of plot
    dt_losses[, layer := -abs(layer)]
    histo_losses[, N := -abs(N)]
    
    # Get chromosome length
    L <- chr_lengths[CHROM==paste0("Chr", tolower(CHR)), LENGTH]
    
    # Constrain segments to the end of the chromosome
    dt_gains[end > L, end := L]
    dt_losses[end > L, end := L]
    
    # Open a new plot
    max_loss_depth <- if (nrow(histo_losses) > 0) histo_losses[, min(N)] - 1 else -1
    max_gain_depth <- if (nrow(histo_gains) > 0) histo_gains[, max(N)] + 1 else 1
    xlim = c(plot_options$MIN.X, plot_options$MAX.X)
    ylim = c(ifelse(CHR %in% c("5", "6", "X"), -3.0e8, -6.2e8), 2e7)
    .plot_open_new_plot_with_dimensions(xlim[1], xlim[2], ylim[1], ylim[2])
    
    # Add a y-axis
    if (CHR %in% c("1", "5")) {
        axis(side = 2, at=seq(0, -(L + 1e8 - L %% 1e8), -1e8),
             labels=abs(seq(0, -(L + 1e8 - L %% 1e8), -1e8)), line = -2)
    }
    
    
    # Plot the chromosome
    chrom_halfwidth = plot_options$CHROM.WIDTH / 2 # Half the width of the chromosome rectangle
    chrom_buffer = plot_options$CHROM.BUFFER # additional buffer space between the chromosome and any plotted elements
    chrom_spacing = chrom_halfwidth + chrom_buffer # Total offset used to position plot elements away from central chromosome
    if (chrom_halfwidth == 0) {
        segments(0, L, lwd = 2, lend = 1,
                 col = plot_options$CHROM.BORDER.COLOUR)
    } else {
        rect(-chrom_halfwidth, 0, chrom_halfwidth, -L,
             border = NA,
             col = plot_options$CHROM.FILL.COLOUR)
        if (nrow (dt_gains) > 0) {
            segments(-chrom_halfwidth, unique(dt_gains[start < L & end <= L, c(-start, -end)]), chrom_halfwidth,
                     col = plot_options$CHROM.STRIPE.COLOUR, lend = 1)
        }
        if (nrow (dt_losses) > 0) {
            segments(-chrom_halfwidth, unique(dt_losses[start < L & end <= L, c(-start, -end)]), chrom_halfwidth,
                     col = plot_options$CHROM.STRIPE.COLOUR, lend = 1)
        }
        rect(-chrom_halfwidth, 0, chrom_halfwidth, -L,
             border = plot_options$CHROM.BORDER.COLOUR,
             col = NA)
        text(0, 2e7, CHR, cex = 1.2, font = 2)
    }
    
    X.SCALE <- plot_options$X.SCALE
    
    # Gains histogram
    if (nrow (histo_gains) > 0) {
        rect(X.SCALE * (chrom_spacing),
             -histo_gains$start,
             X.SCALE * (chrom_spacing + histo_gains$N),
             -histo_gains$end,
             col = plot_options$GAINS.HIST.COLOUR, border = NA)
    }
    
    # Losses histogram
    if (nrow(histo_losses) > 0) {
        rect(X.SCALE * (-chrom_spacing),
             -histo_losses$start,
             X.SCALE * (histo_losses$N - chrom_spacing),
             -histo_losses$end,
             col = plot_options$LOSSES.HIST.COLOUR, border = NA)
    }
    
    ALT.FG.COLOUR   <- "black"
    ALT.BG.COLOUR   <- "black"
    REF.FG.COLOUR   <- "black"
    REF.BG.COLOUR   <- "white"
    UNDEF.FG.COLOUR <- "lightgrey"
    UNDEF.BG.COLOUR <- "lightgrey"
    MIP.POINT.SHAPE <- plot_options$MIP.POINT.SHAPE
    MIP.POINT.SIZE  <- plot_options$MIP.POINT.SIZE
    MIP.POINT.LWD   <- plot_options$MIP.POINT.LWD

    # Draw connections between gain segments
    if (nrow(dt_gains) > 0) {
        unique_gains <- unique(dt_gains, by = c("cnv_id", "sample"))
        connections <- unique_gains[, .(x1=layer, x2=data.table::shift(layer, type="lead"), y1=-end, y2=-data.table::shift(start, type="lead")), by = sample][!is.na(x2)]
        connections[, c("x1", "x2") := .(x_transform(x1, X.SCALE, 2*chrom_halfwidth, chrom_buffer, "gain"),
                                         x_transform(x2, X.SCALE, 2*chrom_halfwidth, chrom_buffer, "gain"))]
        
        for (i in seq_len(connections[, .N])) {
            args <- as.list(connections[i, .(x1,y1,x2,y2)])
            curve <- do.call(make_curve, args)
            draw_curve(curve, col = "grey40", lty=1, lwd = .7)
        }
    
    
        # Gains segments
        if (nrow(unique_gains) > 0) {
            rect(X.SCALE * (unique_gains$layer - 0.8 + chrom_spacing),
                 -unique_gains$start,
                 X.SCALE * (unique_gains$layer - 0.2 + chrom_spacing),
                 -unique_gains$end, col = plot_options$GAINS.SEGMENT.COLOUR,
                 lwd = plot_options$SEGMENT.LWD)
            
            text(x = X.SCALE * (unique_gains$layer - 0. + chrom_spacing),
                 y = unique_gains[, -(start + end)/2],
                 labels = unique_gains[, cnv_id], srt = 90, cex = plot_options$CNV.TEXT.CEX, col = "black")
            #labels = unique_gains[, paste(cnv_id, sample, sep = " ")], srt = 90, cex = plot_options$CNV.TEXT.CEX, col = "black")
        }
    }
    
    # Draw connections between loss segments
    if (nrow(dt_losses) > 0) {
        unique_losses <- unique(dt_losses, by = c("cnv_id", "sample"))
        connections <- unique_losses[, .(x1=layer, x2=data.table::shift(layer, type="lead"), y1=-end, y2=-data.table::shift(start, type="lead")), by = sample][!is.na(x2)]
        connections[, c("x1", "x2") := .(x_transform(x1, X.SCALE, 2*chrom_halfwidth, chrom_buffer, "loss"),
                                         x_transform(x2, X.SCALE, 2*chrom_halfwidth, chrom_buffer, "loss"))]
        
        for (i in seq_len(connections[, .N])) {
            args <- as.list(connections[i, .(x1,y1,x2,y2)])
            curve <- do.call(make_curve, args)
            draw_curve(curve, col = plot_options$CONNECTOR.COLOUR, lty=1, lwd = plot_options$CONNECTOR.LWD)
        }
        
        
        # Losses segments
        if (nrow(unique_losses) > 0) {
            rect(X.SCALE * (unique_losses$layer + 0.8 - chrom_spacing),
                 -unique_losses$start,
                 X.SCALE * (unique_losses$layer + 0.2 - chrom_spacing),
                 -unique_losses$end, col = plot_options$LOSSES.SEGMENT.COLOUR,
                 lwd = plot_options$SEGMENT.LWD)
            text(x = X.SCALE * (unique_losses$layer + 1.0 - chrom_spacing),
                 y = unique_losses[, -(start + end)/2],
                 labels = unique_losses[, cnv_id], srt = 90, cex = plot_options$CNV.TEXT.CEX, col = "black")
            #labels = unique_losses[, paste(cnv_id, sample, sep = " ")], srt = 90, cex = plot_options$CNV.TEXT.CEX, col = "black")
        }
    }
    
    # Gains MIPs
    if ("mip_position" %in% colnames(dt_gains)) {
        points(unique(dt_gains[implicated_allele == "REF", .(X.SCALE * (layer - (0.5 - chrom_spacing)), -mip_position)]),
               pch = MIP.POINT.SHAPE, cex = MIP.POINT.SIZE, lwd = MIP.POINT.LWD, col = REF.FG.COLOUR, bg = REF.BG.COLOUR)
        points(unique(dt_gains[implicated_allele == "ALT", .(X.SCALE * (layer - (0.5 - chrom_spacing)), -mip_position)]),
               pch = MIP.POINT.SHAPE, cex = MIP.POINT.SIZE, lwd = MIP.POINT.LWD, col = ALT.FG.COLOUR, bg = ALT.BG.COLOUR)
        points(unique(dt_gains[!(implicated_allele %in% c("REF", "ALT")), .(X.SCALE * (layer - (0.5 - chrom_spacing)), -mip_position)]),
               pch = MIP.POINT.SHAPE, cex = MIP.POINT.SIZE, lwd = MIP.POINT.LWD, col = UNDEF.FG.COLOUR, bg = UNDEF.BG.COLOUR)
        #points(dt_gains[, .(layer -0.5 + chrom_spacing, -mip_position)], pch = "*", col = "goldenrod")
    }
    
    # Losses MIPs
    if ("mip_position" %in% colnames(dt_losses)) {
        points(unique(dt_losses[implicated_allele == "REF", .(X.SCALE * (layer + (0.5 - chrom_spacing)), -mip_position)]),
               pch = MIP.POINT.SHAPE, cex = MIP.POINT.SIZE, lwd = MIP.POINT.LWD, col = REF.FG.COLOUR, bg = REF.BG.COLOUR)
        points(unique(dt_losses[implicated_allele == "ALT", .(X.SCALE * (layer + (0.5 - chrom_spacing)), -mip_position)]),
               pch = MIP.POINT.SHAPE, cex = MIP.POINT.SIZE, lwd = MIP.POINT.LWD, col = ALT.FG.COLOUR, bg = ALT.BG.COLOUR)
        points(unique(dt_losses[!(implicated_allele %in% c("REF", "ALT")), .(X.SCALE * (layer + (0.5 - chrom_spacing)), -mip_position)]),
               pch = MIP.POINT.SHAPE, cex = MIP.POINT.SIZE, lwd = MIP.POINT.LWD, col = UNDEF.FG.COLOUR, bg = UNDEF.BG.COLOUR)
        #points(dt_losses[, .(layer + 0.5 - chrom_spacing, -mip_position)], pch = "*", col = "goldenrod")
    }
    
    axlim <- if (nrow(dt_losses) > 0) dt_losses[, min(layer)] else -1
    axlim <- axlim + 5 - (axlim %% 5)
    axis(side = 1,
         at = X.SCALE * (seq(axlim, 0, 5) - chrom_spacing),
         labels = abs(seq(axlim, 0, 5)),
         cex.axis = .8, lwd = .5,
         pos = -L - 2e7,
         col = plot_options$LOSSES.SEGMENT.COLOUR,
         col.ticks = plot_options$LOSSES.SEGMENT.COLOUR,
         col.axis = plot_options$LOSSES.SEGMENT.COLOUR)
    
    axlim <- if (nrow(dt_gains) > 0) dt_gains[, max(layer)] else 1
    axlim <- axlim - (axlim %% 5)
    axis(side = 1,
         at = X.SCALE * (seq(0, axlim, 5) + chrom_spacing),
         labels = seq(0, axlim, 5),
         pos = -L - 2e7,
         cex.axis = .8, lwd = 1,
         col = plot_options$GAINS.SEGMENT.COLOUR,
         col.ticks = plot_options$GAINS.SEGMENT.COLOUR,
         col.axis = plot_options$GAINS.SEGMENT.COLOUR)
}

# Plot the table using chromosome map function


#infotable <- as.data.table(readxl::read_excel("infotable_with_idealised_vaf_min_2_cnvs.xlsx"))

logger$warn("Loading data for plotting. This has been subject to MANUAL INSPECTION to select MIP/CNV combinations to retain for plotting.")
logger$warn("These combinations are selected to 1) Pass all filters; 2) Show unambiguously which allele of the MIP is affected by the CNV; 3) Still pass coverage threshold after ambiguous alleles are removed")
logger$warn("Ensure this selection is up to date with the active data.") # up-to-date Nov-12-2019

# Prepare infotable from mips_data
unambig <- copy(mips_data[(FILTERS.PASS) & !is.na(implicated_allele)][implicated_allele %in% c("ALT", "REF", "BORDERLINE_ALT", "BORDERLINE_REF")])
unambig <- filter_coverage_independent_gains_losses(unambig)
unambig <- filter_coverage_account_for_back_mutations(unambig)
infotable <- unambig[CVG.PASS == TRUE]
infotable[, plot := TRUE]

# Manual intervention to remove ALT calls in cvn 9, where majority of calls are REF
infotable[cnv_id == "9" & (grepl("(457|460)$", mip_name)) & implicated_allele == "ALT", plot := FALSE]
# ...remove whole chromosome events when they occur in the same sample as smaller overlapping events
infotable[cnv_id %in% c("599","985"), plot := FALSE]

#infotable <- as.data.table(readxl::read_excel("/Users/kg8/Documents/projects/DFTD_low_coverage/dftd_recurrent_alleles/recurrency_pipeline/intermediate_data/mips_data_08_filtered_flagged.xlsx"))

logger$info("Collapsing data to unique combinations of CNV/MIP/State/Allele")
plotdata <- unique(infotable[(plot)], by = c("mip_name", "cnv_id", "State", "implicated_allele"))
logger$info("Assigning remaining (filtered) borderline REF/ALT calls to nearest REF/ALT state")
plotdata[implicated_allele %like% "BORDERLINE", implicated_allele := sub("BORDERLINE_", "", implicated_allele)]
logger$info("Collapsing plot data to unique combinations of CNV/MIP/State/Allele")
plotdata <- unique(plotdata, by = c("cnv_id", "mip_name", "State", "implicated_allele"))

logger$info("Loading Devil chromosome lengths table")
lengths_filename <- file.path(PIPELINE.DIR, "supporting_data", "devil_chrom_lengths.tsv")
chr.lengths <- fread(lengths_filename)

today <- format(Sys.time(), "%Y_%m_%d")
plotfile <- file.path(PIPELINE.DIR, "intermediate_data", paste0("recurrency_chromosome_map_", today, ".pdf"))
logger$info(sprintf("Making final plot to %s", plotfile))

pdf(plotfile, width = 10, height = 5)
layout(matrix(c(1,2,3,4,1,2,3,4,5,6,7,8), ncol = 4, byrow=T))
par(mar = c(1,0.1,0.5,0.1), oma = c(1, 0, 4, 0))
dt.gains.cache <- list()
dt.losses.cache <- list()
for (CHR in c(1:6, "X")) {
    plotdata.gains <- copy(plotdata[seqnames == CHR & State == "gain" & !is.na(cnv_id)])
    plotdata.losses <- copy(plotdata[seqnames == CHR & State == "loss" & !is.na(cnv_id)])
    plotdata.gains[, end := end - 1]
    plotdata.losses[, end := end - 1]
    plotdata.gains[, width := end - start + 1]
    plotdata.losses[, width := end - start + 1]
    
    if (CHR %in% names(dt.gains.cache)) {
        dt.gains <- dt.gains.cache[[CHR]]
    } else {
        dt.gains <- rectangle_packing4(plotdata.gains, stack_per_sample = TRUE, sortby = c("end", "sample"))
        dt.gains.cache[[CHR]] <- dt.gains
    }
    
    if (CHR %in% names(dt.losses.cache)) {
        dt.losses <- dt.losses.cache[[CHR]]
    } else {
        dt.losses <- rectangle_packing4(plotdata.losses, stack_per_sample = TRUE, sortby = c("end", "sample"))
        dt.losses.cache[[CHR]] <- dt.losses
    }
    
    plot_chromosome_map(dt.gains, dt.losses, chr.lengths,
                        plot_options = list(CHROM.FILL.COLOUR = rgb(170/255, 170/255, 170/255),
                                            CHROM.STRIPE.COLOUR = NA,
                                            GAINS.HIST.COLOUR = NA,
                                            LOSSES.HIST.COLOUR = NA,
                                            X.SCALE = 8,
                                            CHROM.BUFFER = -2.5,
                                            CHROM.WIDTH = 8.0,
                                            MIN.X = -100,
                                            MAX.X = 100,
                                            CNV.TEXT.CEX = 0.3,
                                            MIP.POINT.SHAPE = 22,
                                            MIP.POINT.SIZE = 0.7,
                                            MIP.POINT.LWD = 0.35,
                                            CONNECTOR.LWD = 0.5,
                                            CONNECTOR.COLOUR = "grey50"),
                        chr_override = CHR)
}
plot(c(1,1), c(1, -1), pch = 22, col = "black", bg = c("white", "black"), cex = 2, bty = "n", ylim = c(-5, 5), xaxt = "na", yaxt = "na", ylab = NA, xlab = NA)
text(x = c(1.1, 1.1), y = c(1, -1), labels = c("REF", "ALT"))
mtext("Recurrent gains and losses", outer = TRUE, cex = 1.2)
dev.off()

plotdatafile <- file.path(PIPELINE.DIR, "intermediate_data", paste0("plotdata_", today, ".xlsx"))
loginfo(paste("Saving plotting data to", plotdatafile))
openxlsx::write.xlsx(plotdata, plotdatafile)

# intersect(colnames(infotable), colnames(mips_data))
# 
# chk <- function(N) {
#     losses <- dt.losses.cache[[N]]
#     gains <- dt.gains.cache[[N]]
#     if (is.null(gains) | nrow(gains) == 0) {
#         if (!is.null(losses) & nrow(losses) > 0) {
#             all <- losses
#         }
#         else {
#             return (data.table())
#         }
#     }
#     else if (is.null(losses) | nrow(losses) == 0) {
#         all <- gains
#     }
#     else {
#         all <- rbindlist(list(losses, gains))
#     }
#     dt <- unique(mips_data[(FILTERS.PASS), .(cnv_id, sample, State)])[unique(all[, .(cnv_id, sample, State, seqnames, start, end, layer)]), on = .(cnv_id, sample, State)][order(sample)]
#     #dt[dt[, .N, by = sample][N>1], on = "sample"]
#     dt
# }
