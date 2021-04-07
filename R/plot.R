#' Plot Expected Errors Profile
#'
#' Plots the cumulative expected errors (EE) of a fastq file or
#' collection of fastq files over sequence positions. This
#' function extends \code{\link[dada2]{plotQualityProfile}}
#' as it appears in version 1.12.1 and makes use of its code.
#'
#' The red solid line represents the median cumulative EE, the red
#' dotted lines represent the 25th and 75th percentiles of
#' cumulative EE, and the blue-green line represents the mean
#' cumulative EE.
#'
#' @param fn Full name(s) of input fastq file(s), including directory.
#' @param n (Optional) Default 5e+5. The number of reads to sample when processing fastq files; passed to \code{\link[ShortRead]{qa}}.
#' @param aggregate (Optiona) Default FALSE. If inputting multiple fastq files, whether to aggregate their expected errors into one summary plot. If "aggregate" and "include_quality_profile" are both TRUE, also aggregates the quality profiles into one summary plot.
#' @param logEE (Optional) Default TRUE. Whether to log10-transform the Y-axis of the expected error profile.
#' @param include_quality_profile (Optional) Default FALSE. Whether to include the quality profile in the return value (see below).
#'
#' @return If include_quality_profile = FALSE (default), returns a ggplot of the expected errors. If include_quality_profile = TRUE, returns a length of list two: (1) a ggplot of the expected errors, and (2) a ggplot of the quality profile(s).
#' @export
#'
#' @examples
#' \dontrun{
#' fnFs <- c("sample1_R1.fastq", "sample2_R1.fastq", "sample3_R1.fastq")
#' plotEEProfile(fnFs) # produces a multi-faceted ggplot
#' plotEEProfile(fnFs, aggregate=TRUE) # produces a single-faceted ggplot
#' p <- plotEEProfile(fnFs, aggregate=TRUE, include_quality_profile=TRUE) # produces a list containing two single-faceted ggplots:
#' plot(p[[1]]) # expected errors profile
#' plot(p[[2]]) # quality profile
#' }
#'
#' @importFrom magrittr "%>%"
plotEEProfile <- function(fn, n=5e+5, aggregate=FALSE, logEE=TRUE, include_quality_profile=FALSE) {
  if(length(fn) > 20 & missing(aggregate)) {
    message("More files than recommended for individual plotting. Setting aggregate = TRUE. ",
            "Disable by explicitly setting aggregate = FALSE.")
    aggregate <- TRUE
  }
  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0),
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0),
                       Cum = numeric(0), file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0),
                      rclabel = character(0), rc = numeric(0), file = character(0))
  FIRST <- TRUE
  message("Calculating summary statistics for plotting...")
  progressbar <- txtProgressBar(min = 0, max = length(fn[!is.na(fn)]), style = 3)
  for (i in 1:length(fn[!is.na(fn)])) {
    f <- fn[!is.na(fn)][i]
    srqa <- qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    } else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score,
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score,
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score,
                                                     foo$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count),
               simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s),
                         names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = basename(f))
      FIRST <- FALSE
    } else {
      plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
    }
    n.min <- min(n, rc)
    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)),
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                                       Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/n.min, #instead of /rc
                                       file = basename(f)))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score),
                                     label = basename(f), rclabel = rclabel, rc = rc,
                                     file = basename(f)))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  anndf$minScore <- min(anndf$minScore)

  ee_label0 <- data.frame(ee=2^seq(0,5), label=paste0("maxEE = ", 2^seq(0,5)))

  # If aggregating files
  if (aggregate) {
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf,
                                sum)
    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    means <- rowsum(plotdf.summary$Score * plotdf.summary$Count,
                    plotdf.summary$Cycle)/rowsum(plotdf.summary$Count,
                                                 plotdf.summary$Cycle)
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score,
                                                                             foo$Count, 0.25), simplify = TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score,
                                                                             foo$Count, 0.5), simplify = TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score,
                                                                             foo$Count, 0.75), simplify = TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count),
               simplify = TRUE)
    n.min <- ifelse(anndf$rc < n, anndf$rc, n)
    n.agg <- mean(n.min)
    statdf.summary <- data.frame(Cycle = as.integer(rownames(means)),
                                 Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                                 Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/n.agg)

    if(include_quality_profile==TRUE) {
      p_q <- ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) +
        geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", high = "black") +
        geom_line(data = statdf.summary, aes(y = Mean), color = "#66C2A5") +
        geom_line(data = statdf.summary, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        geom_line(data = statdf.summary, aes(y = Q50), color = "#FC8D62", size = 0.25) +
        geom_line(data = statdf.summary, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        ylab("Quality Score") +
        xlab("Cycle") +
        annotate("text", x = 0, y = 0, label = sprintf("Total reads: %d", sum(anndf$rc)), color = "red", hjust = 0) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        guides(fill = FALSE) +
        facet_wrap(~label) +
        ylim(c(0, NA))
      if (length(unique(statdf$Cum)) > 1) {
        p_q <- p_q +
          geom_line(data = statdf.summary, aes(y = Cum/nrow(anndf)),
                    color = "red", size = 0.25, linetype = "solid") +
          scale_y_continuous(sec.axis = sec_axis(~. * 10, breaks = c(0, 100), labels = c("0%", "100%"))) +
          theme(axis.text.y.right = element_text(color = "red"),
                axis.title.y.right = element_text(color = "red"))
      }
    }

    eedf <- dplyr::mutate(statdf.summary, across(c(starts_with("Q"), Mean), cumulativeExpectedErrors))
    eedf <- tidyr::pivot_longer(eedf, c(starts_with("Q"), Mean), names_to="series", values_to="expected_errors")

    p_ee <- ggplot(eedf) +
      geom_line(aes(x=Cycle, y=expected_errors, col=series, linetype=series, size=series)) +
      scale_color_manual(values=c("#66C2A5", "#FC8D62", "#FC8D62", "#FC8D62")) +
      scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed")) +
      scale_size_manual(values=c(1, 0.5, 0.5, 0.5)) +
      labs(y="Cum. expected errors", col="") +
      guides(col=FALSE, linetype=FALSE, size=FALSE)

    # If not aggregating files
  } else {
    if(include_quality_profile) {
      p_q <- ggplot(data = plotdf, aes(x = Cycle, y = Score)) +
        geom_tile(aes(fill = Count)) +
        scale_fill_gradient(low = "#F5F5F5", high = "black") +
        geom_line(data = statdf, aes(y = Mean), color = "#66C2A5") +
        geom_line(data = statdf, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        geom_line(data = statdf, aes(y = Q50), color = "#FC8D62", size = 0.25) +
        geom_line(data = statdf, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        ylab("Quality Score") + xlab("Cycle") + theme_bw() +
        theme(panel.grid = element_blank()) + guides(fill = FALSE) +
        geom_text(data = anndf, aes(x = 0, label = rclabel,
                                    y = 0), color = "red", hjust = 0) +
        facet_wrap(~file) +
        ylim(c(0, NA))
      if (length(unique(statdf$Cum)) > 1) {
        p_q <- p_q + geom_line(data = statdf, aes(y = Cum), color = "red", size = 0.25, linetype = "solid") +
          scale_y_continuous(sec.axis = sec_axis(~. * 10, breaks = c(0, 100), labels = c("0%", "100%"))) +
          theme(axis.text.y.right = element_text(color = "red"),
                axis.title.y.right = element_text(color = "red"))
      }
    }

    statdf %>%
      tidyr::pivot_longer(c(starts_with("Q"), Mean), names_to="series", values_to="quality") %>%
      group_by(file, series) %>%
      arrange(Cycle) %>%
      dplyr::mutate(expected_errors = cumulativeExpectedErrors(quality)) ->
      eedf

    ee_label <- ee_label0[ee_label0$ee < max(eedf$expected_errors),]

    p_ee <- ggplot(eedf) +
      geom_line(aes(x=Cycle, y=expected_errors, col=series, linetype=series, size=series)) +
      scale_color_manual(values=c("#66C2A5", "#FC8D62", "#FC8D62", "#FC8D62")) +
      scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed")) +
      scale_size_manual(values=c(1, 0.5, 0.5, 0.5)) +
      labs(y="Cum. expected errors", col="") +
      facet_wrap(~file) +
      guides(col=FALSE, linetype=FALSE, size=FALSE)
  }

  ee_label <- ee_label0[ee_label0$ee < max(eedf$expected_errors),]
  if(nrow(ee_label) > 0) {
    p_ee <- p_ee +
      geom_text(data=ee_label, aes(x=0, y=ee, label=label), vjust=-0.2, hjust=0, col="grey50", size=3) +
      geom_hline(data=ee_label, aes(yintercept=ee), lwd=0.5, col="grey50")
  }

  if(logEE) {
    p_ee <- p_ee + scale_y_continuous(trans="log10")
  }

  if(include_quality_profile) {
    return(list(p_ee, p_q))
  } else {
    return(p_ee)
  }
}

cumulativeExpectedErrors <- function(Q) {
  cumsum(10^(-Q/10))
}
