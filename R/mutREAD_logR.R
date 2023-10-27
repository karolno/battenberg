# require(Biobase)
# require(BSgenome.Hsapiens.UCSC.hg19)
# require(viridis)
# require(cowplot)
# require(Rsubread)
# require(ggplot2)
# require(viridis)
# require(cowplot)
# require(GenomicRanges)
# require(gtools)

#' Process mutREAD data
#'
#' This script performs counting of coverage within the individual target regions of mutREAD data for both bam files. It then performs normalisation of the coverage withint these regions. The normalization is done across the specific constant size genomics bins (default is 50 kb). It also outputs some diagnostics figures and raw data.
#'
#' @param binspan The span of the genomic region used to perform averaging of the signal for LogR calculation. It is an integer number. (Default: 5e5L)
#' @param tumourbam path to the location of tumour bam file. This should be absolute path
#' @param normalbam path to the location of normal bam file. This should be absolute path
#' @param ref.sample The name of reference sample
#' @param tumour.sample The name of tumour sample
#' @param directory The output folder. This should be absolute path
#' @param genomebuild Genome version. Can be hg19 or hg38
#' @param bins The location of .rds file contain information about bins to be analysed
#' @param nthreads number of threads to use
#' @param segment Should mutREAD data be segmented using HMMCopy
#' @author Karol Nowicki-Osuch
#' @export
process.mutREAD <- function (binspan=5e5L, tumourbam, normalbam, ref.sample, tumour.sample, directory, bins, genomebuild = "hg19", nthreads = 1, segment = TRUE) {
  # Get all bins to be used for analysis
  # bins <-"~/Dropbox/Postdoc/git/mutREAD/Battenberg_mutREAD/data_files/bins50-800.rds"
  bins<-readRDS(file = bins)
  if (!requireNamespace("HMMcopy", quietly = TRUE)) {
    stop(
      "Package \"HMMcopy\" must be installed to use this function.",
      call. = FALSE
    )
  }

  require("GenomicRanges")
  require("IRanges")
  require("Biobase")
  if (genomebuild == "hg19") {
    require("BSgenome.Hsapiens.UCSC.hg19")
  } else if (genomebuild == "hg38") {
    require("BSgenome.Hsapiens.UCSC.hg38")
  } else {
    stop("This version of genome is not supported")
  }
  # copyNumbersSmooth.cancer <- process.mutREAD.bams(sample = "SLX-15782.C1", bamlocation = "/mnt/data/mutREAD/SLX-15782/bams_alt/", bins = bins, binspan = 5e5)
  # copyNumbersSmooth.normal <- process.mutREAD.bams(sample = "SLX-15782.A2", bamlocation = "/mnt/data/mutREAD/SLX-15782/bams_alt/", bins = bins, binspan = 5e5)
  # Read bam files to get the counts for cancer samples
  copyNumbersSmooth.cancer <- process.mutREAD.bams(sample.name = tumour.sample, bamlocation = tumourbam, bins = bins, directory = directory, genomebuild = genomebuild,  binspan = binspan, segment = segment, nthreads = nthreads)
  # Read bam files to get the counts for the reference sample
  copyNumbersSmooth.normal <- process.mutREAD.bams(sample.name = ref.sample, bamlocation = normalbam, bins = bins, directory = directory, genomebuild = genomebuild, binspan = binspan, segment = FALSE, nthreads = nthreads)

  # Correct the counts by normal sample
  copyNumbersSmooth <- copyNumbersSmooth.cancer
  copyNumbersSmooth[,"corrected.counts"] <- copyNumbersSmooth[,"corrected.counts"] - copyNumbersSmooth.normal[,"corrected.counts"]

  # Created data frame for segmentation and visualisation
  data.df <- as.data.frame(do.call(rbind,base::strsplit(rownames(copyNumbersSmooth), split = ":|-", perl = TRUE)))
  data.df$V1 <- factor(data.df$V1, levels=gtools::mixedsort(unique(data.df$V1)))
  data.df <- cbind(data.df, copyNumbersSmooth[,"corrected.counts"])
  colnames(data.df)<- c("chr", "start", "end",  "copy")
  # Segment the data
  seg.data <- HMMcopy::HMMsegment(data.df, verbose = FALSE)
  rleseg <- rle(paste0(data.df$chr, ":", seg.data$state))
  seg.medians <- rep(seg.data$segs$median, times = rleseg$lengths)
  breakpoints.out <- seg.data$segs[duplicated(seg.data$segs$chr), 1:2]
  colnames(breakpoints.out) <- c("chromosome", "position")
  write.table(breakpoints.out, paste0(directory, "/", tumour.sample, "_breakpoints.after.normal.tab"), quote = FALSE, sep = "\t", row.names = FALSE)
  # write.table(seg.data$segs[,1:3], paste0(directory, "/", tumour.sample, "_breakpoints.after.normal2.tab"), quote = FALSE, sep = "\t", row.names = FALSE)
  copyNumbersSmooth <- cbind(copyNumbersSmooth, "state" = seg.data$state, "state.median" = seg.medians)
  # copyNumbersSmooth[,"state"] <- seg.data$state
  # copyNumbersSmooth[,"state.median"] <- seg.medians

  # Print diagnostics figure
  coverage.plot.state.after.normal <- coverage.plot.mutREAD_2(copyNumbersSmooth, pct.plot = 100, run.median.k = 1) #+ ylim(-2,2)
  ggsave(paste0(directory, "/", tumour.sample, "_mutREAD.combined.coverage.segmented.after.normal.plot.png"), plot = coverage.plot.state.after.normal, device = "png", height = 6, width = 15)

  # add missing data
  # Create intervals of entire genome
  if (genomebuild == "hg19") {
    bins.intervals <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19),
                                 tilewidth=binspan,
                                 cut.last.tile.in.chrom=TRUE)
    seqlevelsStyle(bins.intervals) <- "NCBI"
    bins.intervals <- keepSeqlevels(bins.intervals, c(1:22, "X", "Y"), pruning.mode="tidy")
    names(bins.intervals)<-paste0(seqnames(bins.intervals), ":", ranges(bins.intervals))
  } else {

    bins.intervals <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg38),
                                 tilewidth=binspan,
                                 cut.last.tile.in.chrom=TRUE)
    # seqlevelsStyle(bins.intervals) <- "NCBI"
    bins.intervals <- keepSeqlevels(bins.intervals, paste0("chr",c(1:22, "X", "Y")), pruning.mode="tidy")
    names(bins.intervals)<-paste0(seqnames(bins.intervals), ":", ranges(bins.intervals))
  }

  # Process normalized data into the entire genome (including regions that were not initially counted)
  final.regions.g <- bins.intervals
  out.data <-  GenomicRanges::makeGRangesFromDataFrame(data.df, keep.extra.columns=TRUE)
  hits <-  GenomicRanges::findOverlaps(final.regions.g, out.data, ignore.strand=TRUE, select = "all", minoverlap = 1L)
  hitsByQuery <- as(hits, "List")
  final.regions.g$copynumber <- IRanges::median(IRanges::extractList(out.data$copy, hitsByQuery), na.rm = TRUE)

  #Read data from battenberg default processing of logR
  original.logR<-read.delim(paste0(directory, "/", tumour.sample, "_mutantLogR.tab"))
  original.logR<-cbind(original.logR[,1:2], original.logR[,2:3])
  colnames(original.logR) <- c("chromosome", "start", "end", colnames(original.logR)[4])
  original.logR.g <- GenomicRanges::makeGRangesFromDataFrame(original.logR, keep.extra.columns=TRUE)
  hits2 <- GenomicRanges::findOverlaps(original.logR.g, final.regions.g, ignore.strand=TRUE, select = "all", minoverlap = 1L)
  hitsByQuery2 <- as(hits2, "List")
  elementMetadata(original.logR.g)[,1] <- IRanges::median(IRanges::extractList(final.regions.g$copynumber, hitsByQuery2), na.rm = TRUE)

  for (chromosome in seqlevels(original.logR.g)) {
    tmp.data<-elementMetadata(original.logR.g)[original.logR.g@seqnames == chromosome,1]
    tmp.which<-which(is.na(tmp.data))
    tmp.which2<-which(!is.na(tmp.data))

    for(n in tmp.which) {
      tmp.data[n] <- median(c(tmp.data[rev(tmp.which2[tmp.which2<n])[1]],tmp.data[tmp.which2[tmp.which2>n][1]]), na.rm = TRUE)

    }
    elementMetadata(original.logR.g)[original.logR.g@seqnames == chromosome,1] <- tmp.data
  }

  original.logR[,4]<-elementMetadata(original.logR.g)[,1]
  original.logR <- original.logR[,c(1,2,4)]
  colnames(original.logR)[1:2] <- c("Chromosome", "Position")
  write.table(original.logR, paste0(directory, "/", tumour.sample, "_mutantLogR_gcCorrected.tab"), quote = FALSE, sep = "\t", row.names = FALSE)
}

#' Read bam files and perform counts corrections
#'
#' This function reads the bamfiles and performs counting reads in the target regions. It then performs gc and read length correction
#'
#' @param sample.name The name of reference sample
#' @param bamlocation path to the location of bam files. This should be absolute path
#' @param bins The AnnotatedDataFrame object containing information about bins to be analysed
#' @param binspan The span of the genomic region used to perform averaging of the signal for LogR calculation. It is an integer number. (Default: 5e5L)
#' @param directory The output folder. This should be absolute path
#' @param genomebuild Genome version. Can be hg19 or hg38
#' @param segment Boolean indicating if the data should be segmented (Default: TRUE)
#' @param nthread Number of threads to use
#' @author Karol Nowicki-Osuch
#' @export
process.mutREAD.bams<-function(sample.name, bamlocation, bins, binspan=5e5L, directory, genomebuild = "hg19", segment = TRUE, nthreads = 1) {
  if (!requireNamespace("HMMcopy", quietly = TRUE)) {
    stop(
      "Package \"HMMcopy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop(
      "Package \"Rsubread\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop(
      "Package \"Biobase\" must be installed to use this function.",
      call. = FALSE
    )
  }
  require("GenomicRanges")
  require("Biobase")
  # require("BSgenome.Hsapiens.UCSC.hg19")
  require("viridis")
  require("cowplot")
  if (genomebuild == "hg19") {
    require("BSgenome.Hsapiens.UCSC.hg19")
  } else if (genomebuild == "hg38") {
    require("BSgenome.Hsapiens.UCSC.hg38")
  } else {
    stop("This version of genome is not supported")
  }

  print(paste("Processing file:", sample.name))

  bins2<-bins@data[,1:3]
  # if (genomebuild == "hg19") {
  #   bins2$chromosome<-paste0("chr", bins2$chromosome)
  # }
  colnames(bins2)<-c("Chr", "Start", "End")
  bins2$GeneID<-rownames(bins2)
  bins2$Strand<-"."
  # perform counting
  fcounts <- Rsubread::featureCounts(bamlocation, annot.ext = bins2, fracOverlap = 0.75, minMQS = 37, ignoreDup = TRUE, isPairedEnd = TRUE, requireBothEndsMapped = TRUE, checkFragLength = TRUE, minFragLength = 40,  maxFragLength = 810, autosort = TRUE, nthreads = nthreads, tmpDir = "/tmp")

  corrected.data <- mutREADestimateCorrection(counts.data = fcounts$counts, bins.data = bins, length.bin = 10, maxIter = 10, variables = c("gc", "length"), span = 0.1, nthreads = nthreads)

  write.table(corrected.data, paste0(directory, "/", sample.name, "_mutREAD.region.counts.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)


  all.data.long <- aggregate(cbind(counts, fit) ~ gc + length, data = corrected.data, FUN = mean, na.rm = TRUE)
  # all.data.frequency <- aggregate(length ~ round(gc) + length.bin, data = pData(bins)[,4:10], length)
  # colnames(all.data.frequency) <- c("gc","length","frequency")

  density.raw <- ggplot(all.data.long, aes(x=length, y=gc, z = counts) ) +
    stat_contour(geom = "polygon", aes(fill = after_stat(level))) +
    geom_tile(aes(fill = counts)) +
    stat_contour(bins = 15) +
    xlab("length") +
    ylab("gc") +
    ylim(c(0,100)) +
    scale_fill_viridis(option="magma") +
    guides(fill = guide_colorbar(title = "Mean count")) +
    theme_minimal()
  ggsave(paste0(directory, "/", sample.name, "_mutREAD.densitycounts.plot.png"), plot = density.raw, device = "png", height = 10, width = 12)
  # print(density.raw)

  density.fit <- ggplot(all.data.long, aes(x=length, y=gc, z = fit) ) +
    stat_contour(geom = "polygon", aes(fill = after_stat(level))) +
    geom_tile(aes(fill = fit)) +
    stat_contour(bins = 15) +
    xlab("length") +
    ylab("gc") +
    ylim(c(0,100)) +
    scale_fill_viridis(option="magma") +
    guides(fill = guide_colorbar(title = "Mean fit")) +
    theme_minimal()
  ggsave(paste0(directory, "/", sample.name, "_mutREAD.densityfit.plot.png"), plot = density.fit, device = "png", height = 10, width = 12)

  # print(density.fit)


  length.plot <- ggplot(corrected.data, aes(x = length, y = counts, group = length)) +
    geom_boxplot() +
    scale_y_log10() +
    # ylim(c(0,50)) +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    ggtitle("Total counts vs read lenght") +
    xlab("Read length")+
    theme_minimal()

  ggsave(paste0(directory, "/", sample.name, "_mutREAD.length.plot.png"), plot = length.plot, device = "png", height = 10, width = 12)

  # print(length.plot)

  gc.plot <- ggplot(corrected.data, aes(x = gc, y = counts, group = gc)) +
    geom_boxplot() +
    scale_y_log10() +
    xlim(c(0,100)) +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    ggtitle("Total counts vs GC") +
    xlab("GC")+
    theme_minimal()
  ggsave(paste0(directory, "/", sample.name, "_mutREAD.gc.plot.png"), plot = gc.plot, device = "png", height = 10, width = 12)

  # print(gc.plot)

  if(! all(is.na(corrected.data$mappability))) {
    mappability.plot <- ggplot(corrected.data, aes(x = mappability, y = counts, group = mappability)) +
      geom_boxplot() +
      scale_y_log10() +
      xlim(c(0,100)) +
      stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
      ggtitle("Total counts vs mappability") +
      xlab("Mappability")+
      theme_minimal()
    ggsave(paste0(directory, "/", sample.name, "_mutREAD.mappibility.plot.png"), plot = mappability.plot, device = "png", height = 10, width = 12)
  }
  combined.data <- mutREADcombineData(corrected.data, bins = bins2, n.cores = nthreads, span = binspan, genomebuild = genomebuild )
  # combined.data <- mutREADcorrectData(corrected.data)

  # Center around the mean copy number state
  combined.data[, "corrected.counts"] <- combined.data[, "corrected.counts"] -  median(combined.data[!grepl(pattern = "X|Y", rownames(combined.data), perl = TRUE),"corrected.counts"], na.rm = TRUE)

  coverage.plot <- coverage.plot.mutREAD(combined.data,pct.plot = 100, run.median.k = 1) + geom_hline(yintercept = median(combined.data[,"corrected.counts"], na.rm = TRUE), color = "red")
  # coverage.plot

  ggsave(paste0(directory, "/", sample.name, "_mutREAD.combined.coverage.plot.png"), plot = coverage.plot, device = "png", height = 6, width = 15)



  if (segment) {
    data.df <- as.data.frame(do.call(rbind,base::strsplit(rownames(combined.data), split = ":|-", perl = TRUE)))
    data.df$V1 <- factor(data.df$V1, levels=gtools::mixedsort(unique(data.df$V1)))
    data.df <- cbind(data.df, combined.data[,"corrected.counts"])
    colnames(data.df)<- c("chr", "start", "end",  "copy")

    seg.data <- HMMcopy::HMMsegment(data.df, verbose = FALSE)
    rleseg <- rle(paste0(data.df$chr, ":", seg.data$state))
    seg.medians <- rep(seg.data$segs$median, times = rleseg$lengths)
    combined.data <- cbind(combined.data, state = seg.data$state, state.median = seg.medians)
    # combined.data[,"state.median"] <- seg.medians
    breakpoints.out <- seg.data$segs[duplicated(seg.data$segs$chr), 1:2]
    write.table(breakpoints.out, paste0(directory, "/", sample.name, "_breakpoints.tab"), quote = FALSE, sep = "\t", row.names = FALSE)



    coverage.plot.state <- coverage.plot.mutREAD_2(combined.data, pct.plot = 100, run.median.k = 1) #+ ylim(-2,2)
    # coverage.plot.state
    ggsave(paste0(directory, "/", sample.name, "_mutREAD.combined.coverage.segmented.plot.png"), plot = coverage.plot.state, device = "png", height = 6, width = 15)

  }

  return(combined.data)

}

#' Estimate GC and read length bias in the data
#'
#' This script performs estimation of bias within the data.
#'
#' @param counts.data data.frame contains counts data for each target region
#' @param bins.data The AnnotatedDataFrame object containing information about bins to be analysed
#' @param length.bin The binning value around which region lengths are normaliszed
#' @param span for loess
#' @param family for loess
#' @param maxIter number of interactions used for the reestimation of loess
#' @param method Method used for the calculation of central value during modeling
#' @param pseudo.count Integer value to be added to each target region counts. Useful when log transforming counts
#' @param variable One of c("gc", "mappability", "length")
#' @param genomebuild Genome version. Can be hg19 or hg38
#' @param emove.zeros Boolean value whether targer regions with zero count should be removed from analysis
#' @param log.data Boolean value whether data should be log transformed
#' @param segment Boolean value whether data should be segmented
#' @param ... Other parameters that can be used with loess function
#' @author Karol Nowicki-Osuch
#' @export
mutREADestimateCorrection <- function(counts.data, bins.data, length.bin = 10, span=0.65, family="symmetric",
                                      maxIter=1, cutoff=4.0, method = "mean", pseudo.count = 0,
                                      variables=c("gc", "mappability", "length"), genomebuild = "hg19", remove.zeros = FALSE, log.data = FALSE, segment = TRUE, ...
) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop(
      "Package \"Biobase\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop(
      "Package \"matrixStats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  require("Biobase")
  # require("BSgenome.Hsapiens.UCSC.hg19")
  if (genomebuild == "hg19") {
    require("BSgenome.Hsapiens.UCSC.hg19")
  } else if (genomebuild == "hg38") {
    require("BSgenome.Hsapiens.UCSC.hg38")
  } else {
    stop("This version of genome is not supported")
  }

  require("viridis")
  require("cowplot")
  # counts.data = fcounts$counts
  # bins.data = bins
  # variables=c("gc", "length")
  variables <- match.arg(variables, several.ok=TRUE)

  # Get the data that will be used for modeling
  counts <- counts.data

  # check if values with zeros should be removed
  if (remove.zeros) {
    excluded.data <- counts == 0
  } else {
    excluded.data <- rep(FALSE, times = length(counts))
  }

  # # Correct the counts data by pseudocount value
  working.counts <- counts + pseudo.count

  descriptions <- c(gc="GC content", mappability="mappability", length = "region lenght")
  print(paste("Calculating correction for ",
              paste(descriptions[variables], collapse=" and "), sep = ""))

  # Get the gc data
  gc <- round(pData(bins.data)$gc)

  # Get mappability data
  mappability <- round(pData(bins.data)$mappability)

  # Get the target site length and converted it to the indicated bin size
  length <-pData(bins.data)$length
  length <- floor(length/length.bin)*length.bin

  combined.variables.all <- data.frame(gc=gc, mappability=mappability, length=length)
  combined.variables <- combined.variables.all[!excluded.data,]
  # aggregate the data across all potential conditions
  if (log.data) {
    # work in log space
    working.counts <- log2(counts[!excluded.data])
  } else {
    working.counts <- counts[!excluded.data]
  }

  if (method == "mean") {
    all.data.long<-aggregate(
      x = working.counts,
      by = list(gc=combined.variables$gc, mappability=combined.variables$mappability, length=combined.variables$length)[variables],
      FUN = mean,
      na.rm = TRUE,
    )
  }else if (method == "median") {
    all.data.long<-aggregate(
      x = working.counts,
      by = list(gc=combined.variables$gc, mappability=combined.variables$mappability, length=combined.variables$length)[variables],
      FUN = median,
      na.rm = TRUE,
    )
  }

  rownames(all.data.long) <- as.character(interaction(all.data.long[,variables, drop = FALSE],sep="-"))

  # calculateFits <-  function(i, ...) {


  # Calculate first loess estimate
  print("    Calculating fit...")
  print(paste0("        Calculating iteration ", "1" , " out of ", maxIter, " iterations"))

  l <- loess(formula(paste("x ~", paste(variables, collapse=" * "))),
             data=all.data.long, span=span, family=family, ...)


  fit <- as.vector(predict(l, all.data.long[,variables, drop = FALSE]))
  names(fit) <- rownames(all.data.long)

  if (log.data) {
    residual <- working.counts - fit[as.character(interaction(combined.variables[,variables],sep="-"))]
  } else {
    residual <- (working.counts / fit[as.character(interaction(combined.variables[,variables],sep="-"))]) - 1

  }

  cutoffValue <- cutoff * matrixStats::madDiff(residual, na.rm=TRUE)

  prevGoodBins <- rep(TRUE, length(working.counts))
  goodBins <- !is.na(residual) &
    abs(residual) <= cutoffValue
  iter <- 1


  while(!identical(goodBins, prevGoodBins) && iter < maxIter) {
    print(paste0("        Calculating iteration ", iter + 1, " out of ", maxIter, " iterations"))
    all.data.long2<-aggregate(
      x = working.counts[goodBins],
      by = list(gc=combined.variables$gc[goodBins], mappability=combined.variables$mappability[goodBins], length=combined.variables$length[goodBins])[variables],
      FUN = mean,
      na.rm = TRUE,
    )
    rownames(all.data.long2) <- as.character(interaction(all.data.long2[,variables, drop = FALSE],sep="-"))

    l2 <- loess(formula(paste("x ~", paste(variables, collapse=" * "))),
                data=all.data.long2, span=span, family=family, ...)

    fit2 <- as.vector(predict(l2, all.data.long[,variables, drop = FALSE]))
    names(fit2) <- rownames(all.data.long)


    fit[!is.na(fit2)] <- fit2[!is.na(fit2)]

    if (log.data) {
      residual <- working.counts - fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))]
    } else {
      residual <- (working.counts / fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))]) - 1
    }



    # residual <- log2(working.counts / fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))])
    # residual <- working.counts - fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))]
    prevGoodBins <- goodBins
    goodBins <- !is.na(residual) &
      abs(residual) <= cutoffValue
    iter <- iter + 1
  }

  if (log.data) {
    mean_residual <- aggregate(
      x = working.counts - fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))],
      # x = log2(working.counts / fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))]),
      by = list(gc=combined.variables$gc, mappability=combined.variables$mappability, length=combined.variables$length)[variables],
      FUN = mean,
      na.rm = TRUE,
    )$x
  } else {
    mean_residual <- aggregate(
      # x = working.counts - fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))],
      x = working.counts / fit[as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))],
      by = list(gc=combined.variables$gc, mappability=combined.variables$mappability, length=combined.variables$length)[variables],
      FUN = mean,
      na.rm = TRUE,
    )$x
  }


  if (log.data) {
    residual <- rep(NA, times = length(counts))
    fit.data <- rep(NA, times = length(counts))
    # work in log space
    have.correction <- as.character(interaction(combined.variables.all[,variables, drop = FALSE],sep="-")) %in% as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))
    log2.values <- log2(counts[have.correction])
    log2.values[is.infinite(log2.values)] <- 0
    residual[have.correction] <-
      log2.values - fit[as.character(interaction(combined.variables.all[have.correction,variables, drop = FALSE],sep="-"))]
    fit.data[have.correction] <-
      fit[as.character(interaction(combined.variables.all[have.correction,variables, drop = FALSE],sep="-"))]
  } else {
    residual <- rep(NA, times = length(counts))
    fit.data <- rep(NA, times = length(counts))

    have.correction <- as.character(interaction(combined.variables.all[,variables, drop = FALSE],sep="-")) %in% as.character(interaction(combined.variables[,variables, drop = FALSE],sep="-"))
    residual[have.correction] <-
      counts[have.correction] / fit[as.character(interaction(combined.variables.all[have.correction,variables, drop = FALSE],sep="-"))]
    fit.data[have.correction] <-
      fit[as.character(interaction(combined.variables.all[have.correction,variables, drop = FALSE],sep="-"))]

  }

  outdata <- cbind(counts,residual, fit.data, combined.variables.all)
  colnames(outdata)[1:3] <- c("counts", "residual", "fit")
  attr(outdata, "used.span") <- span
  attr(outdata, "used.family") <- family
  attr(outdata, "correction_variables") <- paste(variables, sep = "_")
  print("Done.")
  return(outdata)
}


#' Estimate GC and read length bias in the data
#'
#' This script performs estimation of bias within the data.
#'
#' @param object Output of mutREADestimateCorrection
#' @param bins The AnnotatedDataFrame object containing information about bins to be analysed
#' @param span the size of bin used for the averaging of signal from target regions
#' @param method method used for calculation of the average signal
#' @param correction.method Correction method used for calculation of signal enrichment in the bin over the signal expected in that bin
#' @param genomebuild Genome version. Can be hg19 or hg38
#' @param keep.empty.bins Boolean whether to keep the bins without any target regions
#' @param n.cores number of cores for paralelisation
#' @author Karol Nowicki-Osuch
#' @export
mutREADcombineData <- function(object, bins, span = 50000L , method = "sum", correction.method = "log.ratio", genomebuild = "hg19", keep.empty.bins = FALSE, n.cores = 1) {
  # object <- corrected.data
  # chr.middle.hg19<-c(
  #   "chr1" = 124300000,
  #   "chr10" = 40300000,
  #   "chr11" = 52900000,
  #   "chr12" = 35400000,
  #   "chr13" = 16000000,
  #   "chr14" = 15600000,
  #   "chr15" = 17000000,
  #   "chr16" = 38200000,
  #   "chr17" = 22200000,
  #   "chr18" = 16100000,
  #   "chr19" = 28500000,
  #   "chr2"  = 93300000,
  #   "chr20" = 27100000,
  #   "chr21" = 12300000,
  #   "chr22" = 11800000,
  #   "chr3"  = 91700000,
  #   "chr4"  = 50700000,
  #   "chr5"  = 47700000,
  #   "chr6"  = 60500000,
  #   "chr7"  = 59100000,
  #   "chr8"  = 45200000,
  #   "chr9"  = 51800000,
  #   "chrX"  = 59500000,
  #   "chrY"  = 11300000
  # )
  # names(chr.middle.hg19)<-gsub("chr", "", names(chr.middle.hg19))
  # require("BSgenome.Hsapiens.UCSC.hg19")
  # lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19)[unique(bins$Chr)]
  # names(lengths)<-gsub("chr", "", names(lengths))

  if (genomebuild == "hg19") {
    require("BSgenome.Hsapiens.UCSC.hg19")
    lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19)[unique(bins$Chr)]
    names(lengths)<-gsub("chr", "", names(lengths))
  } else if (genomebuild == "hg38") {
    require("BSgenome.Hsapiens.UCSC.hg38")
    lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38)[unique(bins$Chr)]
    # names(lengths)<-gsub("chr", "", names(lengths))
  } else {
    stop("This version of genome is not supported")
  }



  # This function merges the bins data
  bin.merger <- function(chr, start, end, data, span, method, correction.method, keep.empty.bins) {
    # chr = 1
    #Check if vector of value were provided for start and end for each chr
    if(length(start) == 1){
      start <- start
    } else {
      start <- start[chr]
    }

    if(length(end) == 1){
      end <- end
    } else {
      end <- end[chr]
    }

    # Get genomic locations
    locations <- as.data.frame(do.call(rbind,base::strsplit(rownames(data), split = ":|-", perl = TRUE)))
    # locations <- as.data.frame(stringr::str_split(rownames(data), ":|-", simplify = T))
    locations[,2] <- as.numeric(locations[,2])
    locations[,3] <- as.numeric(locations[,3])

    # get the data
    working.object <- as.matrix(object[locations[,1] == chr & locations[,2] >= start & locations[,3] <= end ,c("counts", "fit", "gc", "mappability", "length")])
    working.object[,2][working.object[,2] < 0 ] <- 0
    #The locations chromosome that should be merged
    locations <- locations[locations[,1] == chr & locations[,2] >= start & locations[,3] <= end ,]

    # Create empty matrix to store processed counts data
    out.counts <- matrix(data = 0, nrow = 0, ncol = ncol(working.object) + 1)
    colnames(out.counts) <- c(colnames(working.object), "n.targets")

    #This loop will perform mergin within the chromosomes
    # first, I create a set of ranger that will included all data
    for (range.starts in seq(1, locations[nrow(locations),3], by = span) ) {
      # print(range.starts)
      # range.starts = 1
      # get the genomic position 1 nt before the first nucleotide of the next bin
      end.nt <- range.starts + span - 1
      if (end.nt > locations[nrow(locations),3]) {
        end.nt <- locations[nrow(locations),3]
      }


      within.range<-which(locations[,2] >= range.starts & locations[,2] <= end.nt)

      if (length(within.range) == 0 & !keep.empty.bins){

        next

      } else if (length(within.range) == 0 & keep.empty.bins) {

        out.counts <- rbind(out.counts, c(0,0, NA, NA, NA, 0))
        rownames(out.counts)[nrow(out.counts)] <- paste0(chr, ":", range.starts, "-", end.nt )

      } else {

        if (method == "sum") {
          out.counts <- rbind(out.counts, c(matrixStats::colSums2(working.object[within.range, 1:2, drop = FALSE], na.rm = TRUE), # count and fit
                                            matrixStats::colMeans2(working.object[within.range, 3:5, drop = FALSE], na.rm = TRUE), # gc, mappability and length
                                            length(within.range) # number of objects
          )
          )
          rownames(out.counts)[nrow(out.counts)] <- paste0(chr, ":", range.starts, "-", end.nt )

          # out.rowData[nrow(out.rowData)+1,] <- c(
          # matrixStats::colMeans2(rowData.all[within.range,c("gc", "length", "mappability"), drop = FALSE], na.rm = TRUE),
          # length(within.range),
          # matrixStats::colSums2(rowData.all[within.range, fit.data, drop = FALSE], na.rm = TRUE)
          # )
          # rownames(out.rowData)[nrow(out.rowData)] <- paste0(chr, ":", range.starts, "-", end.nt )

        } else {
          stop("This method is not yet ready")
        }

      }
    }

    if( correction.method == "ratio") {
      corrected.counts <- rep(NA, times = nrow(out.counts))
      corrected.counts[out.counts[,2] != 0] <- (out.counts[,1]/out.counts[,2])[out.counts[,2] != 0]


    } else if (correction.method == "log.ratio") {
      corrected.counts <- rep(NA, times = nrow(out.counts))
      corrected.counts[out.counts[,2] != 0] <- log2((out.counts[,1]/out.counts[,2])[out.counts[,2] != 0])
      corrected.counts[out.counts[,1] == 0] <- -3

    } else {
      stop("This correction method is not yet ready")
    }
    out.counts <- cbind(out.counts, corrected.counts)
    return(out.counts)
  }

  if (n.cores == 1) {

    combined.data <-
      c(
        lapply(
          unique(names(lengths)),
          bin.merger,
          start = 1,
          end = lengths,
          data = object,
          span = span,
          method = method,
          correction.method = correction.method,
          keep.empty.bins = keep.empty.bins)
      )

  } else {

    combined.data <-
      c(
        parallel::mclapply(
          unique(names(lengths)),
          bin.merger,
          start = 1,
          end = lengths,
          data = object,
          span = span,
          method = method,
          correction.method = correction.method,
          keep.empty.bins = keep.empty.bins,
          mc.cores = n.cores)
      )
  }
  # print("here")
  combined.data <- do.call(rbind, combined.data)

  return(combined.data)
}

#' Correct the signal.
#'
#' Not in use
#'
#' @author Karol Nowicki-Osuch
#' @noRd
mutREADcorrectData <- function(object, correction.method = "log.ratio") {
  # object <- corrected.data
  object[object[,3] < 0,3] <- 0
  if( correction.method == "ratio") {
    corrected.counts <- rep(NA, times = nrow(object))
    corrected.counts[object[,3] != 0] <- (object[,1]/object[,3])[object[,3] != 0]


  } else if (correction.method == "log.ratio") {
    corrected.counts <- rep(NA, times = nrow(object))
    corrected.counts[object[,3] != 0] <- log2((object[,1]/object[,3])[object[,3] != 0])
    corrected.counts[object[,1] == 0] <- -3

  } else {
    stop("This correction method is not yet ready")
  }
  out.data <- cbind(object[,c(1,3:6)], 1 ,corrected.counts)
  colnames(out.data)[(ncol(out.data)-1):ncol(out.data)] <- c("n.targets", "corrected.counts")
  # object <- cbind(object, corrected.counts)

  # colnames(object)[ncol(object)] <- "corrected"
  return(out.data)

}


#' Calculater Running median
#'
#' @author Karol Nowicki-Osuch
#' @noRd
runmed_data = function(chromosome, data, k=101) {
  data_smoothed = rep(NA, length(data))
  for (chrom in unique(chromosome)) {
    data_smoothed[chromosome==chrom] = runmed(data[chromosome==chrom], k, endrule = "constant")
  }
  return(data_smoothed)
}

#' Plotting function
#'
#' @author Karol Nowicki-Osuch
#' @noRd
coverage.plot.mutREAD = function(object, what="corrected.counts", min.y=-2, max.y=2, pct.plot = 1, horiz.line.dist = 0.5, run.median.k = 101) {
  # object <- combined.data
  require("viridis")
  require("cowplot")
  locations <- as.data.frame(do.call(rbind,base::strsplit(rownames(object), split = ":|-", perl = TRUE)))[,1:2]
  # locations <- as.data.frame(stringr::str_split(rownames(object), ":|-", simplify = T)[,1:2])
  colnames(locations) <- c("Chromosome", "Start")

  locations$Chromosome <- factor(locations$Chromosome, levels=gtools::mixedsort(unique(locations$Chromosome)))
  locations$Start <- as.numeric(locations$Start)
  locations$data_binned = runmed_data(locations$Chromosome, object[,what], k=run.median.k)


  background = data.frame(y=seq(round(min.y),round(max.y), by = horiz.line.dist))
  # plot_title = samplename
  p = ggplot(locations[seq(1, nrow(locations),round(100/pct.plot)),]) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=Start, y=data_binned), alpha=0.5, size=0.5, colour="darkgreen") +
    facet_grid(~Chromosome, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    # geom_smooth() +
    scale_y_continuous(breaks = seq(min.y, max.y, by = 1), limits = c(min.y, max.y)) +
    # ylim(min.y,max.y) +
    #ylab() +
    # ggtitle(plot_title) +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))

  p
}

#' Plotting function 2
#'
#' @author Karol Nowicki-Osuch
#' @noRd
coverage.plot.mutREAD_2 = function(object, min.y=-2, max.y=2, pct.plot = 1, horiz.line.dist = 0.5, run.median.k = 101) {
  require("viridis")
  require("cowplot")
  locations <- as.data.frame(do.call(rbind,base::strsplit(rownames(object), split = ":|-", perl = TRUE)))[,1:2]
  # locations <- as.data.frame(stringr::str_split(rownames(object), ":|-", simplify = T)[,1:2])
  colnames(locations) <- c("Chromosome", "Start")

  locations$Chromosome <- factor(locations$Chromosome, levels=gtools::mixedsort(unique(locations$Chromosome)))
  locations$Start <- as.numeric(locations$Start)
  locations$data_binned = runmed_data(locations$Chromosome, object[,"corrected.counts"], k=run.median.k)


  locations$segments <- object[,"state.median"]
  locations$segments[2:length(locations$segments)] <- ifelse(locations$segments[2:length(locations$segments)] == locations$segments[1:(length(locations$segments)-1)], locations$segments[2:length(locations$segments)], NA)

  states <- c("Deep Loss", "Loss", "Neutral", "Gain", "Amplification", "High Amp")
  locations$state <- states[object[,"state"]]
  locations$state[is.na(locations$segments)] <- NA
  colors.range = c('#9bf542', 'green', 'blue', 'red',
                   'darkred', 'orange')
  names(colors.range) <- states


  background = data.frame(y=seq(round(min.y),round(max.y), by = horiz.line.dist))

  # plot_title = samplename
  p = ggplot(locations[seq(1, nrow(locations),round(100/pct.plot)),]) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=Start, y=data_binned, color = state), alpha=0.5, size=0.5) + #, colour="darkgreen") +
    geom_line(mapping=aes(x=Start, y=segments), size = 2, color = "black") +
    facet_grid(~Chromosome, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    # geom_smooth() +
    scale_y_continuous(breaks = seq(min.y, max.y, by = 1), limits = c(min.y, max.y)) +
    scale_color_manual(values = colors.range) +
    # ylim(min.y,max.y) +
    #ylab() +
    # ggtitle(plot_title) +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))



  p
}

#' Combined breakpoint file
#'
#' @author Karol Nowicki-Osuch
#' @noRd
combine.breakpoints = function(x, out.name) {
  all.breakpoints <- list()
  for(n in x) {
    all.breakpoints[[n]] <- read.delim2(n)
  }
  colnames(all.breakpoints[[2]]) <- colnames(all.breakpoints[[1]])

  all.breakpoints <- do.call(rbind, all.breakpoints)
  all.breakpoints <- all.breakpoints[gtools::mixedorder(paste(all.breakpoints[,1], all.breakpoints[,2], sep = "_")),]
  all.breakpoints <- all.breakpoints[!duplicated(all.breakpoints),]
  write.table(all.breakpoints, paste0(directory, "/", out.name), quote = FALSE, sep = "\t", row.names = FALSE)
}
