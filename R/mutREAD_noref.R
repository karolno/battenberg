#' Process single mutREAD file
#'
#' @author Karol Nowicki-Osuch
#' @noRd
process.mutREAD.single <- function(binspan = 500000L, tumourbam,
                                   tumour.sample, directory, bins, genomebuild = "hg19", nthreads = 1,
                                   segment = FALSE)
{
  bins <- readRDS(file = bins)
  if (!requireNamespace("HMMcopy", quietly = TRUE)) {
    stop("Package \"HMMcopy\" must be installed to use this function.",
         call. = FALSE)
  }
  require("GenomicRanges")
  require("IRanges")
  require("Biobase")
  if (genomebuild == "hg19") {
    require("BSgenome.Hsapiens.UCSC.hg19")
  }
  else if (genomebuild == "hg38") {
    require("BSgenome.Hsapiens.UCSC.hg38")
  }
  else {
    stop("This version of genome is not supported")
  }
  copyNumbersSmooth.cancer <- process.mutREAD.bams(sample.name = tumour.sample,
                                                   bamlocation = tumourbam, bins = bins, directory = directory,
                                                   genomebuild = genomebuild, binspan = binspan, segment = segment,
                                                   nthreads = nthreads)
  # copyNumbersSmooth.normal <- process.mutREAD.bams(sample.name = ref.sample,
  #                                                  bamlocation = normalbam, bins = bins, directory = directory,
  #                                                  genomebuild = genomebuild, binspan = binspan, segment = FALSE,
  #                                                  nthreads = nthreads)
  copyNumbersSmooth <- copyNumbersSmooth.cancer
  # copyNumbersSmooth[, "corrected.counts"] <- copyNumbersSmooth[,
  #                                                              "corrected.counts"] - copyNumbersSmooth.normal[, "corrected.counts"]
  data.df <- as.data.frame(do.call(rbind, base::strsplit(rownames(copyNumbersSmooth),
                                                         split = ":|-", perl = TRUE)))
  data.df$V1 <- factor(data.df$V1, levels = gtools::mixedsort(unique(data.df$V1)))
  data.df <- cbind(data.df, copyNumbersSmooth[, "corrected.counts"])
  colnames(data.df) <- c("chr", "start", "end", "copy")
  seg.data <- HMMcopy::HMMsegment(data.df, verbose = FALSE)
  rleseg <- rle(paste0(data.df$chr, ":", seg.data$state))
  seg.medians <- rep(seg.data$segs$median, times = rleseg$lengths)
  breakpoints.out <- seg.data$segs[duplicated(seg.data$segs$chr),
                                   1:2]
  colnames(breakpoints.out) <- c("chromosome", "position")
  write.table(breakpoints.out, paste0(directory, "/", tumour.sample,
                                      "_breakpoints.after.normal.tab"), quote = FALSE, sep = "\t",
              row.names = FALSE)
  if (segment) {
    copyNumbersSmooth[, "state"] <- seg.data$state
    copyNumbersSmooth[, "state.median"] <- seg.medians
  }
  else {
    copyNumbersSmooth <- cbind(copyNumbersSmooth, state = seg.data$state,
                               state.median = seg.medians)
  }
  coverage.plot.state.after.normal <- coverage.plot.mutREAD_2(copyNumbersSmooth,
                                                                           pct.plot = 100, run.median.k = 1)
  ggsave(paste0(directory, "/", tumour.sample, "_mutREAD.combined.coverage.segmented.after.normal.plot.png"),
         plot = coverage.plot.state.after.normal, device = "png",
         height = 6, width = 15)
  if (genomebuild == "hg19") {
    bins.intervals <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19),
                                 tilewidth = binspan, cut.last.tile.in.chrom = TRUE)
    seqlevelsStyle(bins.intervals) <- "NCBI"
    bins.intervals <- keepSeqlevels(bins.intervals, c(1:22,
                                                      "X", "Y"), pruning.mode = "tidy")
    names(bins.intervals) <- paste0(seqnames(bins.intervals),
                                    ":", ranges(bins.intervals))
  }
  else {
    bins.intervals <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg38),
                                 tilewidth = binspan, cut.last.tile.in.chrom = TRUE)
    bins.intervals <- keepSeqlevels(bins.intervals, paste0("chr",
                                                           c(1:22, "X", "Y")), pruning.mode = "tidy")
    names(bins.intervals) <- paste0(seqnames(bins.intervals),
                                    ":", ranges(bins.intervals))
  }
  final.regions.g <- bins.intervals
  out.data <- GenomicRanges::makeGRangesFromDataFrame(data.df,
                                                      keep.extra.columns = TRUE)
  hits <- GenomicRanges::findOverlaps(final.regions.g, out.data,
                                      ignore.strand = TRUE, select = "all", minoverlap = 1L)
  hitsByQuery <- as(hits, "List")
  final.regions.g$copynumber <- IRanges::median(IRanges::extractList(out.data$copy,
                                                                     hitsByQuery), na.rm = TRUE)
  original.logR <- read.delim(paste0(directory, "/", tumour.sample,
                                     "_mutantLogR.tab"))
  original.logR <- cbind(original.logR[, 1:2], original.logR[,
                                                             2:3])
  colnames(original.logR) <- c("chromosome", "start", "end",
                               colnames(original.logR)[4])
  original.logR.g <- GenomicRanges::makeGRangesFromDataFrame(original.logR,
                                                             keep.extra.columns = TRUE)
  hits2 <- GenomicRanges::findOverlaps(original.logR.g, final.regions.g,
                                       ignore.strand = TRUE, select = "all", minoverlap = 1L)
  hitsByQuery2 <- as(hits2, "List")
  elementMetadata(original.logR.g)[, 1] <- IRanges::median(IRanges::extractList(final.regions.g$copynumber,
                                                                                hitsByQuery2), na.rm = TRUE)
  for (chromosome in seqlevels(original.logR.g)) {
    tmp.data <- elementMetadata(original.logR.g)[original.logR.g@seqnames ==
                                                   chromosome, 1]
    tmp.which <- which(is.na(tmp.data))
    tmp.which2 <- which(!is.na(tmp.data))
    for (n in tmp.which) {
      tmp.data[n] <- median(c(tmp.data[rev(tmp.which2[tmp.which2 <
                                                        n])[1]], tmp.data[tmp.which2[tmp.which2 > n][1]]),
                            na.rm = TRUE)
    }
    elementMetadata(original.logR.g)[original.logR.g@seqnames ==
                                       chromosome, 1] <- tmp.data
  }
  original.logR[, 4] <- elementMetadata(original.logR.g)[,
                                                         1]
  original.logR <- original.logR[, c(1, 2, 4)]
  colnames(original.logR)[1:2] <- c("Chromosome", "Position")
  write.table(original.logR, paste0(directory, "/", tumour.sample,
                                    "_mutantLogR_gcCorrected.tab"), quote = FALSE, sep = "\t",
              row.names = FALSE)
}

#' Process single mutREAD file
#'
#' @author Karol Nowicki-Osuch
#' @noRd
prepare_mutREAD_noref <- function(chrom_names, sample.bam, sample.name, normal.name,
                                  g1000lociprefix, g1000allelesprefix, min_base_qual, min_map_qual,
                                  allelecounter_exe, min.het.prop = 0.9, skip_allele_counting = F, skip_process_mutREAD = F, skip_reconstruct_normal = F, binspan, nthreads = 1, bins, genomebuild)
{
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")

  clp = parallel::makeCluster(nthreads)
  doParallel::registerDoParallel(clp)

  if (!skip_allele_counting) {
    foreach::foreach(i = 1:length(chrom_names), .packages = c("Battenberg", "copynumber",
                                                              "ggplot2", "grid")) %dopar%
      {

        getAlleleCounts(bam.file = sample.bam, output.file = paste(sample.name,
                                                                   "_alleleFrequencies_chr", chrom_names[i], ".txt", sep = ""),
                        g1000.loci = paste(g1000lociprefix, chrom_names[i], ".txt",
                                           sep = ""), min.base.qual = min_base_qual,
                        min.map.qual = min_map_qual, allelecounter.exe = allelecounter_exe)
      }


  }
  mutREAD_noref_baf_logR(sample.name, g1000allelesprefix, chrom_names)

  # Fix chromosome names
  BAF.data <- read.delim(paste0(sample.name, "_mutantBAF.tab"))
  if (!grepl("chr", BAF.data[1,1])) {
    BAF.data[,1] <- paste0("chr",BAF.data[,1])
    write.table(BAF.data, paste0(sample.name, "_mutantBAF.tab"), col.names = T,
                row.names = F, quote = F, sep = "\t")
  }


  LogR.data <- read.delim(paste0(sample.name, "_mutantLogR.tab"))
  if (!grepl("chr", LogR.data[1,1])) {
    LogR.data[,1] <- paste0("chr",LogR.data[,1])
    write.table(LogR.data, paste0(sample.name, "_mutantLogR.tab"), col.names = T,
                row.names = F, quote = F, sep = "\t")
  }

  if (!skip_process_mutREAD) {
    process.mutREAD.single(binspan = binspan, tumourbam = sample.bam,
                           tumour.sample = sample.name,
                           directory = getwd(), nthreads = nthreads, bins = bins,
                           genomebuild = genomebuild, segment = TRUE)


  }

  if (!skip_reconstruct_normal) {
    foreach::foreach(i = 1:length(chrom_names), .packages = c("Battenberg", "copynumber",
                                                              "ggplot2", "grid")) %dopar%
      {
        multi_sample_reconstruct_normal_mutREAD(sample.name = sample.name, chrom = chrom_names[i],
                                                normal.name = normal.name)
      }
  }

  parallel::stopCluster(clp)

}

#' Reconstruct normal sample
#'
#' @author Karol Nowicki-Osuch
#' @noRd
multi_sample_reconstruct_normal_mutREAD <- function(sample.name, chrom,
                                                    normal.name) {
  sample.read <- read.delim(paste0(sample.name, "_alleleFrequencies_chr",chrom, ".txt"), header = TRUE)
  sample.read[!grepl(pattern = "chr", sample.read[,1]),1] <- paste0("chr",sample.read[!grepl(pattern = "chr", sample.read[,1]),1] )
  colnames(sample.read)[1] <- "#CHR"
  if (file.exists(paste0(normal.name, "_alleleFrequencies_chr",chrom, ".txt"))) {
    combined.sample.read <- read.delim(paste0(normal.name, "_alleleFrequencies_chr",chrom, ".txt"), header = TRUE)
    colnames(combined.sample.read)[1] <- "#CHR"
    combined.sample.read[,3:7] <- combined.sample.read[,3:7] + sample.read[,3:7]

  } else {
    combined.sample.read <- sample.read
  }
  write.table(sample.read, paste0(sample.name, "_alleleFrequencies_chr",chrom, ".txt"), col.names = T, row.names = F,
              quote = F, sep = "\t")

  write.table(combined.sample.read, paste0(normal.name, "_alleleFrequencies_chr",chrom, ".txt"), col.names = T, row.names = F,
              quote = F, sep = "\t")
}

# This one is used
run_haplotyping_noref_mutREAD <- function(chrom, samplename, normalname, ismale, imputeinfofile,
                                          problemloci, impute_exe, min_normal_depth, chrom_names, g1000allelesprefix,
                                          externalhaplotypeprefix = NA, use_previous_imputation = F,
                                          snp6_reference_info_file = NA, heterozygousFilter = NA,
                                          usebeagle = FALSE, beaglejar = NA, beagleref = NA, beagleplink = NA,
                                          beaglemaxmem = 10, beaglenthreads = 1, beaglewindow = 40,
                                          beagleoverlap = 4, javajre = "java", iter = 12, burn.in = 3, imp.segment = 6, phase.states = 280)
{
  previoushaplotypefile <- list.files(pattern = paste0("_impute_output_chr",
                                                       chrom, "_allHaplotypeInfo.txt"))[1]
  if (use_previous_imputation & !is.na(previoushaplotypefile)) {
    print(paste0("Previous imputation results found, copying info from",
                 previoushaplotypefile, " to flip alleles"))
    currenthaplotypefile <- paste(samplename, "_impute_output_chr",
                                  chrom, "_allHaplotypeInfo.txt", sep = "")
    if (previoushaplotypefile != currenthaplotypefile) {
      file.copy(from = previoushaplotypefile, to = paste(samplename,
                                                         "_impute_output_chr", chrom, "_allHaplotypeInfo.txt",
                                                         sep = ""))
    }
  }
  else {
    if (file.exists(paste(samplename, "_alleleFrequencies_chr",
                          chrom, ".txt", sep = ""))) {
      generate.impute.input.mutREAD.noref(chrom = chrom,
                                             sample.allele.counts.file = paste(samplename,
                                                                                 "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
                                             normal.allele.counts.file = paste(normalname,
                                                                               "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
                                             output.file = paste(samplename, "_impute_input_chr",
                                                                 chrom, ".txt", sep = ""),
                                          g1000allelesprefix = g1000allelesprefix,
                                             is.male = ismale, problemLociFile = problemloci, heterozygousFilter = heterozygousFilter, min.normal.depth = min_normal_depth,
                                             useLociFile = NA)
    }
    else {
      stop("Germline calling is currently on WGS data only - SNP array data is not sufficiently dense to detect all germline CNVs")
    }
    if (usebeagle) {
      imputeinputfile <- paste(samplename, "_impute_input_chr",
                               chrom, ".txt", sep = "")
      vcfbeagle <- convert.impute.input.to.beagle.input(imputeinput = imputeinputfile,
                                                        chrom = ifelse(grepl("chr", chrom), chrom, paste0("chr", chrom)))
      vcfbeagle_path <- paste(samplename, "_beagle5_input_chr",
                              chrom, ".txt", sep = "")
      outbeagle_path <- paste(samplename, "_beagle5_output_chr",
                              chrom, ".txt", sep = "")
      writevcf.beagle(vcfbeagle, filepath = vcfbeagle_path)
      run.beagle5_mutREAD(beaglejar = beaglejar, vcfpath = vcfbeagle_path,
                          reffile = beagleref, outpath = outbeagle_path,
                          plinkfile = beagleplink, maxheap.gb = beaglemaxmem,
                          nthreads = beaglenthreads, window = beaglewindow,
                          overlap = beagleoverlap, javajre = javajre, iter = iter, burn.in = burn.in, imp.segment = imp.segment, phase.states = phase.states)
      outfile <- paste(samplename, "_impute_output_chr",
                       chrom, "_allHaplotypeInfo.txt", sep = "")
      vcfout <- paste(outbeagle_path, ".vcf.gz", sep = "")
      writebeagle.as.impute(vcf = vcfout, outfile = outfile)
    }
    else {
      run.impute(inputfile = paste(samplename, "_impute_input_chr",
                                   chrom, ".txt", sep = ""), outputfile.prefix = paste(samplename,
                                                                                       "_impute_output_chr", chrom, ".txt", sep = ""),
                 is.male = ismale, imputeinfofile = imputeinfofile,
                 impute.exe = impute_exe, region.size = 5000000,
                 chrom = chrom)
      combine.impute.output(inputfile.prefix = paste(samplename,
                                                     "_impute_output_chr", chrom, ".txt", sep = ""),
                            outputfile = paste(samplename, "_impute_output_chr",
                                               chrom, "_allHaplotypeInfo.txt", sep = ""),
                            is.male = ismale, imputeinfofile = imputeinfofile,
                            region.size = 5000000, chrom = chrom)
      unlink(paste(samplename, "_impute_output_chr",
                   chrom, ".txt*K.txt*", sep = ""))
    }
  }
  allelefrequenciesfile <- paste0(samplename, "_alleleFrequencies_chr",
                                  chrom, ".txt")
  print(allelefrequenciesfile)
  print(file.exists(allelefrequenciesfile))
  if (file.exists(allelefrequenciesfile)) {
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix,
                                                              chrom, ".vcf"))) {
      print("Adding in the external haplotype blocks")
      GetChromosomeBAFs(chrom = chrom, SNP_file = allelefrequenciesfile,
                        haplotypeFile = paste(samplename, "_impute_output_chr",
                                              chrom, "_allHaplotypeInfo.txt", sep = ""),
                        samplename = samplename, outfile = paste(samplename,
                                                                   "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt",
                                                                   sep = ""), chr_names = chrom_names, minCounts = min_normal_depth)
      plot.haplotype.data(haplotyped.baf.file = paste(samplename,
                                                                   "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt",
                                                                   sep = ""), imageFileName = paste(samplename,
                                                                                                    "_chr", chrom, "_heterozygousData_noExt.png",
                                                                                                    sep = ""), samplename = samplename, chrom = chrom,
                                       chr_names = chrom_names)
      input_known_haplotypes(chrom = chrom, chrom_names = chrom_names,
                             imputedHaplotypeFile = paste0(samplename,
                                                           "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                             externalHaplotypeFile = paste0(externalhaplotypeprefix,
                                                            chrom, ".vcf"))
    }
    GetChromosomeBAFs(chrom = chrom, SNP_file = paste(samplename,
                                                      "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
                      haplotypeFile = paste(samplename, "_impute_output_chr",
                                            chrom, "_allHaplotypeInfo.txt", sep = ""), samplename = samplename,
                      outfile = paste(samplename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt",
                                      sep = ""), chr_names = chrom_names, minCounts = 0)
  }
  else {
    stop("Germline calling is only on WGS data - SNParray data not sufficiently dense")
  }
  plot.haplotype.data(haplotyped.baf.file = paste(samplename,
                                                               "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt",
                                                               sep = ""), imageFileName = paste(samplename, "_chr",
                                                                                                chrom, "_heterozygousData.png", sep = ""), samplename = samplename,
                                   chrom = chrom, chr_names = chrom_names)
}


#' Run beagle with additional settings
#'
#' @author Karol Nowicki-Osuch
#' @noRd
run.beagle5_mutREAD <- function (beaglejar, vcfpath, reffile, outpath, plinkfile, nthreads = 1,
                                 window = 40, overlap = 4, maxheap.gb = 10, javajre = "java", iter = 12, burn.in = 3, imp.segment = 6, phase.states = 280)
{
  cmd <- paste0(javajre, " -Xmx", maxheap.gb, "g", " -Xms",
                maxheap.gb, "g", " -XX:+UseParallelGC", " -jar ", beaglejar,
                " gt=", vcfpath, " ref=", reffile, " out=", outpath,
                " map=", plinkfile, " nthreads=", nthreads, " window=",
                window, " overlap=", overlap, " burnin=", burn.in ," iterations=", iter ," imp-segment=", imp.segment ," phase-states=", phase.states, " impute=false")
  EXIT_CODE = system(cmd, wait = T)
  stopifnot(EXIT_CODE == 0)
}

#' Prepare impute file
#'
#' @author Karol Nowicki-Osuch
#' @noRd
generate.impute.input.mutREAD.noref <- function (chrom, sample.allele.counts.file, normal.allele.counts.file,
                                                    output.file, is.male, problemLociFile = NA, g1000allelesprefix,
                                                    useLociFile = NA, heterozygousFilter = 0.1, min.normal.depth)
{
  known_SNPs <- read.delim(paste0(g1000allelesprefix, chrom, ".txt"))
  # Allele counting only takes place over known SNPs so there is no need to filter for "knownSNPs"
  normal_snp_data = read.table(normal.allele.counts.file,comment.char = "", check.names = FALSE,
                               sep = "\t", header = T, stringsAsFactors = F)
  # a copy of this file
  # write.table(normal_snp_data, paste0(normal.allele.counts.file, ".all.txt"),
  #                              sep = "\t", quote = FALSE)
  # normal_snp_data <- normal_snp_data[normal_snp_data$Good_depth > min.normal.depth,]


  if ((problemLociFile != "NA") & (!is.na(problemLociFile))) {
    problemSNPs = read.table(problemLociFile, header = T,
                             sep = "\t", stringsAsFactors = F)
    problemSNPs = problemSNPs$Pos[problemSNPs$Chr == chrom]
    badIndices = match(known_SNPs$position, problemSNPs)
    known_SNPs = known_SNPs[is.na(badIndices), ]
    rm(problemSNPs, badIndices)
  }
  if ((useLociFile != "NA") & (!is.na(useLociFile))) {
    goodSNPs = read.table(useLociFile, header = T, sep = "\t",
                          stringsAsFactors = F)
    goodSNPs = goodSNPs$pos[goodSNPs$chr == chrom]
    goodIndices = match(known_SNPs$position, goodSNPs)
    known_SNPs = known_SNPs[!is.na(goodIndices), ]
    rm(goodSNPs, goodIndices)
  }
  indices = match(known_SNPs$position, normal_snp_data[, 2])
  known_SNPs <-  known_SNPs[!is.na(indices),]

  snp_data = read.table(sample.allele.counts.file, comment.char = "", check.names = FALSE,
                        sep = "\t", header = T, stringsAsFactors = F)
  # write.table(snp_data, paste0(sample.allele.counts.file, ".all.txt"),
  #             sep = "\t", quote = FALSE)
  # snp_data <- snp_data[snp_data$POS %in% normal_snp_data$POS,]

  nucleotides = c("A", "C", "G", "T")
  ref_indices = known_SNPs[, 2] + 2
  alt_indices = known_SNPs[, 3] + 2
  alt.count <- as.numeric(normal_snp_data[cbind(indices,
                                                alt_indices)])
  ref.alt.cov <- (as.numeric(normal_snp_data[cbind(indices,
                                    alt_indices)]) + as.numeric(normal_snp_data[cbind(indices,
                                                                                      ref_indices)]))
  BAFs = alt.count/ref.alt.cov
  # write.table(normal_snp_data, normal.allele.counts.file,
  #             sep = "\t", quote = FALSE)
  # write.table(snp_data, sample.allele.counts.file,
  #             sep = "\t", quote = FALSE)
  BAFs[is.nan(BAFs)] = 0


  rm(ref_indices, alt_indices,
     normal_snp_data)
  minBaf = min(heterozygousFilter, 1 - heterozygousFilter)
  maxBaf = max(heterozygousFilter, 1 - heterozygousFilter)
  genotypes = array(0, c(sum(!is.na(indices)), 3))
  genotypes[BAFs <= minBaf, 1] = 1
  genotypes[BAFs > minBaf & BAFs < maxBaf, 2] = 1
  genotypes[BAFs >= maxBaf, 3] = 1
  snp.names = paste("snp", 1:sum(!is.na(indices)), sep = "")
  out.data = cbind(snp.names,
                   paste(chrom, known_SNPs[,1], sep = ":"),
                   known_SNPs[,1],
                   nucleotides[known_SNPs[, 2]],
                   nucleotides[known_SNPs[, 3]],
                   genotypes)[ref.alt.cov > min.normal.depth, ]
  head(out.data)
  write.table(out.data, file = output.file, row.names = F,
              col.names = F, quote = F)
  print(chrom)
  if (is.na(as.numeric(gsub("chr", "",chrom)))) {
    sample.g.file = paste(dirname(output.file), "/sample_g.txt",
                          sep = "")
    sample_g_data = data.frame(ID_1 = c(0, "INDIVI1"), ID_2 = c(0,
                                                                "INDIVI1"), missing = c(0, 0), sex = c("D", 2))
    write.table(sample_g_data, file = sample.g.file, row.names = F,
                col.names = T, quote = F)
  }
}


#' Prepare impute file
#'
#' @author Karol Nowicki-Osuch
#' @noRd
mutREAD_noref_baf_logR <- function (sample.name, g1000alleles.prefix, chrom_names)
{
  AC = list()
  AL = list()
  MaC = list()
  OHET = list()
  for (chr in chrom_names) {

    ac = read.table(paste0(sample.name, "_alleleFrequencies_chr",
                           chr, ".txt"), stringsAsFactors = F)
    ac = ac[order(ac$V2), ]
    AC[[chr]] = ac
    print(length(AC))
    al = read.table(paste0(g1000alleles.prefix, chr, ".txt"),
                    header = T, stringsAsFactors = F)
    AL[[chr]] = al
    print(length(AL))
    ref = al$a0
    ref_df = data.frame(pos = 1:nrow(al), ref = ref + 2)
    REF = ac[cbind(ref_df$pos, ref_df$ref)]
    alt = al$a1
    alt_df = data.frame(pos = 1:nrow(al), alt = alt + 2)
    ALT = ac[cbind(alt_df$pos, alt_df$alt)]
    mac = data.frame(ref = REF, alt = ALT)
    mac$depth = as.numeric(mac$ref) + as.numeric(mac$alt)
    mac$baf = as.numeric(mac$alt)/as.numeric(mac$depth)
    o = cbind(al, mac)
    names(o) = c("Position", "a0", "a1", "ref", "alt", "depth",
                 "baf")
    MaC[[chr]] = o
    # ohet = o[which(o$baf >= 0.1 & o$baf <= 0.9 & o$depth >
    #                  10), ]
    # ohet$Position2 = c(ohet$Position[2:nrow(ohet)], 2 *
    #                      ohet$Position[nrow(ohet)] - ohet$Position[nrow(ohet) -
    #                                                                  1])
    # ohet$Position_dist = ohet$Position2 - ohet$Position
    # ohet$Position_dist_percent = ohet$Position_dist/max(ohet$Position_dist)
    # OHET[[chr]] = ohet
    print(paste("chromosome", chr, "file read"))
  }
  MAC = data.frame()
  for (chr in chrom_names) {
    MaC_CHR = data.frame(chr = chr, MaC[[chr]])
    MAC = rbind(MAC, MaC_CHR)
    print(chr)
  }
  names(MAC) = c("chr", "position", "a0", "a1", "ref", "alt",
                 "coverage", "baf")
  print(head(MAC))
  print(dim(MAC))
  MAC$logr = log2(MAC$coverage/mean(MAC$coverage, na.rm = TRUE))
  MACC = MAC[which(!is.na(MAC$baf)), ]
  print(nrow(MAC) - nrow(MACC))
  BAF = data.frame(Chromosome = MACC$chr, Position = MACC$pos,
                   germline = MACC$baf)
  names(BAF)[names(BAF) == "germline"] <- sample.name
  BAF = BAF[order(BAF$Chromosome, BAF$Position), ]

  write.table(BAF, paste0(sample.name, "_mutantBAF.tab"), col.names = T,
              row.names = F, quote = F, sep = "\t")
  rm(BAF)
  LogR = data.frame(Chromosome = MACC$chr, Position = MACC$pos,
                    germline = MACC$logr)
  names(LogR)[names(LogR) == "germline"] <- sample.name
  LogR = LogR[order(LogR$Chromosome, LogR$Position), ]

  write.table(LogR, paste0(sample.name, "_mutantLogR.tab"), col.names = T,
              row.names = F, quote = F, sep = "\t")
  rm(MAC)
  rm(MaC)
  rm(MACC)
}




