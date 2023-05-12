#' Prepare mutREAD data for haplotype construction
#'
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, constructing BAF and logR
#' and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file
#' @param normalbam Full path to the normal BAM file
#' @param tumourname Identifier to be used for tumour output files
#' @param normalname Identifier to be used for normal output files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes alleles reference files
#' @param g1000prefix Prefix path to the 1000 Genomes SNP reference files
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allelecounter_exe Path to the allele counter executable (can be found in $PATH)
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param nthreads The number of paralel processes to run
#' @param skip_allele_counting Flag, set to TRUE if allele counting is already complete (files are expected in the working directory on disk)
#' @param skip_allele_counting_normal Flag, set to TRUE from the second sample onwards for multisample case (Default: FALSE)
#' @param binspan The span of the genomic region used to perform averaging of the signal for LogR calculation. It is an integer number. (Default: 5e5L)
#' @param bins The location of .rds file contain information about bins to be analysed
#' @param genomebuild Genome version. Can be hg19 or hg38
#' @author Karol Nowicki-Osuch
#' @export
prepare_mutREAD = function(chrom_names, tumourbam, normalbam, tumourname, normalname, g1000allelesprefix, g1000prefix, gccorrectprefix,
                           repliccorrectprefix, min_base_qual, min_map_qual, allelecounter_exe, min_normal_depth, nthreads, skip_allele_counting, skip_allele_counting_normal = F,
                           genomebuild = "hg19",
                           binspan = 5e5,
                           bins = NA

) {

  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")

  if (!skip_allele_counting) {

    clp = parallel::makeCluster(nthreads)
    doParallel::registerDoParallel(clp)

    # Obtain allele counts for 1000 Genomes locations for both tumour and normal
    foreach::foreach(i=1:length(chrom_names)) %dopar% {
      getAlleleCounts(bam.file=tumourbam,
                      output.file=paste(tumourname,"_alleleFrequencies_chr", chrom_names[i], ".txt", sep=""),
                      g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                      min.base.qual=min_base_qual,
                      min.map.qual=min_map_qual,
                      allelecounter.exe=allelecounter_exe)

      if (!skip_allele_counting_normal) {
        getAlleleCounts(bam.file=normalbam,
                        output.file=paste(normalname,"_alleleFrequencies_chr", chrom_names[i], ".txt",  sep=""),
                        g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                        min.base.qual=min_base_qual,
                        min.map.qual=min_map_qual,
                        allelecounter.exe=allelecounter_exe)
      }
    }
    # Kill the threads
    parallel::stopCluster(clp)
  }

  # Obtain BAF and LogR from the raw allele counts
  getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste(tumourname,"_alleleFrequencies_chr", sep=""),
                  normalAlleleCountsFile.prefix=paste(normalname,"_alleleFrequencies_chr", sep=""),
                  figuresFile.prefix=paste(tumourname, "_", sep=''),
                  BAFnormalFile=paste(tumourname,"_normalBAF.tab", sep=""),
                  BAFmutantFile=paste(tumourname,"_mutantBAF.tab", sep=""),
                  logRnormalFile=paste(tumourname,"_normalLogR.tab", sep=""),
                  logRmutantFile=paste(tumourname,"_mutantLogR.tab", sep=""),
                  combinedAlleleCountsFile=paste(tumourname,"_alleleCounts.tab", sep=""),
                  chr_names=chrom_names,
                  g1000file.prefix=g1000allelesprefix,
                  minCounts=min_normal_depth,
                  samplename=tumourname,
                  strip.chr=FALSE)

  # Perform GC and fragment length correction using mutREAD pipeline
  process.mutREAD(binspan=binspan,
                  tumourbam = tumourbam,
                  normalbam = normalbam,
                  ref.sample = normalname,
                  tumour.sample = tumourname,
                  directory = getwd(),
                  nthreads = nthreads,
                  bins = bins,
                  genomebuild = genomebuild)

  print(paste0("Successfully completed processing of mutREAD data for sample: ", tumourname,
               " Reference samples was: ", normalname))

}
