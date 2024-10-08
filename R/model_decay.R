#' filter annotation for plotting coverages
#'
#' @description Removes isoform annotations that could produce ambigious reads, such as isoforms
#' that only differ by the 5' / 3' end. This could be useful for plotting average
#' coverage plots.
#'
#' @importFrom txdbmaker makeTxDbFromGFF makeTxDbFromGRanges
#' @importFrom rtracklayer import
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges strand
#' @importFrom BiocGenerics start end
#'
#' @param annotation path to the GTF annotation file, or the parsed GenomicRanges
#' object.
#' @param keep string, one of 'tss_differ' (only keep isoforms that all differ
#' by the transcription start site position), 'tes_differ' (only keep those that
#' differ by the transcription end site position), 'both' (only keep those that
#' differ by both the start and end site), or 'single_transcripts' (only keep
#' genes that contains a sinlge transcript).
#' @return GenomicRanges of the filtered isoforms
#' @examples
#' filtered_annotation <- filter_annotation(
#'   system.file("extdata", "rps24.gtf.gz", package = 'FLAMES'), keep = 'tes_differ')
#' filtered_annotation
#'
#' @md
#' @export
filter_annotation <- function(annotation, keep = "tss_differ") {
  if (is.character(annotation)) {
    annotation <- annotation |>
      txdbmaker::makeTxDbFromGFF() |>
      GenomicFeatures::transcripts()
  } else {
    annotation <- annotation |>
      txdbmaker::makeTxDbFromGRanges() |>
      GenomicFeatures::transcripts()
  }

  unique_fn <- function(x, keep) {
    if (keep == "tss_differ") {
      return(!duplicated(GenomicRanges::start(x)) & !duplicated(GenomicRanges::start(x), fromLast = TRUE))
    }
    if (keep == "tes_differ") {
      return(!duplicated(GenomicRanges::end(x)) & !duplicated(GenomicRanges::end(x), fromLast = TRUE))
    }
    if (keep == "both") {
      return(unique_fn(x, "tss_differ") & unique_fn(x, "tes_differ"))
    }
  }

  return(annotation[unique_fn(annotation, keep)])
}

#' plot read coverages
#'
#' @description Plot the average read coverages for each length bin or a
#' perticular isoform
#'
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicAlignments readGAlignments seqnames
#' @importFrom GenomicRanges width strand granges coverage
#' @importFrom Rsamtools ScanBamParam
#' @importFrom tidyr as_tibble pivot_longer
#' @importFrom dplyr filter mutate group_by summarize_at summarise across
#' @importFrom ggplot2 ggplot geom_line aes
#' @importFrom stats weighted.mean
#'
#' @param isoform string vector, provide isoform names to plot the coverage for the
#' corresponding isoforms, or provide NULL to plot average coverages for each
#' length bin
#' @param length_bins, numeric vector to specify the sizes to bin the isoforms by
#' @param bam, path to the BAM file (aligning reads to the transcriptome), or
#' the (GenomicAlignments::readGAlignments) parsed GAlignments object
#' @param weight_fn "read_counts" or "sigmoid", determins how the transcripts
#' should be weighted within length bins.
#' @return a ggplot2 object of the coverage plot(s)
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- 'https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', paste(file_url, 'fastq/sample1.fastq.gz', sep = '/')))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, 'genome.fa', paste(file_url, 'SIRV_isoforms_multi-fasta_170612a.fasta', sep = '/')))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, 'annot.gtf', paste(file_url, 'SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf', sep = '/')))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#' minimap2_realign(
#'   config = jsonlite::fromJSON(
#'     system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")),
#'   fq_in = fastq1,
#'   outdir = outdir
#' )
#' plot_coverage(bam = file.path(outdir, 'realign2transcript.bam'))
#' @md
#' @export
plot_coverage <- function(bam, isoform = NULL,
    length_bins = c(0, 1, 2, 5, 10, Inf), weight_fn = "read_counts") {

  transcript_info <- transcript_coverage(bam, isoform, length_bins, weight_fn)

  if (!is.null(isoform)) {
    p <- transcript_info |>
      tidyr::as_tibble(rownames = "transcript") |>
      tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
      dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
      ggplot2::ggplot(aes(x = x, y = coverage, color = transcript)) + geom_line()
    return(p)
  }

  mean_coverage <- transcript_info |>
    dplyr::group_by(length_bin) |>
    dplyr::summarise(dplyr::across(paste0("coverage_", 1:100), ~ stats::weighted.mean(.,
      w = weight)))

  p <- mean_coverage |>
    tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
    dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
    ggplot2::ggplot(aes(x = x, y = coverage, color = length_bin)) + geom_line()

  return(p)
}


transcript_coverage <- function(bam, isoform = NULL,
    length_bins = c(0, 1, 2, 5, 10, Inf), weight_fn = "read_counts") {

  if (!is(bam, "GAlignments")) {
    bam <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))
  }

  if (!is.null(isoform)) {
    bam <- bam[GenomicAlignments::seqnames(bam) %in% isoform]
  }

  read_counts <- table(GenomicAlignments::seqnames(bam))
  transcript_names <- names(read_counts)
  transcript_info <- data.frame(tr_length = GenomeInfoDb::seqlengths(bam)[transcript_names],
    read_counts = as.data.frame(read_counts[transcript_names])$Freq)
  transcript_info$length_bin <- cut(transcript_info$tr_length / 1000, length_bins)

  cover <- bam |>
    GenomicRanges::granges() |>
    GenomicRanges::coverage() |>
    sapply(function(x) {
      x[round(seq(1, length(x), length.out = 100), 0)] |>
        as.integer()
    }) |>
    subset(select = transcript_names) |>
    t() |>
    as.data.frame()
  colnames(cover) <- paste0("coverage_", 1:100)
  cover <- cover[transcript_names, ]

  if (weight_fn == "sigmoid") {
    weight_fn <- function(mat, read_counts) {
      sigmoid <- function(x) {
        exp(x) / (exp(x) + 1)
      }
      sigmoid((read_counts - 2000) / 500)
    }
  } else if (weight_fn == "read_counts") {
    weight_fn <- function(mat, read_counts) {read_counts}
  }

  cover <- cover / transcript_info$read_counts # scale by read counts
  transcript_info <- cbind(transcript_info, cover)
  transcript_info$weight <- weight_fn(mat, transcript_info$read_counts)
  return(transcript_info)
}
