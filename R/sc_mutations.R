#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches
#' @importFrom dplyr group_by mutate ungroup
variant_count_tb <- function(bam_path, seqname, pos, indel, barcodes, verbose = TRUE) {
  # allele by barcode matrix (value: read count)
  tryCatch(
    {
      variant_count_matrix( # throws Rcpp::exception when no reads at pos
        bam_path = bam_path,
        seqname = seqname, pos = pos, indel = indel, barcodes = barcodes, verbose = verbose
      ) |>
        tibble::as_tibble(rownames = "allele") |>
        # pivot to long format: allele, barcode, allele_count
        tidyr::pivot_longer(
          cols = -tidyselect::matches("allele"),
          values_to = "allele_count", names_to = "barcode"
        ) |>
        dplyr::group_by(barcode) |>
        dplyr::mutate(
          cell_total_reads = sum(allele_count),
          pct = allele_count / cell_total_reads,
          pos = pos, seqname = seqname
        ) |>
        dplyr::ungroup()
    },
    error = function(e) {
      if (inherits(e, "Rcpp::exception") & conditionMessage(e) == "Failed to fetch an alignment") {
        message(paste0("No reads found at ", seqname, ":", pos, " in ", bam_path))
        message("Returning empty tibble")
        return(tibble::tibble())
      } else {
        stop(e)
      }
    }
  )
}

#' Variant count for single-cell data
#'
#' Count the number of reads supporting each variants at the given positions for each cell.
#'
#' @importFrom BiocParallel bplapply bpmapply MulticoreParam
#' @importFrom dplyr mutate select bind_rows
#'
#' @param bam_path character(1) or character(n): path to the bam file(s) aligned to the
#' reference genome (NOT the transcriptome! Unless the postions are also from the transcriptome).
#' @param seqnames character(n): chromosome names of the postions to count alleles.
#' @param positions integer(n): positions, 1-based, same length as seqnames. The positions to count alleles.
#' @param indel logical(1): whether to count indels (TRUE) or SNPs (FALSE).
#' @param barcodes character(n) when bam_path is a single file, or list of character(n)
#' when bam_path is a list of files paths. The cell barcodes to count alleles for.
#' Only reads with these barcodes will be counted.
#' @param threads integer(1): number of threads to use. Maximum number of threads is
#' the number of bam files * number of positions.
#' @return A tibble with columns: allele, barcode, allele_count, cell_total_reads, pct, pos, seqname.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#' minimap2_align( # align to genome
#'   config = jsonlite::fromJSON(
#'     system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")
#'   ),
#'   fa_file = genome_fa,
#'   fq_in = system.file("extdata", "fastq", "demultiplexed.fq.gz", package = "FLAMES"),
#'   annot = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   outdir = outdir
#' )
#' snps_tb <- sc_mutations(
#'   bam_path = file.path(outdir, "align2genome.bam"),
#'   seqnames = c("chr14", "chr14"),
#'   positions = c(1260, 2714), # positions of interest
#'   indel = FALSE,
#'   barcodes = read.delim(
#'     system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'     header = FALSE)$V1
#' )
#' head(snps_tb)
#' snps_tb |>
#'   dplyr::filter(pos == 1260) |>
#'   dplyr::group_by(allele) |>
#'   dplyr::summarise(count = sum(allele_count)) # should be identical to samtools pileup
#' @export
sc_mutations <- function(bam_path, seqnames, positions, indel = FALSE, barcodes, threads = 1) {
  stopifnot(
    "seqnames not the same length as positions" =
      length(seqnames) == length(positions)
  )

  if (length(bam_path) == 1) {
    # single bam file, parallelize over positions
    stopifnot(
      "barcodes must be a character vector" =
        is.character(barcodes)
    )
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got 1 bam file, parallelizing over each position ..."))
    variants <- tryCatch({
      BiocParallel::bpmapply(
        function(seqname, pos) {
          variant_count_tb(bam_path, seqname, pos, indel, barcodes, verbose = FALSE)
        },
        seqname = seqnames, pos = positions, SIMPLIFY = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE
        )
      )
    }, error = identity)
  } else {
    # multiple bam files, parallelize over bam files
    stopifnot(
      "barcodes must be a list of character vectors, same length as bam_path" =
        is.list(barcodes) & length(barcodes) == length(bam_path)
    )
    # data frame of all combinations between (seqname, pos) and (bam_path, barcodes)
    args_grid <- expand.grid(
      mutation_index = seq_along(positions),
      bam_index = seq_along(bam_path),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        seqname = seqnames[mutation_index],
        pos = positions[mutation_index],
        sample_bam = bam_path[bam_index],
        sample_barcodes = barcodes[bam_index]
      ) |>
      dplyr::select(-mutation_index, -bam_index)

    message(paste0(format(Sys.time(), "%H:%M:%S "), "Multi-threading over bam files x positions ..."))
    variants <- tryCatch({
      BiocParallel::bpmapply(
        function(sample_bam, seqname, pos, sample_barcodes) {
          variant_count_tb(sample_bam, seqname, pos, indel, sample_barcodes, verbose = FALSE) |>
            dplyr::mutate(bam_file = sample_bam)
        },
        sample_bam = args_grid$sample_bam, seqname = args_grid$seqname,
        pos = args_grid$pos, sample_barcodes = args_grid$sample_barcodes,
        SIMPLIFY = FALSE, BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE
        )
      )
    }, error = identity)
  }

  if (inherits(variants, "error")) {
    warning("Error occurred in `sc_mutations`, returning error object")
    return(variants)
  }

  message(paste0(format(Sys.time(), "%H:%M:%S "), "Merging results ..."))
  variants <- variants |>
    dplyr::bind_rows()
  return(variants)
}


extract_nt <- function(ref, seqname, pos) {
  mapply(function(seqname, pos) {
    as.character(ref[[seqname]][pos])
  }, seqname, pos)
}

#' @importFrom BiocParallel bplapply bpmapply MulticoreParam
#' @importFrom dplyr mutate pull
homopolymer_pct <- function(ref, seqname, pos, include_alt = FALSE, n = 3, threads = 1) {
  pcts <- tryCatch({
    BiocParallel::bpmapply(
      function(seqname, pos, include_alt) {
        if (pos == 1 | pos == length(ref[[seqname]])) {
          return(TRUE) # variant at the ends should not be considered
        }
        start <- max(1, pos - n)
        end <- min(length(ref[[seqname]]), pos + n)
        if (include_alt) {
          ref[[seqname]][start:end] |>
            as.character() |>
            strsplit("") |>
            unlist() |>
            table() |>
            as.data.frame() |>
            dplyr::mutate(pct = Freq / sum(Freq)) |>
            dplyr::pull(pct) |>
            max()
        } else {
          # exclude the position itself
          ref[[seqname]][c(start:(pos - 1), (pos + 1):end)] |>
            as.character() |>
            strsplit("") |>
            unlist() |>
            table() |>
            as.data.frame() |>
            dplyr::mutate(pct = Freq / sum(Freq)) |>
            dplyr::pull(pct) |>
            max()
        }
      },
      seqname, pos, include_alt, SIMPLIFY = TRUE,
      BPPARAM = BiocParallel::MulticoreParam(
        workers = threads, stop.on.error = TRUE, progressbar = TRUE
      )
    )
  }, error = identity)
  if (inherits(pcts, "error")) {
    warning("Error occurred while running `homopolymer_pct`, returning error object")
  }
  return(pcts)
}

# WIP: too much sequencing errrors / splice sites
# find variants in a single grange
#' @importFrom Rsamtools pileup PileupParam ScanBamParam
#' @importFrom dplyr select mutate group_by ungroup filter
#' @importFrom S4Vectors mcols
find_variants_grange <- function(bam_path, reference, gene_grange, min_nucleotide_depth,
                                 names_from) {
  # read bam file
  mutations <- Rsamtools::pileup(bam_path,
    pileupParam = Rsamtools::PileupParam(
      max_depth = .Machine$integer.max - 1, min_base_quality = 0, min_mapq = 0,
      min_nucleotide_depth = min_nucleotide_depth, min_minor_allele_depth = 0,
      distinguish_strands = FALSE, distinguish_nucleotides = TRUE,
      ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = TRUE,
      left_bins = NULL, query_bins = NULL, cycle_bins = NULL
    ),
    scanBamParam = Rsamtools::ScanBamParam(which = gene_grange)
  ) |>
    dplyr::select(-which_label) |>
    # move insertion postions to the base at which insertion start
    # like in IGV
    # don't count insertion in sum
    dplyr::mutate(
      pos = ifelse(nucleotide == "+", pos + 1, pos),
      counts_no_ins = ifelse(nucleotide != "+", count, 0)
    ) |>
    dplyr::group_by(seqnames, pos) |>
    dplyr::mutate(sum = sum(counts_no_ins)) |>
    dplyr::ungroup() |>
    dplyr::select(-counts_no_ins) |>
    dplyr::mutate(
      freq = count / sum,
      ref = factor(extract_nt(reference, seqnames, pos))
    ) |>
    dplyr::filter(as.character(nucleotide) != as.character(ref))

  if (nrow(mutations) == 0) {
    return(mutations)
  } else {
    mutations$bam_path <- bam_path
    if (names_from %in% colnames(S4Vectors::mcols(gene_grange))) {
      mutations$region <- S4Vectors::mcols(gene_grange)[, names_from]
    } else {
      mutations$region <- NA # no gene name / gap
    }
    return(mutations)
  }
}

#' bulk variant identification
#'
#' Treat each bam file as a bulk sample and identify variants against the reference
#'
#' Each bam file is treated as a bulk sample to perform pileup and identify variants.
#' You can run \code{sc_mutations} with the variants identified with this function
#' to get single-cell allele counts. Note that reference genome FASTA files may have
#' the chromosome names field as `>chr1 1` instead of `>chr1`. You may need to remove
#' the trailing number to match the chromosome names in the bam file, for example with
#' \code{names(ref) <- sapply(names(ref), function(x) strsplit(x, " ")[[1]][1])}.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames seqlengths seqinfo
#' @importFrom GenomicRanges gaps
#' @importFrom BiocParallel bplapply bpmapply MulticoreParam
#' @importFrom dplyr bind_rows mutate
#'
#' @param bam_path character(1) or character(n): path to the bam file(s) aligned to the
#' reference genome (NOT the transcriptome!).
#' @param reference DNAStringSet: the reference genome
#' @param annotation GRanges: the annotation of the reference genome. You can load
#' a GTF/GFF annotation file with \code{anno <- rtracklayer::import(file)}.
#' @param min_nucleotide_depth integer(1): minimum read depth for a position to be
#' considered a variant.
#' @param threads integer(1): number of threads to use. Threading is done over each
#' annotated region and (if \code{annotated_region_only = FALSE}) unannotated gaps for
#' each bam file.
#' @param homopolymer_window integer(1): the window size to calculate the homopolymer
#' percentage. The homopolymer percentage is calculated as the percentage of the most
#' frequent nucleotide in a window of \code{-homopolymer_window} to \code{homopolymer_window}
#' nucleotides around the variant position, excluding the variant position itself.
#' Calculation of the homopolymer percentage is skipped when \code{homopolymer_window = 0}.
#' This is useful for filtering out Nanopore sequencing errors in homopolymer regions.
#' @param annotated_region_only logical(1): whether to only consider variants outside
#' annotated regions. If \code{TRUE}, only variants outside annotated regions will be
#' returned. If \code{FALSE}, all variants will be returned, which could take significantly
#' longer time.
#' @param names_from character(1): the column name in the metadata column of the annotation
#' (\code{mcols(annotation)[, names_from]}) to use for the \code{region} column in the output.
#' @return A tibble with columns: seqnames, pos, nucleotide, count, sum, freq, ref, region,
#' homopolymer_pct, bam_path The homopolymer percentage is calculated as the percentage of the
#' most frequent nucleotide in a window of \code{homopolymer_window} nucleotides around
#' the variant position, excluding the variant position itself.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' genome_fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")
#' minimap2_align( # align to genome
#'   config = jsonlite::fromJSON(
#'     system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")),
#'   fa_file = genome_fa,
#'   fq_in = system.file("extdata", "fastq", "demultiplexed.fq.gz", package = "FLAMES"),
#'   annot = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   outdir = outdir
#' )
#' variants <- find_variants(
#'   bam_path = file.path(outdir, "align2genome.bam"),
#'   reference = genome_fa,
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   min_nucleotide_depth = 4
#' )
#' head(variants)
#' @export
find_variants <- function(bam_path, reference, annotation, min_nucleotide_depth = 100,
                          homopolymer_window = 3, annotated_region_only = FALSE,
                          names_from = "gene_name", threads = 1) {
  if (is.character(reference)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading reference ..."))
    reference <- Biostrings::readDNAStringSet(reference)
    # get rid of the `1` from >chr1 1
    names(reference) <- sapply(names(reference), function(x) strsplit(x, " ")[[1]][1])
  }
  if (is.character(annotation)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading annotation ..."))
    annotation <- rtracklayer::import(annotation) |>
      (\(x) x[S4Vectors::mcols(x)$type == "gene", ])()
  }

  if (!annotated_region_only) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Adding unannotated gaps ..."))
    # parsed annotation might not have seqlengths in seqinfo
    if (!all(as.character(GenomeInfoDb::seqnames(annotation)) %in% names(reference))) {
      warning("Some seqnames in annotation not found in reference")
      annotation <- annotation[as.character(GenomeInfoDb::seqnames(annotation)) %in% names(reference)]
    }
    GenomeInfoDb::seqinfo(annotation) <-
      reference[GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(annotation))] |>
      Biostrings::seqinfo()
    if (any(is.na(GenomeInfoDb::seqlengths(annotation)))) {
      stop("Missing seqlengths in seqinfo of annotation")
    }
    annotation <- c(annotation, GenomicRanges::gaps(annotation))
  }

  if (length(bam_path) == 1) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got 1 bam file, parallelizing over each region ..."))
    variants <- tryCatch({
      BiocParallel::bplapply(
        sapply(seq_along(annotation), function(x) annotation[x]),
        function(grange) {
          find_variants_grange(bam_path, reference, grange, min_nucleotide_depth, names_from)
        },
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE # interactive()
        )
      )}, error = identity)
  } else {
    # multi-threading over bam files x granges
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got multiple bam files, preparing for multi-threading ..."))
    args_grid <- expand.grid(
      grange = seq_along(annotation), # does not work on GRanges directly
      bam = bam_path,
      stringsAsFactors = FALSE
    )

    message(paste0(format(Sys.time(), "%H:%M:%S "), "Multi-threading over bam files x ranges ..."))
    variants <- tryCatch({
      BiocParallel::bpmapply(
        function(bam, grange) {
          find_variants_grange(bam, reference, annotation[grange], min_nucleotide_depth, names_from)
        },
        bam = args_grid$bam, grange = args_grid$grange,
        SIMPLIFY = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE # interactive()
        )
      )
    }, error = identity)
  }

  if (inherits(variants, "error")) {
    warning("Error occurred in `find_variants`, returning error object")
    return(variants)
  }

  message(paste0(format(Sys.time(), "%H:%M:%S "), "Merging results ..."))
  variants <- dplyr::bind_rows(variants)

  if (homopolymer_window > 1) {
    if (nrow(variants) == 0) {
      variants <- variants |>
        dplyr::mutate(homopolymer_pct = numeric(0))
    } else {
      message(paste0(format(Sys.time(), "%H:%M:%S "), "Calculating homopolymer percentages ..."))
      variants$homopolymer_pct <- homopolymer_pct(
        reference, variants$seqnames, variants$pos,
        include_alt = FALSE, n = homopolymer_window, threads = threads
      )
    }
  }

  return(variants)
}

#' Relative mutation positions within the gene body
#'
#' Given a set of mutations and a gene annotation, calculate the relative position of each mutation
#' within the gene body. The gene annotation must have the following types: "gene" and "exon".
#' The gene annotation must be for one gene only. The mutations must be within the gene region.
#' The function will merge overlapping exons and calculate the relative position of each mutation
#' within the gene body, excluding intronic regions.
#'
#' @importFrom GenomicRanges GRanges findOverlaps reduce
#' @importFrom IRanges IRanges start end width
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom magrittr equals
#' @importFrom BiocGenerics strand
#'
#' @param mutations either the tibble output from \code{find_variants} or a GRanges object.
#' Make sure to filter it for only the gene of interest.
#' @param annotation_grange GRanges: the gene annotation. Must have the following types: "gene" and "exon".
#' @param verbose logical(1): whether to print messages.
#' @return A numeric vector of relative positions of each mutation within the gene body. Ranging from
#' 0 (start of the gene) to 1 (end of the gene).
relative_mutation_positions_single <- function(mutations, annotation_grange, verbose = TRUE) {

  if (is.data.frame(mutations)) {
    # verify that all mutations are in the same gene
    if (!length(unique(mutations$region)) == 1) {
      errMsg <- sprintf(
        "Incorrent number of unique values in `mutations` for column `region`: %s",
        paste0(unique(mutations$region), collapse = ", ")
      )
      stop(errMsg)
  }
   
    # convert to GRanges of mutations 
    mutations <- GenomicRanges::GRanges(
      seqnames = mutations$seqnames,
      ranges = IRanges::IRanges(start = mutations$pos, end = mutations$pos)
    )
  } else {
    stopifnot("muations must be either a GRanges or a data.frame object" = is(mutations, "GRanges"))
  }

  # verify that the gene annotation is for one gene only
  if (!length(unique(S4Vectors::mcols(annotation_grange)$gene_id)) == 1) {
    errMsg <- sprintf(
      "Incorrent number of gene_id(s) in `annotation_grange`: %s",
      paste0(unique(S4Vectors::mcols(annotation_grange)$gene_id), collapse = ", ")
    )
    stop(errMsg)
  }

  # verify that all mutations are within the gene region
  mutations_ok <- GenomicRanges::findOverlaps(
    mutations,
    subset(annotation_grange, type == "gene")
  ) |>
    S4Vectors::queryHits() |>
    length() |>
    magrittr::equals(length(mutations))
  if (!mutations_ok) {
    stop("Some mutations are not within the gene region")
  }

  # merge overlapping exons
  merged_exons <- annotation_grange |>
    subset(type == "exon") |>
    GenomicRanges::reduce()

  overlaps <- GenomicRanges::findOverlaps(mutations, merged_exons)
  mutations_in_exons <- mutations[S4Vectors::queryHits(overlaps)]
  if (verbose) {
    message(
      sprintf(
        "Found %d mutations in exons out of at total of %d mutations in %s",
        length(mutations_in_exons), length(mutations),
        unique(S4Vectors::mcols(annotation_grange)$gene_name)
      )
    )
  }
  # no mutations in exons
  if (length(mutations_in_exons) == 0) {
    return(numeric(0))
  }

  cumulative_exon_lengths <- cumsum(IRanges::width(merged_exons))

  positions <-
    mapply(
      function(mutation_id, exon_idx) {
        mutation <- mutations_in_exons[mutation_id]
        exon_end <- IRanges::end(merged_exons[exon_idx])
        exon_cumulative_length <- cumulative_exon_lengths[exon_idx]
        relative_position <- (exon_cumulative_length - exon_end + IRanges::start(mutation)) /
          sum(IRanges::width(merged_exons))
        return(relative_position)
      },
      mutation_id = seq_along(mutations_in_exons), exon_idx = S4Vectors::subjectHits(overlaps),
      SIMPLIFY = TRUE
    )

  # flip the positions if the gene is on the negative strand
  flip <- subset(annotation_grange, type == "gene") |>
    BiocGenerics::strand() |>
    as.character() |>
    magrittr::equals("-")
  if (flip) {
    positions <- 1 - positions
  }

  return(positions)
}

#' Relative mutation positions within the gene body
#'
#' Given a set of mutations and gene annotation, calculate the relative position of each mutation
#' within the gene body they are in.
#'
#' @param mutations either the tibble output from \code{find_variants}. It must have columns \code{seqnames},
#' \code{pos}, and a third column for specifying the gene id or gene name. The mutation must be within the gene region.
#' @param annotation Either path to the annotation file (GTF/GFF) or a GRanges object of the gene annotation.
#' @param bin logical(1): whether to bin the relative positions into 100 bins.
#' @param by character(1): the column name in the annotation to match with the gene annotation.
#' E.g. \code{c("region" = "gene_name")} to match the `region` column in the mutations with the
#' `gene_name` column in the annotation.
#' @param threads integer(1): number of threads to use.
#' @return If \code{bin = FALSE}, a list of numeric vectors of relative positions of each mutation within the gene body.
#' If \code{bin = TRUE}, a numeric vector of length 100 of the number of mutations in each bin.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' genome_fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")
#' minimap2_align( # align to genome
#'   config = jsonlite::fromJSON(
#'     system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")),
#'   fa_file = genome_fa,
#'   fq_in = system.file("extdata", "fastq", "demultiplexed.fq.gz", package = "FLAMES"),
#'   annot = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   outdir = outdir
#' )
#' variants <- find_variants(
#'   bam_path = file.path(outdir, "align2genome.bam"),
#'   reference = genome_fa,
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   min_nucleotide_depth = 4
#' )
#' positions <- 
#'  relative_mutation_positions(
#'    mutations = variants,
#'    annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
#'  )
#' @importFrom rtracklayer import
#' @importFrom dplyr mutate filter
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom S4Vectors mcols
#' @export
relative_mutation_positions <- function(mutations, annotation, bin = FALSE, by = c("region" = "gene_name"), threads = 1){
  if (!length(by) == 1) {
    stop("by must be a character vector of length 1")
  }
  if (is.null(names(by))) {
    names(by) <- by
  }

  if (is.character(annotation)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading annotation ..."))
    annotation_grange <- rtracklayer::import(annotation)
  } else {
    annotation_grange <- annotation
  }

  mutations <- mutations |> 
    dplyr::mutate(region = mutations[, names(by)]) |>
    dplyr::filter(!is.na(region))
  mutations_split <- split(mutations, mutations$region)

  coverages <- tryCatch({
    BiocParallel::bplapply(
      names(mutations_split),
      function(x) {
        annot_i <- subset(annotation_grange, S4Vectors::mcols(annotation_grange)[, by] == x)
        pos <- relative_mutation_positions_single(mutations_split[[x]], annot_i, verbose = FALSE)
        if (!bin) {
          return(pos)
        }
        pos <- round(pos * 100)
        coverage <- sapply(1:100, function(i) sum(pos == i))
        return(coverage)
      },
      BPPARAM = BiocParallel::MulticoreParam(
        workers = threads, stop.on.error = TRUE, progressbar = TRUE
      )
    )
  }, error = identity)
  if (inherits(coverages, "error")) {
    warning("Error occurred in `mutations_coverage`, returning error object")
    return(coverages)
  }

  if (!bin) {
    names(coverages) <- names(mutations_split)
    return(coverages)
  }

  # sum up all coverages
  total_coverage <- do.call(rbind, coverages) |>
    colSums()
  return(total_coverage)
}
