#' Pipeline for Bulk Data
#'
#' @md
#' @description Semi-supervised isofrom detection and annotation for long read data.
#' This variant is meant for bulk samples. Specific parameters relating to
#' analysis can be changed either through function arguments, or through a
#' configuration JSON file.
#' @inherit sc_long_pipeline details description
#'
#' @return if \code{do_transcript_quantification} set to true, \code{bulk_long_pipeline} returns a SummarizedExperiment object, containing a count
#' matrix as an assay, gene annotations under metadata, as well as a list of the other
#' output files generated by the pipeline. The pipeline also outputs a number of output
#' files into the given \code{outdir} directory. These output files generated by the pipeline are:
#' \itemize{
#'  \item{transcript_count.csv.gz}{ - a transcript count matrix (also contained in the SummarizedExperiment)}
#'  \item{isoform_annotated.filtered.gff3}{ - isoforms in gff3 format (also contained in the SummarizedExperiment)}
#'  \item{transcript_assembly.fa}{ - transcript sequence from the isoforms}
#'  \item{align2genome.bam}{ - sorted BAM file with reads aligned to genome}
#'  \item{realign2transcript.bam}{ - sorted realigned BAM file using the transcript_assembly.fa as reference}
#'  \item{tss_tes.bedgraph}{ - TSS TES enrichment for all reads (for QC)}
#' }
#' if \code{do_transcript_quantification} set to false, nothing will be returned
#'
#' @inheritParams sc_long_pipeline
#'
#' @seealso
#' [sc_long_pipeline()] for single cell data,
#' [SummarizedExperiment()] for how data is outputted
#'
#' @example inst/examples/pipeline_example.R
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData rowData<- colData<- rowRanges rowRanges<-
#' @importFrom utils read.csv read.table
#' @importFrom dplyr group_by summarise_at slice_max filter
#' @importFrom magrittr "%>%"
#' @importFrom BiocGenerics cbind colnames rownames start end
#' @importFrom Rsamtools indexBam
#' @export
bulk_long_pipeline <-
    function(annotation,
             fastq,
             outdir,
             genome_fa,
             minimap2 = NULL,
             k8 = NULL,
             config_file = NULL) {
        checked_args <- check_arguments(
            annotation,
            fastq,
            genome_bam = NULL,
            outdir,
            genome_fa,
            config_file
        )

        config <- checked_args$config

        # create output directory if one doesn't exist
        if (!dir.exists(outdir)) {
            cat("Output directory does not exists: one is being created\n")
            dir.create(outdir)
            print(outdir)
        }

        if (utils::file_test("-d", fastq)) {
            fastq_files <- file.path(fastq, list.files(fastq))
            fastq_files <- fastq_files[grepl("\\.(fastq|fq)(\\.gz)?$", fastq_files) & utils::file_test("-f", fastq_files)]
        } else if (utils::file_test("-f", fastq)) {
            fastq_files <- fastq
        } else {
            stop("fastq must be a valid path to a folder or a FASTQ file")
        }
        samples <- gsub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_files))

        using_bam <- FALSE
        genome_bam <- file.path(outdir, paste0(samples, "_", "align2genome.bam"))
        if (all(utils::file_test("-f", genome_bam))) {
            cat("Found all corresponding '[sample]_align2genome.bam' files, will skip initial alignment.\n")
            using_bam <- TRUE
            config$pipeline_parameters$do_genome_alignment <- FALSE
            if (!all(
                utils::file_test("-f", file.path(outdir, paste0(samples, "_", "align2genome.bam.bai"))) |
                utils::file_test("-f", file.path(outdir, paste0(samples, "_", "align2genome.bam.csi"))))) {
                for (bam in genome_bam) {
                    Rsamtools::indexBam(bam)
                }
            }
        }

        cat("#### Input parameters:\n")
        cat(jsonlite::toJSON(config, pretty = TRUE), "\n")
        cat("gene annotation:", annotation, "\n")
        cat("genome fasta:", genome_fa, "\n")
        cat("input fastq files:", gsub("$", "\n", fastq_files))
        cat("output directory:", outdir, "\n")
        cat("minimap2 path:", minimap2, "\n")
        cat("k8 path:", k8, "\n")

        if (config$pipeline_parameters$do_genome_alignment) {
            cat("#### Aligning reads to genome using minimap2\n")
            for (i in 1:length(samples)) {
                cat(paste0(c("\tAligning sample ", samples[i], "...\n")))
                minimap2_align(
                    config,
                    genome_fa,
                    fastq_files[i],
                    annotation,
                    outdir,
                    minimap2,
                    k8,
                    prefix = samples[i],
                    threads = config$pipeline_parameters$threads
                )
            }
        } else {
            cat("#### Skip aligning reads to genome\n")
        }

        # find isofroms
        if (config$pipeline_parameters$do_isoform_identification) {
            find_isoform(annotation, genome_fa, genome_bam, outdir, config)
        }

        # realign to transcript
        if (config$pipeline_parameters$do_read_realignment) {
            cat("#### Realign to transcript using minimap2\n")
            for (i in 1:length(samples)) {
                cat(paste0(c("\tRealigning sample ", samples[i], "...\n")))
                minimap2_realign(config, fastq_files[i], outdir, minimap2, prefix = samples[i], 
                                 threads = config$pipeline_parameters$threads)
            }
        } else {
            cat("#### Skip read realignment\n")
        }

        # quantification
        if (config$pipeline_parameters$do_transcript_quantification) {
            cat("#### Generating transcript count matrix\n")
            quantify_transcript(annotation = annotation, outdir = outdir, config = config, pipeline = "bulk")

            out_files <- list(
                "annotation" = annotation,
                "genome_fa" = genome_fa,
                "counts" = file.path(outdir, "transcript_count.csv.gz"),
                "isoform_annotated" = file.path(outdir, "isoform_annotated.filtered.gff3"),
                "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
                "align_bam" = genome_bam,
                "realign2transcript" = file.path(outdir, list.files(outdir))[grepl("realign2transcript\\.bam$", list.files(outdir))],
                "tss_tes" = file.path(outdir, "tss_tes.bedgraph"),
                "outdir" = outdir
            )
            load_genome_anno <- rtracklayer::import(annotation, feature.type = c("exon", "utr"))

            se <- generate_bulk_summarized(out_files, load_genome_anno = load_genome_anno)
            return(se)
        } else {
            cat("#### Skip transcript quantification\n")
        }
    }
