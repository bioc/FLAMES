% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minimap2_align.R
\name{minimap2_realign}
\alias{minimap2_realign}
\title{Minimap2 re-align reads to transcriptome}
\usage{
minimap2_realign(
  config,
  fq_in,
  outdir,
  minimap2,
  samtools = NULL,
  prefix = NULL,
  minimap2_args,
  sort_by,
  threads = 1
)
}
\arguments{
\item{config}{Parsed list of FLAMES config file}

\item{fq_in}{File path to the fastq file used as a query sequence file}

\item{outdir}{Output folder}

\item{minimap2}{Path to minimap2 binary}

\item{samtools}{path to the samtools binary, required for large datasets since \code{Rsamtools} does not support \code{CSI} indexing}

\item{prefix}{String, the prefix (e.g. sample name) for the outputted BAM file}

\item{minimap2_args}{vector of command line arguments to pass to minimap2}

\item{sort_by}{String, If provided, sort the BAM file by this tag instead of by position.}

\item{threads}{Integer, threads for minimap2 to use, see minimap2 documentation for details,
FLAMES will try to detect cores if this parameter is not provided.}
}
\value{
a \code{data.frame} summarising the reads aligned
}
\description{
Uses minimap2 to re-align reads to transcriptome
}
\examples{
outdir <- tempfile()
dir.create(outdir)
annotation <- system.file('extdata', 'rps24.gtf.gz', package = 'FLAMES')
genome_fa <- system.file('extdata', 'rps24.fa.gz', package = 'FLAMES')
fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
fastq <- system.file('extdata', 'fastq', 'demultiplexed.fq.gz', package = 'FLAMES')
minimap2_realign(
  config = jsonlite::fromJSON(
    system.file("extdata", "config_sclr_nanopore_3end.json", package = 'FLAMES')
  ),
  fq_in = fastq,
  outdir = outdir
)
}
\seealso{
[minimap2_align()]
}
