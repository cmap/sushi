% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_df_from_fastq.R
\name{write_df_from_fastq}
\alias{write_df_from_fastq}
\title{write_df_from_fastq}
\usage{
write_df_from_fastq(
  forward_read_fastq_files,
  index_1_files,
  index_2_files,
  write_interval = NA
)
}
\arguments{
\item{forward_read_fastq_files}{\itemize{
\item vector of fastq file paths
}}

\item{index_1_files}{\itemize{
\item vector of fastq file paths
}}

\item{index_2_files}{\itemize{
\item vector of fastq file paths
}}

\item{write_interval}{\itemize{
\item integer for how often a temp count file is written, NA by default.
}}
}
\value{
\itemize{
\item cumulative_count_df A data.frame of readcounts by index_1, index_2 and forward_read_cl_barcode
}
}
\description{
takes PRISM miseq or hiseq fastqs and returns a dataframe
with the number of barcodes from each cell line in each well
}
