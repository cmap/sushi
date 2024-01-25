# Sushi
Data Processing Pipeline for PRISM Sequencing.

## Install

Package can be installed using `devtools` by running the following command:

```r
devtools::install_github('https://github.com/cmap/sushi')
```

Note that the package requires `ShortRead`from BioconductR which may not be automatically installed. 
Instructions for installing `ShortRead` can be found on the [BioconductR page](https://bioconductor.org/packages/release/bioc/html/ShortRead.html)

## Setup

Each PrismSeq project should have its own folder within which to process data and save output files. For the purposes of this tutorial, we will assume you are working with the __example_project__ folder. To begin, your project folder must contains a _sample_meta.csv_ file and a __fastq__ folder with 3 files per sample (index1, index2, read).

A template for the _sample_meta.csv_ file may be found [here](https://docs.google.com/spreadsheets/d/1t0Avob53rSio4qcxb5QFqnFjRSZC3wnb/edit?usp=sharing&ouid=112283500068607320752&rtpof=true&sd=true). Create a copy of the _Sequencing Metadata Template [internal]_ tab, fill it out according to the specifications in the instructions tab, and export the metadata sheet as a csv to your project directory.

Please contact prism@broadinsitute.org for access, if needed. 

## Running the pipeline

In order to run the PrismSeq processing pipeline, you must have copies of the follow three metadata files. For the purposes of this tutorial, we will assume they are saved in the __metadata__ folder.
1. _cell_line_meta.csv_
2. _cell_set_meta.csv_
3. _CB_meta.csv_

For additional information on the outputs of the PrismSeq processing pipeline see the notes and FAQs [here](https://docs.google.com/document/d/1sHpkXQzzFu63QbXYc4W3_YmBCqBikQKXmj57yBJFiCw/edit?usp=sharing). 

Please contact prism@broadinstitue.org to request metadata files and for access to the notes and FAQs document, if needed.

### Using R functions

The PrismSeq processing pipeline may run in R using the following functions.

1. Generate a read count table from fastq files

The function _write_df_from_fastq_ reads in fastq files associated with the project and organizes the reads into a table or dataframe saved as _raw_counts.csv_.

For assistance in generating fastq files from BCL files please contact prism@broadinstitute.org, if needed.

```r
# define invariant strings to identify fastq files
index_1 = "_I1_"
index_2 = "_I2_"
barcode_suffix = "_R1_"

read_directory_contents = c("example_project/fastq/") %>% 
  purrr::map(list.files, full.names = T) %>%
  purrr::reduce(union)

barcode_read_files = read_directory_contents %>%
  purrr::keep(stringr::str_detect, fixed(barcode_suffix)) %>%
  sort()

index_1_files = read_directory_contents %>%
  purrr::keep(stringr::str_detect, fixed(index_1)) %>%
  sort()

index_2_files = read_directory_contents %>%
  purrr::keep(stringr::str_detect, fixed(index_2)) %>%
  sort()

raw_counts = write_df_from_fastq(forward_read_fastq_files = barcode_read_files, 
                                 index_1_files = index_1_files, 
                                 index_2_files = index_2_files,
                                 write_interval = NA) # define write_interval to intermitently write out raw_counts file

raw_counts %>% write.csv("example_project/raw_counts.csv", row.names=F, quote=F)
```

2. Filter raw read counts

The function _filter_raw_reads_ filters _raw_counts.csv_ for only the reads associated with the project as defined in the submitted sample meta.

```r
# load relevant metadata
sample_meta = read.csv("example_project/sample_meta.csv")
cell_line_meta = read.csv("metadata/cell_line_meta.csv")
cell_set_meta = read.csv("metadata/cell_set_meta.csv") 
CB_meta = read.csv("metadata/CB_meta.csv")

filtered_counts = filter_raw_reads(raw_counts, 
                                   sample_meta, 
                                   cell_line_meta, 
                                   cell_set_meta, 
                                   CB_meta, 
                                   id_cols=c('cell_set', 'treatment', 'dose','dose_unit','day','bio_rep','tech_rep')) # change the id_cols parameter to designate which metadata column uniquely define each profile

QC_table = filtered_counts$qc_table
annotated_counts= filtered_counts$annotated_counts
filtered_counts = filtered_counts$filtered_counts

QC_table %>% write.csv("example_project/QC_table.csv", row.names=F, quote=F)
annotated_counts %>% write.csv("example_project/annotated_counts.csv", row.names=F, quote=F)
filtered_counts %>% write.csv("example_project/filtered_counts.csv", row.names=F, quote=F)
```

3. [Optional] Normalize filtered counts

The function _normalize_ fits a linear model to a set of control barcodes (usually 10) supplied at specific doses. This model is then used to calculate the equivalent dose of cell line barcodes detected in each well. The normalized values are in the units of control barcode doses, and are not interpretable on their own.

This function can only be used if control barcodes are included in the run. If control barcodes were not used, this function can be skipped.

```r
normalized_counts = normalize(filtered_counts, 
                              CB_meta$Name)

normalized_counts %>% write.csv("example_project/normalized_counts.csv", row.names=F, quote=F)
```

4. Generate log-fold change values

The function _l2fc_ computes the log2 fold change between treatments and negative controls. This is done by collecting the negative control and collapsing across all replicates. The negative controls are then joined with the treatments to calculate log2 fold change values for each treatment replicate. If _normalize_ was called, then log2 fold change values should be generated from normalized counts. If _normalize_ was not called, log2 fold change values can be calculated from the read counts in _filtered_counts.csv_ by adjusting the _count_col_name_ parameter. 

The example code below assues that the previous step, _normalize_, was performed.

```r
l2fc = compute_l2fc(normalized_counts,
                    control_type = "negcon",
                    sig_cols=c('cell_set', 'treatment', 'dose','dose_unit','day'),
                    count_col_name="normalized_n") # change based on whether you are running compute_l2fc on normalized_counts or filtered_counts

l2fc %>% write.csv("example_project/l2fc.csv", row.names=F, quote=F)
```

5. Collapse replicates

This function collapses the biological replicates of the treatments to get a single value for each unique cell line and treatment combination. 

```r
collapsed_values = collapse_counts(l2fc)

collapsed_values %>% write.csv("example_project/collapsed_values.csv", row.names=F, quote=F)
```

6. [Optional] Generate QC images

This function outputs several files and images that can be used to assess the quality of the run. 

```r
QC_images(sample_meta,
          annotated_counts,
          filtered_counts,
          cell_set_meta,
          out = "example_project/") 
```

### Using command line tools

The PrismSeq processing tools may alternatively be run through the command line. For the purposes of this tutorial, we will assume that the PrismSeq processing scripts are saved in the __tools/scripts/__ folder.

`Rscript [toolname] --help` will provide information about arguments for specified tools.

The same considerations as above apply to the following tools.

1. Generate a read count table from fastq files

```
Rscript tools/scripts/fastq2readcounts.R --fastq fastq\ --index_1 _I1_ --index_2 _I2_ --barcode_suffix _R1_
```

2. Filter raw read counts

```
Rscript tools/scripts/filter_counts.R --raw_counts raw_counts.csv --sample_meta sample_meta.csv --id_cols cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep --cell_line_meta metadata/cell_line_meta.csv --cell_set_meta metadata/cell_set_meta.csv --CB_meta metadata/CB_meta.csv
```

3. [Optional] Normalize filtered counts

```
Rscript tools/scripts/CBnormalize.R --filtered_counts filtered_counts.csv --CB_meta metadata/CB_meta.csv
```

4. Generate log-fold change values

```
Rscript tools/scripts/compute_l2fc.R --normalized_counts normalized_counts.csv --control_type negcon --count_col_name normalized_n
```

5. Collapse replicates

```
Rscript tools/scripts/collapse_replicates.R --lfc l2fc.csv
```

6. [Optional] Generate QC images

```
Rscript filteredCounts_QC.R --sample_meta sample_meta.csv --annotated_counts annotated_counts.csv --filtered_counts filtered_counts.csv --cell_set_meta metadata/cell_set_meta.csv
```

## Other

This document was last update on Jan 25, 2024
