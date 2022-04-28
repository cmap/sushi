# sushi
Data Processing Pipeline for PRISM Sequencing

## Install

Package can be installed using `devtools` by running the following command:

```r
devtools::install_github('https://github.com/cmap/sushi')
```

Note that the package requires `ShortRead`from BioconductR which may not be automatically installed. 
Instructions for installing `ShortRead` can be found on the [BioconductR page](https://bioconductor.org/packages/release/bioc/html/ShortRead.html)

## Setup

Each PrismSeq project should have its own folder within which to process data and save output files. For the purposes of this tutorial, we will assume you are working with the __example_project__ folder. To begin, your project folder must contains a _sample_meta.csv_ file and a __fastq__ folder with 3 files per sample (index1, index2, read).

A template for the _sample_meta.csv_ file may be found here: https://docs.google.com/spreadsheets/d/1t0Avob53rSio4qcxb5QFqnFjRSZC3wnb/edit?usp=sharing&ouid=112283500068607320752&rtpof=true&sd=true. Please contact cmapa@broadinsitute.org for access, if needed. 
Create a copy of the _Sequencing Metadata Template [internal]_ tab, fill it out according to the specifications in the instructions tab, and export the metadata sheet as a csv to your project directory.

## Run
In order to run the PrismSeq processing pipeline, you must have copies of the follow three metadata files. For the purposes of this tutorial, we will assume they are saved in the __metadata__ folder.
1. _cell_line_meta.csv_
2. _cell_set_meta.csv_
3. _CB_meta.csv_
Please contact cmapa@broadinstitue.org to request metadata files.

For additional information on the outputs of the PrismSeq processing pipeline see the notes and FAQs here: https://docs.google.com/document/d/1sHpkXQzzFu63QbXYc4W3_YmBCqBikQKXmj57yBJFiCw/edit?usp=sharing. Please contact cmapa@broadinstitute.org for access, if needed.

### R functions
The PrismSeq processing pipeline may be run in R using the following functions.

1. Generate a read count table from fastq files
For assitance in generating fastq files from BCL files please contact anup@broadinstitute.org, if needed.

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

```r
# load relevant metadata
sample_meta = read.csv("example_project/sample_meta.csv")
cell_line_meta = read.csv("metadata/cell_line_meta.csv")
cell_set_meta = read.csv("metadata/cell_set_meta.csv") 
CB_meta = read.csv("metadata/CB_meta.csv")

filered_counts = filter_raw_reads(raw_counts, 
                                  sample_meta, 
                                  cell_line_meta, 
                                  cell_set_meta, 
                                  CB_meta, 
                                  id_cols=c('cell_set', 'treatment', 'dose','dose_unit','day','bio_rep','tech_rep')) # change the id_cols parameter to designate which metadata column uniquely define each profile

QC_table = filtered_counts$qc_table
filtered_counts = filtered_counts$filtered_counts

QC_table %>% write.csv("example_project/QC_table.csv", row.names=F, quote=F)
filtered_counts %>% write.csv("example_project/filtered_counts.csv", row.names=F, quote=F)
```

3. [Optional] Generate QC images
```r
filteredCounts_QC(filtered_counts, 
                  cell_set_meta) 
```

4. [Optional] Run DEseq on filtered counts
Your project must include control barcodes in order for DESeq to run.

```r
# define DESeq input values. Change sample_cols and id_cols to reflect which columns uniquely define each sample and signature
sample_cols = c("cell_set", "treatment", "dose", "dose_unit", "day", "bio_rep")
sig_cols = c("cell_set", "treatment", "dose", "dose_unit", "day")
control_type = "negcon"

l2fc_DEseq = data.frame()
for(cs in unique(filtered_counts$cell_set)) {
  subset = filtered_counts %>% 
    filter(cell_set==cs)
  hold = run_DE(subset,
                sample_cols,
                sig_cols,
                control_type)
  l2fc_DEseq = l2fc_DEseq %>% 
    rbind(hold)
}

l2fc_DEseq %>% write.csv("example_project/l2fc_DEseq.csv", row.names=F, quote=F)
```

5. [Optional] Normalize filtered counts
Your project must include control barcodes if you want to normalize your counts.
```r
normalized_counts = normalize(filtered_counts, 
                              CB_meta$Name)

normalized_counts %>% write.csv("example_project/normalized_counts.csv", row.names=F, quote=F)
```

6. Generate log-fold change values
If you normalized your counts, log-fold change values should be generated from normalized counts, otherwise they should be generated from basic filtered counts. Here we assume that you have generated normalized counts as in the previous step.

```r
l2fc = compute_l2fc(normalized_counts,
                    control_type = "negcon",
                    sig_cols=c('cell_set', 'treatment', 'dose','dose_unit','day'))

l2fc %>% write.csv("example_project/l2fc.csv", row.names=F, quote=F)
```

7. Collapse replicates
```r
collapsed_values = collapse_counts(l2fc)

collapsed_values %>% write.csv("example_project/collapsed_values.csv", row.names=F, quote=F)
```

8. Run biomarker analysis
This function uses data which is stored on taiga. If you are a member of the Broad Institute you can install taigr, the taiga client for R, by following the instructions here: https://github.com/broadinstitute/taigr. 
This function can take several hours to run.

```r
biomarkers = generate_biomarkers(collapsed_values)

lin_table = biomarkers$lin_table
rf_table = biomarkers$rf_table
disc_table = biomarkers$disc_table

lin_table %>% write.csv("example_project/lin_table.csv", row.names=F, quote=F)
rf_table %>% write.csv("example_project/rf_table.csv", row.names=F, quote=F)
disc_table %>% write.csv("example_project/disc_table.csv", row.names=F, quote=F)
```

### Command line tools
The PrismSeq processing tools may alternatively be run through the command line. For the purposes of this tutorial, we will assume that the PrismSeq processing scripts are saved in the __tools/scripts/__ folder.

`Rscript [toolname] --help` will provide information about arguments for specified tools.

The same considerations as above apply to the following tools.

1. Generate a read count table from fastq files
`Rscript tools/scripts/fastq2readcounts.R --fastq fastq\ --index_1 _I1_ --index_2 _I2_ --barcode_suffix _R1_`

2. Filter raw read counts
`Rscript tools/scripts/filter_counts.R --raw_counts raw_counts.csv --sample_meta sample_meta.csv --id_cols cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep --cell_line_meta metadata/cell_line_meta.csv --cell_set_meta metadata/cell_set_meta.csv --CB_meta metadata/CB_meta.csv`

3. [Optional] Generate QC images
`Rscript filteredCounts_QC.R --filtered_counts filtered_counts.csv --cell_set_meta metadata/cell_set_meta.csv`

4. [Optional] Run DEseq on filtered counts
`Rscript tools/scripts/call_run_DE.R --filtered_counts filtered_counts.csv --sample_cols cell_set,treatment,dose,dose_unit,day,bio_rep --sig_cols cell_set,treatment,dose,dose_unit,day --control_type negcon`

5. [Optional] Normalize filtered counts
`Rscript tools/scripts/CBnormalize.R --filtered_counts filtered_counts.csv --CB_meta metadata/CB_meta.csv`

6. Generate log-fold change values
`Rscript tools/scripts/compute_l2fc.R --normalized_counts normalized_counts.csv --control_type negcon`

7. Collapse replicates
`Rscript tools/scripts/collapse_replicates.R --lfc l2fc.csv`

8. Run biomarker analysis
`Rscript tools/scripts/generate_biomarkers.R --collapsed_values collapsed_values.csv`








