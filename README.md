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

## PRISM EPS sequencing data processing

In contrast to the conventional MTS assay which uses Luminex detection, the Extended Assay (EPS) uses sequencing-based detection. 

![EPS Pipeline (1)](https://github.com/cmap/sushi/assets/125501149/0d783769-1e6a-4077-8c4a-e61fb40ebbd1)

Raw data for EPS is thus first processed using SUSHI (yellow), the PRISM sequencing data processing pipeline, to obtain normalized counts and log fold-change data. 

This data is then passed to the MTS pipeline (orange) for downstream processing such as dose-response fitting and biomarker analysis via the EPS conversion modules (green). Documentation for the MTS pipeline can be found [here](https://github.com/cmap/dockerized_mts/tree/master). 

The processed data and visualizations are then available on the PRISM Portal (blue). 

## Running the Sushi pipeline

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

## Running the EPS conversion modules

![EPS Pipeline (1)](https://github.com/cmap/sushi/assets/125501149/728ca2b7-8304-44cc-8aef-304dea9c396b)

When processing data from the Extended Day assay, the EPS conversion modules are also run. The 2 EPS conversion modules (green) perform separate tasks:

1. Quality metric calculations ( _EPS_QC.R_ )

2. Adding or renaming the column headers of the SUSHI output for downstream compatibility ( _seq_to_mts.py_ )


__Files generated by the EPS pipeline that are shared with collaborators__

1. _EPS_QC_table.csv_: Contains a table with day 10 count information and quality control metrics for each cell line indicating if the data passed raw count QC.

2. _LEVEL3_NORMALIZED_COUNTS.csv_: Contains normalized counts in EPS format, generated by running a conversion module to align headers to our downstream EPS pipeline.

3. _LEVEL4_LFC.csv_: Contains log-fold change data in EPS format, generated by running a conversion module to align headers to our downstream EPS pipeline.

4. _LEVEL5_LFC.csv_: Contains collapsed log fold change data in EPS format, generated by running a conversion module to align headers to our downstream EPS pipeline.


__Intermediate files generated by the EPS pipeline that are for internal use (not shared)__

1. _normalized_counts.csv_: Contains normalized counts, generated by fitting a linear model to a set of 10 control barcodes. Used to generate LEVEL3_NORMALIZED_COUNTS.csv.

2. _l2fc.csv_: Contains log-fold change values for each biological replicate you submitted, as compared to your annotated negative control. Used to generate LEVEL4_LFC.csv.

3. _collapsed_values.csv_: Contains median-collapsed log-fold change values.Used to generate LEVEL5_LFC.csv.

4. _compound_key.csv_: Contains an overview of the perturbations, projects, and doses. Used by EPS conversion modules.

5. _inst_info.csv_: Contains a comprehensive summary of compound, dose, cell line, and experimental information. Used by EPS conversion modules.


## Conversion of Sushi column headers into EPS format 
Column header conversions are outlined below for counts and LFC data, QC, and compound key files. Column headers which are present in only select levels of data are specified in the Description column (i.e. _LEVEL4_LFC_, _LEVEL5_LFC_, etc). Headers which are renamed during the SUSHI to EPS Data Conversion show the original name in _SUSHI Column_ and the rename in _EPS Column_. Headers which are generated during the SUSHI to EPS Data Conversion or downstream in the pipeline have no entries for _SUSHI Column_.


__Column names in LEVEL3_NORMALIZED_COUNTS, LEVEL4_LFC, and LEVEL5_LFC:__

| __SUSHI Column__ | __EPS Column__ | __Description__ |
|:-------------|:-----------|:------------|
| bio_rep      | replicate  | Plate replicate |
| cb_intercept | cb_intercept | Intercept of the control barcode normalization |
| CCLE_name    | ccle_name  | Cell line name |
| cell_set     | cell_set   | Cell set information |
| control_barcode | control_barcode | Presence of control barcode in a given well |
| control_MAD_QC | control_MAD_QC | TRUE if the variance across the biological replicates in the negative control is below threshold of ~1.66 (LEVEL4_LFC) |
| control_mad_sqrtN | control_mad_sqrtN | Median Absolute Deviation over the square root of the number of negative controls (LEVEL4_LFC) |
| control_median_n | control_median_n | Median collapsed n across the biological replicates of the negative controls (LEVEL4_LFC) |
| control_median_normalized_n | control_median_normalized_n | Median normalized n across the biological replicates of the negative control (LEVEL4_LFC) |
| counts_flag | counts_flag | This flags entries where the collapsed raw count in the negative control is below a threshold as “negcon<__”. (LEVEL4_LFC) |
| DepMap_ID | depmap_id | DepMap cell line ID |
|   | feature_id | Culture concatenated with cell line name |
| flag | flag | This column is filled with “Missing” if the number of reads n is zero. The column is filled with “low counts” if the number of reads n is below the count threshold. This threshold defaults to 40. |
| log2_dose | log2_dose | In reference to the control barcode |
| log2_n | log2_n | log2(n+pseudocount). A pseudocount (default of 20) is added so that missing cell lines can be carried through the normalization. |
| log2_normalized_n | log2_normalized_n | Log normalized count including pseudocount |
| mean_n | mean_n | Mean of n across technical reps (LEVEL4_LFC) |
| mean_normalized_n | mean_normalized_n | Mean of normalized_n across technical reps (LEVEL4_LFC) |
| n | n | Number of Reads |
| Name | Name | Control barcode name |
| normalized_n | normalized_n | Normalized read count |
| num_bio_reps | num_bio_reps | Number of biological replicates of treatment that were collapsed. (LEVEL5_LFC) |
| num_ctrl_bio_reps | num_ctrl_bio_reps | Number of biological replicates that were collapsed in the negative controls. (LEVEL4_LFC) |
| num_tech_reps | num_tech_reps | Number of technical replicates that were collapsed. (LEVEL4_LFC) |
| dose | pert_dose | Perturbation dose (numeric) |
| dose_unit | pert_dose_unit | Perturbation dose units |
| treatment | pert_id | Uppercase of `treatment` values. Perturbation Broad ID. |
|   | pert_idose | Perturbation dose with units |
| treatment | pert_iname | Perturbation name (e.g. AZ-628) |
| day | pert_itime | Assay length with units |
|   | pert_plate | Sample compound plate |
| day | pert_time | Assay length |
|   | pert_time_unit | Assay length units |
| trt_type | pert_type | Perturbation type (e.g. negative control) ('poscon': 'trt_poscon', 'negcon': 'ctl_vehicle') |
|   | Perturbation vehicle (e.g. DMSO) |
| pcr_plate | pcr_plate | PCR plate information |
| pcr_well | pcr_well | PCR well information |
| pcr_well | pert_well | Sample well |
|   | pool_id | Assay pool |
| prism_cell_set | culture | Cell line culture (e.g. PR500) |
| profile_id | profile_id | Concatenation of cell set, replicate, compound, and dose information. |
| treatment | prc_id | Uppercase of `treatment` values. Perturbation Broad ID. |
| x_project_id | x_project_id | Project name (e.g. Validation Compounds) |
| screen | screen | Screen type and iteration |
| tech_rep | tech_rep | Technical replicate |
| trt_MAD_QC | trt_MAD_QC | TRUE if the variance across the biological replicates in the negative control is below a threshold of ~1.66. (LEVEL5_LFC) |
| trt_mad_sqrtN | trt_mad_sqrtN | Median Absolute Ddeviation over the square root of the number of treatment biological replicates. (LEVEL5_LFC) |
| trt_median_n | trt_median_n | Median collapsed n across the biological replicates of the treatments. (LEVEL5_LFC) |
| trt_median_normalized_n | trt_median_normalized_n | Median normalized n across the biological replicates of the treatments. (LEVEL5_LFC) |
|   | feature_id | Unique feature ID |
| l2fc | LFC | log2 fold-change relative to control vehicle (LEVEL4_LFC) |
| median_l2fc | LFC | median log2 fold-change relative to control vehicle (LEVEL5_LFC) |


__Column names in EPS QC file:__

| __Column__              | __Description__                                               |
|:--------------------|:----------------------------------------------------------|
| CCLE_name           | Cancer cell line encyclopedia cell line name             |
| DepMap_ID           | DepMap cell line ID                                       |
| cell_set            | Cell set information                                      |
| day                 | Day at end of experiment (0, 6, 10, etc)                  |
| med_log_raw_counts  | Median log raw counts                                     |
| med_log_norm_counts | Median log normalized counts                              |
| med_raw_counts      | Median raw counts                                         |
| pass_raw_count_qc   | True/False value indicating whether raw counts passed QC  |
| pert_plate          | Pert plate information                                    |


__Column names in compound key file:__

| __Column__        | __Description__                                     |
|:--------------|:------------------------------------------------|
| pert_iname    | Compound                                        |
| pert_id       | Compound ID                                     |
| pert_plate    | Pert plate                                      |
| x_project_id  | Project code                                    |
| pert_dose_1   | Number of unique doses per compound per plate  |


## Other

This document was last update on Feb 25, 2024
