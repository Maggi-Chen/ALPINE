# ALPINE: AAV integration evaluation using targeted long-read sequencing data

This repo represents ALPINE pipeline for AAV integration evaluation using targeted long-read
sequencing data from both PacBio HiFi and Oxford Nanopore platforms. It contains Docker files, CWL workflows, multi-platform launcher code, and other scripts used for
on-target AAV integration evaluation.

## Table of Contents

- [Overview](#overview)
- [Steps](#steps)
  - [Per-sample AAV analysis](#per-sample-aav-analysis)
    - [Read filtering](#read-filtering)
    - [Read alignment](#read-alignment)
    - [Read classification](#read-classification)
    - [Read counting](#read-counting)
  - [Merge classification table](#merge-classification-table)
- [Inputs](#inputs)
  - [Sequencing FASTQ](#sequencing-fastq)
  - [Reference sequences](#reference-sequences)
  - [AAV vector config](#aav-vector-config)
  - [Primer sequences](#primer-sequences)
  - [Sample name](#sample-name)
  - [Optional inputs](#optional-inputs)
  - [Inputs for merge classification table step](#inputs-for-merge-classification-table-step)
- [Test Data and Examples](#test-data-and-examples)
- [Outputs](#outputs)
  - [Per-sample outputs](#per-sample-outputs)
    - [Filtered FASTQ](#filtered-fastq)
    - [Sorted BAM](#sorted-bam)
    - [Read annotation](#read-annotation)
    - [Single-sample classification table](#single-sample-classification-table)
  - [Merged classification table](#merged-classification-table)
  - [ALPINE Output Classification Categories](#alpine-output-classification-categories)
- [How to run pipeline](#how-to-run-pipeline)
- [Update workflows in SBG and Arvados](#update-workflows-in-sbg-and-arvados)
- [Contacts](#contacts)

## Overview

ALPINE (AAV Long-read Integration Analysis Pipeline) evaluates on-target AAV integration efficiency by analyzing targeted long-read sequencing data from PacBio HiFi and Oxford Nanopore platforms.

## Steps

![plot](./workflow_pipeline.png)

The ALPINE analysis pipeline mainly includes two steps: per-sample AAV analysis step and merge count table step. Per-sample
AAV analysis workflow is run on all samples individually to classify reads into multiple categories and count reads in each
category. After all samples are processed through per-sample workflow, merge count table step merges results from all
samples together into a merged table.

### Per-sample AAV analysis

Per-sample AAV analysis step includes four parts: Read filtering, Read alignment, Read classification, and Read
counting.  
Input files of per-sample AAV analysis include sequencing reads (FASTQ), reference sequences (FASTA), vector config file
(TSV); and outputs include filtered sequencing reads (FASTQ), alignment file (BAM), read classification results (TXT),
and read count table (TXT).

#### Read filtering

Input reads are first filtered to ensure only complete and high-quality sequencing reads are used for later analysis.  
Filtering criteria include:  

__a. Read completeness.__
Complete sequencing reads are expected to include PCR primer sequences on both ends when PCR amplification is performed
for targeted sequencing. Sequencing reads without primer sequences on either end are filtered out due to truncation.  
Considering high sequencing error rate of long-read data, a 10-bp sliding window is used to check presence/absence of
primer sequence. The primer sequence is considered as present if any sequential 10 bases of the primer is present in the
first/last 100 bp of the sequencing read.
(e.g., for a 25-bp forward primer sequence, 16 sub-sequences of 10bp (sub-sequence from base 1-10, 2-11, 3-12, ..., and
16-25) will be used to check if the subsequence is present in first 100-bp of the sequencing read. If any one of these
16 subsequences is present, then this sequencing read is considered to include forward primer at beginning.)
Reverse complementary sequences of primer sequences will also be checked in case read is from anti-sense strand.
Sequencing reads containing only 1 side of primer sequence or neither side will be filtered out.  
If only forward or only reverse primer sequence is provided, then read filtering step will only check primer for one
side and skip the other side. If no primer sequence is provided (i.e. when no PCR is conducted, or user wants to disable
PCR primer-based filtering), read filtering step will skip checking for presence of PCR primers and all sequencing reads
will pass this filter.
  
__b. Base quality.__
Average base quality of a sequencing read is calculated based on the Phred quality score for all bases from the input
FASTQ file. Sequencing reads with average base quality below input "Minimal average base quality" (default varies by platform: PacBio HiFi = 30, Nanopore = 10)
are filtered out due to low quality.

#### Read alignment

Filtered sequencing reads are then aligned to the reference sequences using minimap2. The alignment preset is automatically selected based on the sequencing platform: "-x map-hifi" for PacBio HiFi data and "-x map-ont" for Oxford Nanopore data.
Read alignment file is then sorted and indexed using Samtools.

#### Read classification

![plot](./workflow_classification.png)

Read alignment results are used to classify reads into 10 categories. Overall, reads aligned to HDR/ITR sequences will
be classified based on the length of ITR and HDR sequences in the alignment; and reads aligned to WT will go through
variant calling process to be assigned as Unmodified or into a variant category. To achieve most accurate read
classification results, some additional processes were implanted, including re-aligning unmapped sequence, re-aligning
INS sequence, and re-aligning extra-short HDR reads to WT.

Major processes of read classification includes:  
__a. Process alignments on HDR/ITR sequences.__ For reads aligned to ITR sequences, if the alignment includes at
least 10bp on ITR sequences either on left or right side, this read will be assigned as ITR-integration or Non-HDR-with-ITR
(Truncated ITR-integration) depending on if the ITR integration length is above the input "Truncated HDR/ITR Threshold"
percentage (default is 95%) of expected full-length ITR-integration length. If the alignment includes less than 10bp of ITR
sequences on both side, it will be treated as aligned to HDR sequences as described in next paragraph.  
For reads aligned to HDR sequences, if it includes at least 10bp of transgene/HDR sequence, this read will be
assigned as HDR/Non-HDR-without-ITR depending on of integration length (cutoff is the same with ITR/Truncated-ITR cutoff).  
If the read alignment contains less than 10bp of HDR sequence, it will be re-aligned against only the WT sequence from
input reference and then re-processed as a read aligned to WT as described below in #c.  
  
__b. Process alignments on WT sequence.__ For reads aligned to wildtype (WT) sequence, if there are >100bp unaligned
sequences on either end of the alignment, the unaligned sequences will be extracted and re-aligned as described in #d.
Variant calling is performed on the read alignments to detect SNPs, indels (large and small), inversions, and
duplications within a configurable window (default 20bp) on both sides of the cleavage site.
If there is no variants reported, the read will be assigned as Unmodified.  
If only one type of variant is reported (e.g. only a 10bp DEL is reported for ReadA, two SNPs are reported for ReadB),
the read will be assigned to the corresponding variant category (e.g. ReadA -> DEL-small, ReadB -> SNP).  
If more than one type of variants are reported, the largest variant will be used to classify read category (e.g. ReadC
has a 90bp INS, two SNPs, and a 40bp DEL, ReadC -> INS).  
If a read carries insertion longer than 50bp, the insertion sequence will be extracted and realigned against HDR/ITR
sequences as described in #e.  
  
__c. Re-aligning extra-short HDR reads to WT.__ When a read is aligned to HDR/ITR sequences but the integrated HDR/ITR
sequence is shorter than 10bp, it will be too short to be assigned as Non-HDR-without-ITR, and the entire read will be
re-aligned against WT sequence to detect variants. The alignment against WT will then be processed as "aligned to WT"
as described in #b.  
  
__d. Re-aligning unmapped sequence.__ For reads aligned to WT with >100bp unmapped sequence (present as soft or hard
clip) on either end of the alignment, we first check all other supplementary alignments of the same read to make sure
these unmapped sequences are not mapped to any reference sequences. A new reference sequence file containing only HDR
and ITR sequences from original input reference is created and used to re-align these unmapped sequences using
minimap2 with platform-specific presets ("map-hifi" for PacBio, "map-ont" for Nanopore).
If the unmapped sequence can be aligned to HDR/ITR, the new alignment result will then be processed as "aligned to
HDR/ITR" as described in #a.
If the unmapped sequence still cannot be aligned to HDR/ITR, the read will be processed as "aligned to WT" as described
in #b.
  
__e. Re-aligning INS sequence.__ If a large (>50bp) insertion is reported for a read, re-aligning is performed to check
if inserted sequence is actually coming from AAV vector.
For reads containing large INS, the inserted sequence is extracted (e.g. a 1500bp INS is reported for ReadD, this 1500bp
sequence will be extracted, not the entire read) and re-aligned to the original input reference file.
If the inserted sequence is aligned to HDR/ITR sequence, this new alignment result will then be processed as "aligned
to HDR/ITR" as described in #a.
Otherwise, this read will be classified as INS-large.  

#### Read counting

Read counting step takes read classification result table (readname_$(SAMPLE_NAME).txt) as input and counts number of
reads in each read category.
A customized Python script is used to read the input table, count occurrence of each read category in second column, and
write a tab-delimited text file to list number of reads in each category.

### Merge classification table

Merge classification table merges the individual read count tables into a single table. A customized Python script is
used to read in all input read count tables, merge the tables according to column name, and write a tab-delimited text
file as output.  
If different vector config files are used in per-sample step (i.e. SampleA has AAV1 and SampleB has AAV2), merged table
will have more columns than individual count table to include extra columns for HDR/ITR of other AAVs. In this case, 0
will be assigned to these new columns for samples without this AAV (i.e. SampleA will have 0 reads for HDR-AAV2,
Non-HDR-with-ITR-AAV2, and Non-HDR-without-ITR-AAV2; same for SampleB in AAV1-related columns).

## Inputs

Inputs of the per-sample analysis workflow include five required inputs and optional inputs. Inputs of merge step
includes a list of read annotation results from per-sample step for all samples.

### Sequencing FASTQ

Long-read sequencing data of the sample from PacBio HiFi or Oxford Nanopore platforms. File can be in FASTQ or zipped FASTQ.GZ format. Currently only 1 FASTQ/FASTQ.GZ
file is allowed as input. If there are multiple FASTQ files for one sample, merge them into one FASTQ file before
running this pipeline.

### Reference sequences

A reference sequence file used for read alignment and read classification. There are three types of reference sequences:
WT, HDR and ITR, representing differen outcomes of AAV integration:

![plot](./example_ref.png)

- WT - Wildtype sequence covering the cleavage site. Sequence is usually identical with the human reference genome such
as GRCh38 or GRCh37.
- HDR - HDR knock in of AAV vector, containing trans gene sequence from AAV between the cleavage site.
- ITR - ITR integration of AAV vector. Containing ITR sequence and HA sequence on both side of the trans gene sequence.

This reference file should contain at least wildtype (WT) sequence and ITR sequence for each AAV vector. It is good to
include HDR sequence as well.  
  
An example reference sequence FASTA file containing sequences shown in above plot is like:

> \>ref_WT  
> ATCGATCG...  
> \>ref_HDR_AAV1  
> ATCGATCG...  
> \>ref_ITR_AAV1  
> ATCGATCG...  
> \>ref_ITR_AAV2  
> ATCGATCG...

### AAV vector config

A config file to list detailed information about AAV vectors. A header line must be present, and then each row
represents one sequence in the reference sequence file. 9 columns are required for the config file:

- AAV_Vector - AAV vector name.
- Ref_Name - Reference sequence name (should be identical with the sequence name from reference FASTA file).
- Ref_Type - Type of reference sequence. Choice from HDR, WT, ITR.
- ITR1_start - Starting position of first ITR element on left side.
- ITR1_end - Ending position of first ITR element on left side.
- gene_start - Starting position of trans gene.
- gene_end - Ending position of trans gene.
- ITR2_start - Ending position of second ITR element on right side.
- ITR2_end - Starting position of second ITR element on right side.

Hints:

- All positions in the config file should be 1-based.
- Use NA if any column is not applicable, i.e., use NA for ITR start/end for HDR.
- For WT sequence, use any existing AAV vector name in the AAV_Vector column (e.g. use AAV1 for WT sequence in the
example below). Do not use a new name (such as WT or NA), which will mislead the program to consider it as a new AAV
vector.
- Order of columns must follow the order as listed above (see example below). Order of rows does not matter.
  
An example vector file matched with previous reference sequence file is like:

\#AAV_Vector | Ref_Name | Ref_Type | ITR1_start | ITR1_end | gene_start | gene_end | ITR2_start | ITR2_end
:-|:-|:-|:-|:-|:-|:-|:-|:-
AAV1 | ref_WT | WT | NA | NA | 1000 | 1001 | NA | NA
AAV1 | ref_HDR_AAV1 | HDR | NA | NA | 1000 | 2000 | NA | NA
AAV1 | ref_ITR_AAV1 | ITR | 1000 | 1100 | 1300 | 2300 | 2500 | 2600
AAV2 | ref_ITR_AAV2 | ITR | 1000 | 1100 | 1300 | 3000 | 3200 | 3300

### Primer sequences

Primer sequence of the forward and reverse primers used in read filtering step. Both primer sequences should be in
5'->3' direction.  
If no primer sequence is provided as input, the filter step will not filter read based on presence of primers, and only
discards reads with low quality.

### Sample name

A string of Sample Name is used as prefix to name output files. Do not include space or special characters in the name.

### Optional inputs

- **Data Type** - Sequencing platform type. Choices: "pacbio-hifi" or "nanopore". This parameter determines the minimap2 alignment preset and platform-specific thresholds. Default value is "pacbio-hifi".
- **Minimal average base quality** - Cutoff used in read filtering step. Sequencing reads with average Phred base quality below this cutoff will be filtered out. Default based on data type: PacBio HiFi (30), Nanopore (10).
- **Truncated HDR/ITR Threshold** - Truncated HDR/ITR threshold. Input value should be a percentage between 0-1. Reads containing HDR/ITR integration with length shorter than this percentage of full length integration will be classified as Non-HDR with/without ITR. Default value is 0.99.
- **Variant Window** - Size of the window (in bp) around the cleavage site for variant detection. Default value is 20bp (±20bp from cleavage site).
- **Primer Check Length** - Length in base pairs from each end of the read to search for PCR primer sequences. Default value is 50bp.
- **Minimum ITR Length** - Minimum ITR sequence coverage required for ITR classification. Default value is 10bp.
- **5' Homology Arm Sequence** - 5' Homology Arm sequence in 5'->3' direction for HDR validation.
- **3' Homology Arm Sequence** - 3' Homology Arm sequence in 5'->3' direction for HDR validation.
- **HA Sequence Match Ratio** - Homology arm sequence identity threshold for HDR validation. Default based on data type: PacBio HiFi (0.98), Nanopore (0.90).
- **Left ITR Sequence** - Left ITR sequence in 5'->3' direction for ITR detection.
- **Right ITR Sequence** - Right ITR sequence in 5'->3' direction for ITR detection.
- **Seed Size** - Seed size for blastn alignment of ITR sequences. Default value is 15.
- **Percent Identity** - Minimum percent identity to match ITR sequences for blastn. Default value is 80.

### Inputs for merge classification table step

Inputs for merge classification table step are the single-sample read classification files for all samples.

- Read classification file - read classification output file from the per-sample workflow. Default file names are like
"$(SAMPLE_NAME)_read_classification.txt". Select this output file for all samples in the study so that a merged table
containing all samples will be generated for easier comparison and plot generation.

For detailed explanations of all output classification categories, see the [ALPINE Output Classification Categories](#alpine-output-classification-categories) section.

## Test Data and Examples

The `testdata/` folder contains example input and output files to help users understand the expected file formats and test the ALPINE pipeline:

### Example Input Files
- **`example_input.test_sample.fastq`** - Sample FASTQ file with long-read sequencing data
- **`example_ref.reference_trac.fa`** - Reference sequence file containing WT, HDR, and ITR sequences for TRAC locus
- **`example_ref.config_trac.txt`** - Vector configuration file specifying AAV vector details and genomic coordinates
- **`example_input.yaml`** - CWL workflow input parameter file showing how to configure pipeline runs

### Example Output Files
- **`example_output.readname_test_sample.txt`** - Per-read classification results with read names, categories, and variant sizes
- **`example_output.test_sample_read_classification.txt`** - Summary count table showing number of reads in each classification category

### Usage
These test data files serve multiple purposes:
- **Format templates**: Use as examples for preparing your own input files with correct formatting
- **Pipeline testing**: Run ALPINE with these files to verify proper installation and functionality
- **Parameter reference**: The YAML file demonstrates how to specify optional parameters for different experimental setups

To run ALPINE with the test data:
```bash
cwl-runner CWL/aav_per_sample_workflow.cwl testdata/example_input.yaml
```

## Outputs

### Per-sample outputs

#### Filtered FASTQ

Filtered FASTQ files are named as "$(SAMPLE_NAME).fastq.gz" by default. This is the output file from the read filtering
step, and it contains sequencing reads with primer sequences on both ends and with average base quality higher than a
threshold.

#### Sorted BAM

Sorted BAM files are named as "$(SAMPLE_NAME).bam" by default. This is the output file from sort alignment step. The
BAM file contains read alignments against input reference file, and is sorted based on coordinates and indexed. It
could be used for visualization using IGV or to check individual read alignments.

#### Read annotation

Read annotation files are named as "readname_$(SAMPLE_NAME).txt" by default. This is an intemediate output file from
read classification step.  
There is no header line for this file. It contains three columns:

- Read name - Name of the read.
- Read category - The category a read is classified into. See details in Steps - Read classification section.
- Size - This is the size of the variant or the integration length of HDR/ITR. No value is provided for Unmodified
category.  

An example is like:

sequence_read_1 | DEL-small | 5
:-|:-|:-
__sequence_read_2__ | __HDR-AAV1__ | __1000__
__sequence_read_3__ | __Non-HDR-without-ITR-AAV2__ | __1200__
  
ALPINE classifies reads into multiple categories including unmodified reads, CRISPR-induced variants (deletions, insertions, SNPs, inversions, duplications), AAV integration events (HDR, Non-HDR with/without ITR), and other integration outcomes. See the ALPINE Output Classification Categories section below for detailed descriptions.

#### Single-sample classification table

Single-sample classification table files are named as "$(SAMPLE_NAME)_read_classification.txt" by default. This table
includes number of reads classified into each category.

An example table is like:

\#Sample | DEL-large | DEL-small | INS-large | INS-small | SNP | INV | DUP | Unmodified | HDR-AAV1 | Non-HDR-with-ITR-AAV1 | Non-HDR-without-ITR-AAV1 | HDR-AAV2 | Non-HDR-with-ITR-AAV2 |  Non-HDR-without-ITR-AAV2 
:-|:-|:-|:-|:-|:-|:-|:-|:-|:-|:-|:-|:-|:-|:-
SampleA | 10 | 1 | 0 | 2 | 3 | 7 | 5 | 7 | 1000 | 10 | 7  | 0 | 0 | 1 

### Merged classification table

Merged classification table file is the merged table of single-sample classification tables. It contains the number of
reads in each class for all samples.
  
An example table is like:

|Sample|DEL-large|DEL-small|INS-large|INS-small|SNP|INV|DUP|Unmodified|HDR-AAV1|Non-HDR-with-ITR-AAV1|Non-HDR-without-ITR-AAV1|HDR-AAV2|Non-HDR-with-ITR-AAV2|Non-HDR-without-ITR-AAV2|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
SampleA|10|1|0|2|3|7|5|7|1000|10|8|0|0|1
SampleB|2|4|3|1|5|8|4|1|0|0|1000|20|10|6

### ALPINE Output Classification Categories

ALPINE classifies each sequencing read into one of the following categories based on the detected genomic alterations (numerical thresholds mentioned below represent default settings and are configurable parameters):

#### Unmodified Outcomes
- **Unmodified**: Reads that align perfectly to the wildtype reference sequence with no detectable variants within the variant calling window. These represent cells that were not edited by the CRISPR system.
- **Unmodified-with-SNP**: Reads that contain only single nucleotide polymorphisms (SNPs) or substitutions compared to the reference. These likely represent natural genetic variation or sequencing errors rather than CRISPR-induced modifications.

#### Transgene Integration Events
- **HDR**: Homology-directed repair events where the AAV vector template has been successfully integrated with high fidelity. These reads show coverage of HDR template sequence with high sequence identity to the expected template, with minimal (<10bp) ITR sequence content. These represent the desired precise gene editing outcome.
- **Non-HDR-with-ITR**: Vector integration events that contain AAV inverted terminal repeat (ITR) sequences (>10bp coverage) but lack proper HDR template integration. These represent partial or imprecise vector integration events.
- **Non-HDR-without-ITR**: Vector integration events that show evidence of AAV vector sequences without ITR content, but lack proper HDR template integration (<98% identity) or have truncated HDR integration (<99% transgene length). These represent imprecise integration events.

#### CRISPR-Induced Variants
- **DEL-small**: Small deletions (<50bp) at the target site, typically resulting from imprecise non-homologous end joining (NHEJ) repair of CRISPR-induced double strand breaks.
- **DEL-large**: Large deletions (≥50bp) at the target site, representing more extensive NHEJ-mediated repair or loss of genomic material between multiple cut sites.
- **INS-small**: Small insertions (<50bp) at the target site, usually resulting from NHEJ repair mechanisms that insert random nucleotides during double strand break repair.
- **INS-large**: Large insertions (≥50bp) at the target site that do not match the provided AAV vector sequence. These may result from integration of contaminating DNA, complex rearrangements, or large NHEJ-mediated insertions.
- **INV**: Inversions detected at the target site, where a genomic segment has been reversed in orientation, typically resulting from NHEJ repair between two CRISPR cut sites.
- **DUP**: Duplications detected at the target site, where genomic segments have been duplicated, often arising from complex NHEJ repair mechanisms.

#### Unresolved Classifications
- **Unclassified**: Reads that could not be confidently classified into any of the above categories. These include unmapped reads, reads with complex alignment patterns, and reads not spanning the cleavage site.

## How to run pipeline

### On SevenBridges

After creating a project on SBG, go to the Apps page of the reference project
and copy workflow "Transgene Integration Long-Read Launcher" into your working project.
Provide sample sheet and other inputs to run launcher app. Launcher will run per-sample workflow for each sample,
and run merge step to generate  merged read classification table.

### On Arvados

After creating a project on SBG, go to reference project and run "Transgene Integration Long-Read Launcher" workflow.
Specify the working project as the project where workflow will run and provide sample sheet
and other inputs. Launcher will run per-sample workflow for each sample, and run merge step to generate
merged read classification table.

## Update workflows in SBG and Arvados

__(For Admin only)__ Makefile is used in this repo to simplify workflow/launcher updates on SBG and Arvados after
modifications have been made to CWL/Docker files.
To push updated workflows:

```sh
# Clone the latest repo
git clone https://github.com/Maggi-Chen/ALPINE
cd ALPINE
# Push updated workflows to SBG
make sevenbridges-push-workflow
# Push updated workflows to Arvados
make arvados-push-workflows
```

To push updated launcher code:

```sh
# Push updated launcher to SBG
make sevenbridges-push-launcher
# Push updated launcher to Arvados
make arvados-push-launcher
```

## Contacts

For questions about this pipeline, please contact:

- [Maggie Chen](mailto:maggiustc@gmail.com?subect=%5BPangenome%20Pipeline%20Question%5D)
- [Xing Huang Gao](mailto:xinghuang.gao@bms.com?subect=%5BPangenome%20Pipeline%20Question%5D)
