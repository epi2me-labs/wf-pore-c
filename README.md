# Pore-c Workflow

Workflow for analysing Pore-c data for chromatin conformation capture.



## Introduction

Pore-C is an end-to-end workflow unique to Oxford Nanopore which combines chromatin conformation capture (3C) with direct, long nanopore sequencing reads. With nanopore reads, long-range, multi-way contact information can be obtained. 

This workflow can be used for the following:

* Pre-processing a reference genome or draft assembly to generate auxiliary files used in downstream analyses.
* Creating virtual digests of Pore-c reads.
* Filtering the raw reads to remove any that might break downstream tools.
* Align virtually digested reads against a reference genome.
* Processing results to filter spurious alignments, detect ligation junctions and assign fragments.
* Outputting aligned, sorted and annotated BAM files.
* Generating a contact map, which shows the intensity of the physical interaction between two genome regions.
* Create output files for downstream analysis in the following formats.
  - [Pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)
  - [Cooler format](https://mirnylab.github.io/cooler/)
  - [Hic format](https://github.com/aidenlab/juicer/wiki/)



## Compute requirements

Recommended requirements:

+ CPUs = 64
+ Memory = 128GB

Minimum requirements:

+ CPUs = 8
+ Memory = 32GB

Approximate run time: 12 hours for 100GB input BAM using the recommended resources, this will vary depending on number of monomers found per read.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-pore-c --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-pore-c
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-pore-c/wf-pore-c-demo.tar.gz
tar -xzvf wf-pore-c-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-pore-c \
	--bam 'wf-pore-c-demo/porec_test.concatemers.bam' \
	--chunk_size 100 \
	--cutter 'NlaIII' \
	--hi_c \
	--mcool \
	--paired_end \
	--paired_end_maximum_distance 200 \
	--paired_end_minimum_distance 100 \
	--phased_vcf 'wf-pore-c-demo/porec_test.phased_variants.vcf.gz' \
	--ref 'wf-pore-c-demo/porec_test.fasta' \
	--vcf 'wf-pore-c-demo/porec_test.phased_variants.vcf.gz' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example


<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or unaligned BAM files as input.

The FASTQ and BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bam | string | An unaligned BAM file containing Pore-C concatemer sequences. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| sample_sheet | string | A CSV file used to map barcodes to sample aliases and optionally provide per-sample parameters. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Optionally, a `cutter` column can contain the name of the enzyme used per sample (see the `--cutter` parameter for more details) and a `vcf` column can be used to provide a phased VCF file per sample if you require haplotagged alignments. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| ref | string | A FASTA file containing the reference genome to map against. |  |  |
| vcf | string | An optional phased VCF file that will be used to haplotag alignments. |  |  |
| cutter | string | The enzyme used in the restriction digest. | Any enzyme from the Biopython restriction dictionary can be used. See `https://github.com/biopython/biopython/blob/master/Bio/Restriction/Restriction_Dictionary.py`. This can also be defined per sample: see `--sample_sheet` parameter. | NlaIII |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |
| hi_c | boolean | Output a Hi-C formatted file; will convert pairs format to a Hi-C (`.hic`) file which will be compatible with [juicer](https://github.com/aidenlab/juicer). | Load this file with [Juice box](https://www.aidenlab.org/juicebox/) for an alternative contact map visualisation. | False |
| bed | boolean | Output a BED file of the paired-end BAM alignments for use with downstream tools. Setting this to true will also trigger creation of the paired-end BAM. | Will use the paired-end BAM to create a BED file compatible with downstream tools including scaffolding tool [Yahs](https://github.com/c-zhou/yahs). | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| chunk_size | integer | Process input in chunks of this number of reads. | To reduce per-process memory requirements for large datasets, process the inputs in chunks of reads. Set to 0 to process entire dataset in one go. | 20000 |
| threads | integer | Set maximum number of threads to use for more intense processes (limited by config executor cpus). We recommend a minimum of 4, but if available 19. |  | 4 |


### Pore-C Tools Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| minimap2_settings | string | The minimap2 settings for mapping monomers. |  | -x map-ont |
| max_monomers | integer | The maximum number of monomers allowed for a read to be included in downstream analysis. |  | 250 |
| coverage | boolean | Calculate restriction-fragment coverage using mosdepth. |  | False |
| summary_json | boolean | Output pore-c-py annotation summary in json format. |  | True |


### Chromunity Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| chromunity | boolean | Create parquet files for Chromunity. | See the chromunity documentation for further details 'https://github.com/mskilab/chromunity'. | False |
| chromunity_merge_distance | integer | Merge colinear alignments separated by less than this base pair distance into a single monomer. |  | -1 |


### 4DN files Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| pairs | boolean | Create a 4DN-format pairs file (also calculate stats). | Outputs a directory with a pairs stats report and a pairs file which can be used for downstream anaylsis. | False |
| pairtools_chunksize | integer | Number of pairs to be processed in each chunk in the prepare_hic process which uses the pairtools dedup tool. | Reduce for lower memory footprint. Below 10,000 performance starts suffering significantly. | 100000 |
| mcool | boolean | Create a multi-resolution cooler file. Will output the cooler formatted file which you can load with cooler. | See 'https://open2c.github.io/cooler' for more details. | False |
| cool_bin_size | integer | The bin size of the cooler output file in base pairs. | See 'https://open2c.github.io/cooler' for more details. | 1000 |
| mcool_resolutions | string | The resolutions of the mcool file in pixels (see cooler documentation for details). | Comma-separated list of target resolutions. Use suffixes B or N to specify a progression: B for binary (geometric steps of factor 2), N for nice (geometric steps of factor 10 interleaved with steps of 2 and 5). This is the equivalent of the `--resolutions` flag in cooler; see an example here 'https://cooler.readthedocs.io/en/latest/cli.html'. | 1000,2000,5000N |


### Paired-end BAM Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| paired_end | boolean | Create mock paired-end BAM files. |  | False |
| filter_pairs | boolean | Filter paired-end reads using minimum and maximum distance parameters. |  | False |
| paired_end_minimum_distance | integer | Remove trans/cis pairs separated by a distance shorter than this. |  | -1 |
| paired_end_maximum_distance | integer | Remove trans/cis pairs separated by a distance greater than this. |  | -1 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-template-report.html | Report for all samples. | aggregated |
| Per file read stats | ./ingress_results/reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Per read stats | ./ingress_results/reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | ./ingress_results/reads/fastcat_stats/run_ids | List of run ID's present in reads. | aggregated |
| Meta map json | ./ingress_results/reads/metamap.json | Meta data used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | ./ingress_results/reads/{{ alias }}.fastq.gz | Per-sample reads concatenated in to one fastq file. | per-sample |
| Coordinate-sorted Bam | ./bams/{{ alias }}.cs.bam | Coordinate-sorted Bam. | per-sample |
| Coordinate-sorted Bam Index | ./bams/{{ alias }}.cs.bam.bai | Coordinate-sorted Bam Index. | per-sample |
| Name-sorted Bam | ./bams/{{ alias }}.ns.bam | Name-sorted Bam. | per-sample |
| Pairs file | ./pairs/{{ alias }}.pairs.gz | This file contains contact information in a human-readable tabular format, and can be used with downstream tools. See [Pairtools documentation](https://pairtools.readthedocs.io/en/latest/formats.html#pairs) for full specification. | per-sample |
| Pairs summary stats file | ./pairs/{{ alias }}.pairs.stats.txt | Summary statistics of the pairs file. See this [overview](https://pairtools.readthedocs.io/en/latest/stats.html) for a full specification. | per-sample |
| Pairs summary report | ./pairs/{{ alias }}.pairs.stats.html | Pairs html report with result including an interactive contact map and statistics. See [pairsqc documentation](https://github.com/4dn-dcic/pairsqc) for further details. | per-sample |
| Multi-resolution cool file | ./cooler/{{ alias }}.mcool | Multi-resolution cool `.mcool` file which can be used with downstream tools to provide a high resolution genomic interaction matrix. See [Cool tools documentation](https://github.com/open2c/cooltools) for details on downstream analysis. | per-sample |
| Paired-end BAM | ./paired_end/{{ alias }}.ns.bam | Mock paired end BAM. | per-sample |
| Chromunity parquet files. | ./chromunity | Chromunity directory with parquet files which can be used with the Chromunity package. Chromunity enables the nomination and statistical evaluation of high order interactions. See [Chromunity documentation](http://mskilab.com/chromunity/tutorial.html) for further details. | per-sample |
| Fragments BED | ./paireds/fragments.bed | File with the DNA fragments created from the virtual digest. | per-sample |
| Hi-C for contact map | ./hi-c/{{ alias }}.hic | File which can be loaded into the [Juice box tool](https://www.aidenlab.org/juicebox/) for an alternative contact map visualisation. | per-sample |
| Filtered out reads | ./filtered_out/{{ alias }}.bam | BAM file containing any reads that were filtered out at the digest step and not included in the analysis. | per-sample |




## Pipeline overview

### 1. Concatenate input files and generate per read stats.

This workflow accepts FASTQ or unaligned BAM as input. [Fastcat or Bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2. Index reference

The input reference genome is indexed with [Minimap2](https://github.com/lh3/minimap2).

### 3. Split input file

The reads are indexed in chunks for parallel processing using the `chunk_size` parameter which is defaulted to 10,000.

### 4. Digest Reads

Chimeric Pore-C reads are digested using the [Pore-c-py](https://github.com/epi2me-labs/pore-c-py) python package. The enzyme provided to the `cutter` parameter will be used by the Pore-c-py package to find the corresponding sequence using the [Biopython](https://biopython.org/) restriction enzymes library. Any reads containing more than `max_monomers` (default: 250) will be excluded at this stage as they are assumed to have been created in error.

### 5. Align Reads

The monomers are then aligned with the reference genome using Minimap2.

### 6. Annotate

The Pore-c-py package is then used again to filter spurious alignments, detect ligation junctions and assign chimeric fragments. The aligned segmnets will be used to generate a "walk" which enumerates the alignment coordinates of the monomers comprising the chimeric read and this is used to annotate the alignments.

### 7. Output BAMS

The Pore-c-py will output the tagged alignments in a name sorted and coordinate sorted BAM. If the `paired_end` parameter is selected a mock paired end bam will also be output, this is for use with downstream tools such as [Pairtools](https://github.com/open2c/pairtools). At this stage if the [Chromunity](https://github.com/mskilab-org/chromunity) parameter is set to true the annotate script will also output the parquet files required for us with the downstream Chromunity tool.

### 8. Haplotag Alignments

If a phased VCF is provided using the `vcf` parameter the output BAM will be haplotagged using [Whatshap](https://github.com/whatshap/whatshap).

### 9. Merge BAMS

The outputs BAM's from each of the split chunks will be merged and sorted per sample using [Samtools](https://www.htslib.org/doc/samtools-merge.html).

### 10. Coverage is calculated

If the `coverage` parameter is set to true [Mosdepth](https://github.com/brentp/mosdepth) is used to find coverage across the input reference genome.

### 11. Additional output formats for downstream analysis

The workflow will output several formats that can be used with downstream tools. 

+ [Pairtools](https://github.com/open2c/pairtools) is used to create pairs format file and html report which contains a contact map and other statistics. Use the `pairs` parameter for the workflow to generate this output.

+ [Cooler](https://github.com/open2c/cooler) is used to output cooler format for use with cooler, a multi-resolution contact map. Use the `mcool` parameter to generate this output.

+ [Juicer](https://github.com/aidenlab/juicer) tools is used to create `.hic` format file which can be used for visualising the file which can be loaded into the [Juice box tool](https://www.aidenlab.org/juicebox/) for an alternative contact map visualisation. Use the `hi_c` parameter to generate this output. 



## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ If the workflow breaks with a memory error, try running the workflow again with a reduced chunk size parameter.



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->
* Does the workflow have support for a scaffolding tool? * - Currently we do not support any scaffolding tool but you may like to try [Yahs](https://academic.oup.com/bioinformatics/article/39/1/btac808/6917071).

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-pore-c/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



