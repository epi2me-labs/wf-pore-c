### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bam | string | An unaligned BAM file containing Pore-C concatemer sequences | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| sample_sheet | string | A CSV file used to map barcodes to sample aliases and optionally provide per-sample parameters. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Optionally, a `cutter` column can contain the name of the enzyme used per sample (see the `--cutter` parameter for more details) and a `vcf` column can be used to provide a phased VCF file per sample if you require haplotagged alignments. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| ref | string | A fasta file containing the reference genome to map against |  |  |
| vcf | string | An optional phased VCF file that will be used to haplotag alignments |  |  |
| cutter | string | The enzyme used in the restriction digest. | Any enzyme from the Biopython restriction dictionary can be used. See `https://github.com/biopython/biopython/blob/master/Bio/Restriction/Restriction_Dictionary.py` | NlaIII |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |
| hi_c | boolean | Output a Hi-C formatted file; will convert pairs format to a Hi-C (`.hic`) file which will be compatible with [juicer](https://github.com/aidenlab/juicer) | Load this file with [Juice box](https://www.aidenlab.org/juicebox/) for an alternative contact map visualisation. | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| chunk_size | integer | Process input in chunks of this number of reads. | To reduce per-process memory requirements for large datasets, process the inputs in chunks of reads. Set to 0 to process entire dataset in one go. | 10000 |
| threads | integer | Set maximum number of threads to use for more intense processes (limited by config executor cpus). We recommend a minimum of 4, but if available 19. |  | 4 |


### Pore-C Tools Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| minimap2_settings | string | The minimap2 settings for mapping monomers |  | -x map-ont |
| coverage | boolean | Calculate restriction-fragment coverage using mosdepth |  | False |
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
| mcool | boolean | Create a multi-resolution cooler file. Will output the cooler formatted file which you can load with cooler. | see 'https://open2c.github.io/cooler' for more details. | False |
| cool_bin_size | integer | The bin size of the cooler output file in base pairs. | see 'https://open2c.github.io/cooler' for more details. | 1000 |
| mcool_resolutions | string | The resolutions of the mcool file in pixels (see cooler documentation for details). | Comma-separated list of target resolutions. Use suffixes B or N to specify a progression: B for binary (geometric steps of factor 2), N for nice (geometric steps of factor 10 interleaved with steps of 2 and 5). This is the equivalent of the `--resolutions` flag in cooler; see an example here 'https://cooler.readthedocs.io/en/latest/cli.html'. | 1000,2000,5000N |


### Paired-end BAM Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| paired_end | boolean | Create mock paired-end BAM files. |  | False |
| filter_pairs | boolean | Filter paired end reads using minimum and maximum distance parameters. |  | False |
| paired_end_minimum_distance | integer | Remove trans pairs and cis- pairs separated by a distance shorter than this |  | -1 |
| paired_end_maximum_distance | integer | Remove trans pairs and cis- pairs separated by a distance greater than this |  | -1 |


### Misc

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |


