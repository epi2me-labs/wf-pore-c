Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-template-report.html | Report for all samples. | aggregated |
| Per file read stats | ingress_results/reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Per read stats | ingress_results/reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | ingress_results/reads/fastcat_stats/run_ids | List of run ID's present in reads. | aggregated |
| Meta map json | ingress_results/reads/metamap.json | Meta data used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | ingress_results/reads/{{ alias }}.fastq.gz | Per-sample reads concatenated in to one fastq file. | per-sample |
| Coordinate-sorted Bam | bams/{{ alias }}.cs.bam | Coordinate-sorted Bam. | per-sample |
| Coordinate-sorted Bam Index | bams/{{ alias }}.cs.bam.bai | Coordinate-sorted Bam Index. | per-sample |
| Name-sorted Bam | bams/{{ alias }}.ns.bam | Name-sorted Bam. | per-sample |
| Pairs file | pairs/{{ alias }}.pairs.gz | This file contains contact information in a human-readable tabular format, and can be used with downstream tools. See [Pairtools documentation](https://pairtools.readthedocs.io/en/latest/formats.html#pairs) for full specification. | per-sample |
| Pairs summary stats file | pairs/{{ alias }}.pairs.stats.txt | Summary statistics of the pairs file. See this [overview](https://pairtools.readthedocs.io/en/latest/stats.html) for a full specification. | per-sample |
| Pairs summary report | pairs/{{ alias }}.pairs.stats.html | Pairs html report with result including an interactive contact map and statistics. See [pairsqc documentation](https://github.com/4dn-dcic/pairsqc) for further details. | per-sample |
| Multi-resolution cool file | cooler/{{ alias }}.mcool | Multi-resolution cool `.mcool` file which can be used with downstream tools to provide a high resolution genomic interaction matrix. See [Cool tools documentation](https://github.com/open2c/cooltools) for details on downstream analysis. | per-sample |
| Paired-end BAM | paired_end/{{ alias }}.ns.bam | Mock paired end BAM. | per-sample |
| Chromunity parquet files. | chromunity | Chromunity directory with parquet files which can be used with the Chromunity package. Chromunity enables the nomination and statistical evaluation of high order interactions. See [Chromunity documentation](http://mskilab.com/chromunity/tutorial.html) for further details. | per-sample |
| Fragments BED | paireds/fragments.bed | File with the DNA fragments created from the virtual digest. | per-sample |
| Hi-C for contact map | hi-c/{{ alias }}.hic | File which can be loaded into the [Juice box tool](https://www.aidenlab.org/juicebox/) for an alternative contact map visualisation. | per-sample |
| Filtered out reads | filtered_out/{{ alias }}.bam | BAM file containing any reads that were filtered out at the digest step and not included in the analysis. | per-sample |
