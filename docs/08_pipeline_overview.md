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