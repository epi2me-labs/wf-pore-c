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