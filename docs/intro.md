## Introduction

Pore-C is an end-to-end workflow unique to Oxford Nanopore which combines chromatin conformation capture (3C) with direct, long nanopore sequencing reads. With nanopore reads, long-range, multi-way contact information can be obtained.  Pore-C can be used to scaffold and improve the quality of existing genome assemblies to generate chromosome-scale assemblies, discover unique insights into the higher-order genome organisation of a species of interest, and reveal epigenetic information for a more comprehensive understanding of gene regulation. Find out more about the workflow [here](https://nanoporetech.com/about-us/news/pore-c-complete-end-end-workflow-chromatin-conformation-capture-published-nature)

This nextflow workflow will create virtual digests of the genome using the [pore-c-py package](https://github.com/epi2me-labs/pore-c-py), align the resulting monomers against a reference genome with [minimap2](https://github.com/lh3/minimap2) and then filter spurious alignments, detect ligation junctions and assign fragments. The resulting BAM alignment file can be used with downstream tools. 

Optionally the workflow can output a [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) report which uses [pairtools](https://github.com/open2c/pairtools).
