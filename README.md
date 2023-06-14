# Workflow Pore-C

This repository contains a [nextflow](https://www.nextflow.io/) workflow
to convert data generated with the [Pore-C](https://nanoporetech.com/resource-centre/porec) protocol to outputs in a variety of standard file formats used by downstream tools.




## Introduction

Pore-C is an end-to-end workflow unique to Oxford Nanopore which combines chromatin conformation capture (3C) with direct, long nanopore sequencing reads. With nanopore reads, long-range, multi-way contact information can be obtained.  Pore-C can be used to scaffold and improve the quality of existing genome assemblies to generate chromosome-scale assemblies, discover unique insights into the higher-order genome organisation of a species of interest, and reveal epigenetic information for a more comprehensive understanding of gene regulation. Find out more about the workflow [here](https://nanoporetech.com/about-us/news/pore-c-complete-end-end-workflow-chromatin-conformation-capture-published-nature)

This nextflow workflow will create virtual digests of the genome using the [pore-c-py package](https://github.com/epi2me-labs/pore-c-py), align the resulting monomers against a reference genome with [minimap2](https://github.com/lh3/minimap2) and then filter spurious alignments, detect ligation junctions and assign fragments. The resulting BAM alignment file can be used with downstream tools. 

Optionally the workflow can output a [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) report which uses [pairtools](https://github.com/open2c/pairtools).




## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
to provide isolation of the required software. Both methods are automated
out-of-the-box provided either docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-pore-c --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* Coordinate sorted Bam and index
* Name sorted Bam

Alignments in the Bams are annotated/tagged with a "walk", which simply
enumerates the alignment coordinates of the monomers comprising the concatemer.
  
Optional outputs:
* Pairs file `.pairs.gz` containing the contact information in a tabular format, which is human readable and can be used with downstream tools. See [Pairtools documentation](https://pairtools.readthedocs.io/en/latest/formats.html#pairs) for full specification.
* Pairs `.stats` file with summary statistics. See this [overview](https://pairtools.readthedocs.io/en/latest/stats.html) for a full specification.
* Pairs html report with result including an interactive contact map and statistics. See [pairsqc documentation](https://github.com/4dn-dcic/pairsqc) for further details.
* Multi-resolution cool `.mcool` file which can be used with downstream tools to provide a high resolution genomic interaction matrix. See [Cool tools documentation](https://github.com/open2c/cooltools) for details on downstream analysis. 
* Chromunity directory with parquet files which can be used with the Chromunity package. Chromnity enables the nomination and statistical evaluation of high order interactions. See [Chromunity documentation](http://mskilab.com/chromunity/tutorial.html) for further details.
* `Fragments.bed` file with the DNA fragments created from the virtual digest.





## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [pore-c-py](https://github.com/epi2me-labs/pore-c-py)
* [minimap2](https://github.com/lh3/minimap2)
* [pairtools](https://github.com/open2c/pairtools)
* [Pore-C](https://nanoporetech.com/resource-centre/porec): multi-contact, chromosome conformation capture for both genome-wide and targeted analyses.
* Visit the [Nanopore Community](https://community.nanoporetech.com/info_sheets/restriction-enzyme-pore-c/v/rpc_s1015_v1_revf_12nov2019) for guidance through every step of the Pore-C workflow. Optimised Nanopore workflows have been developed: blood samples, [plants](https://nanoporetech.com/sites/default/files/s3/literature/plant-pore-c-workflow.pdf), and cell culture, animal, or insect tissues – the method is simple and scalable.


### References

Deshpande, A.S., Ulahannan, N., Pendleton, M. et al. Identifying synergistic high-order 3D chromatin conformations from genome-scale nanopore concatemer sequencing. Nat Biotechnol (2022). https://doi.org/10.1038/s41587-022-01289-z: https://www.nature.com/articles/s41587-022-01289-z  