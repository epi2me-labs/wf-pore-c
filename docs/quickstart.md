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

