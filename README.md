<p align="center">
    <img height="150" src="https://uni-muenster.sciebo.de/s/I0dmJFFdq4UUYmk/download">
</p>

`wg-blimp` (Whole Genome BisuLfIte sequencing Methylation analysis Pipeline) can be utilised to analyse WGBS data. It performs alignment, qc, methylation calling, DMR calling, segmentation and annotation using a multitude of tools. First time using `wg-blimp`? We recommend having a look at our [step-by-step guide](https://github.com/MarWoes/wg-blimp/wiki/Tutorial).

## Requirements
To run `wg-blimp` you need a UNIX environment that contains a [Bioconda](http://bioconda.github.io/) setup.

## Installation

### Recommended
It is advised to install `wg-blimp` through Bioconda. It is also recommended to install `wg-blimp` in a fresh environment, as it has many dependencies that may conflict with other packages, for this you can use:

```
conda create -n wg-blimp wg-blimp python=3.6.7 r-base=3.6.1
```

### From source
You can also install `wg-blimp` from source using
```
python setup.py install
```
Using this installation method requires you to make sure all external tools are installed (such as bwa-meth).

## Running wg-blimp

### WGBS pipeline

`wg-blimp` is a cli wrapper for the WGBS pipeline implemented using [Snakemake](http://snakemake.readthedocs.io/). In general, a pipeline config us fed to the Snakemake workflow and the corresponding tools are called. However, `wg-blimp` also provides some commands to ease creation of config files, or working without config files altogether.

The command `wg-blimp run-snakemake` will run the pipeline with its default parameters. Make sure to set the `--cores` and `--genome-build` options appropriately. This command will also internally create a `config.yaml` file containing all parameters used for the analysis.

However, in case the default configurations are not sufficient, users can provide their own configurations. The commands `wg-blimp create-config` and `wg-blimp run-snakemake-from-config` can be used for this purpose.

`wg-blimp` will attempt to match .fastq files to sample names by searching for sample names in .fastq file names. By default Illumina naming conventions are expected, e.g. for a samples _test1_ the .fastq files should be named as follows:
```
test1_L001_R1_001.fastq.gz
test1_L001_R1_002.fastq.gz
test1_L001_R2_001.fastq.gz
test1_L001_R2_002.fastq.gz
test1_L002_R1_001.fastq.gz
test1_L002_R1_002.fastq.gz
test1_L002_R2_001.fastq.gz
test1_L002_R2_002.fastq.gz
```

If names derive from this pattern, users can adjust the regular expression to match in the config file's `rawsuffixregex` entry.

The folder structure created by `wg-blimp run-snakemake` will look as follows:

* alignment - contains all bam/bai files
* dmr - contains dmr files by different callers
* logs - each pipeline step deposits its logs here
* methylation - methylation bedgraph files
* qc - multiqc and other qc related files
* raw - text files describing which fastq files have been used for each sample
* segmentation - methylome segments (UMRs/LMRs/PMDs) as computed by MethylSeekR
* config.yaml - configuration file used for the analysis

It is recommended to check the *raw* folder if all samples contain the correct raw fastq source files.
When in doubt, `wg-blimp` also allows for explicit association of samples and read files by setting `sample_fastq_csv` in the configuration file.
An example csv file could look as follows (column names must be set to `sample`, `forward` and `reverse`):
```
sample,forward,reverse
sample1,/my/path/sample1_L1_1.fq.gz,/my/path/sample1_L1_2.fq.gz
sample1,/my/path/sample1_L2_1.fq.gz,/my/path/sample1_L2_2.fq.gz
sample2,/my/path/sample2_L1_1.fq.gz,/my/path/sample2_L1_2.fq.gz
sample3,/my/path/sample3_L1_1.fq.gz,/my/path/sample3_L1_2.fq.gz
```

### Shiny GUI

You can use the command `wg-blimp run-shiny` to load one or more project config files into a shiny GUI for easier access.

## Example

Some example `.fastq`can be found on [Sciebo](https://uni-muenster.sciebo.de/s/7vpqRSEATYcvlnP). You can use the command
```
wg-blimp run-snakemake <folder-with-fastqs> <reference.fa> simulated1,simulated2 simulated3,simulated4 <output-folder>
```

Please not that the pipeline commands also allow a `--use-sample-files` option so sample groups can be loaded from text files instead of comma separates files.


## Config parameters

The following entries are used for running the Snakemake pipeline and may be specified in the `config.yaml` files:

| Key | Value |
| --- | ----- |
| *annotation_allowed_biotypes* | Only genes with this biotype will be annotated in the DMR table |
| *annotation_min_mapq* | When annotating coverage, only use reads with a minimum mapping quality |
| *bsseq_local_correct* | Use local correction for bsseq DMR calling. Usually, setting this to FALSE will increase the number of calls. |
| *cgi_annotation_file* | Gzipped csv file used for cg island annotation. |
| *computing_threads* | Number of processors a single job is allowed to use. Remember to use `--cores` parameter for Snakemake. |
| *dmr_tools* | Tools to use for DMR calling. Available: `bsseq`, `camel`, `metilene`
| *gene_annotation_file* | File used for genetic annotation. |
| *group1* | Samples in first group for DMR analysis |
| *group2* | Samples in second group for DMR analysis |
| *io_threads* | IO intensive tools virtually reserve this many cores (while actually using only one) to reduce file system IO load. |
| *methylation_rate_on_chromosomes* | Compute methylation rates for these chromosome during qc |
| *methylseekr_cgi_genome* | Reference genome to use for MethylSeekR CGI queries. |
| *methylseekr_fdr_cutoff* | FDR cutoff for MethylSeekR segmentation. |
| *methylseekr_methylation_cutoff* | Methylation cutoff for MethylSeekR segmentation. |
| *methylseekr_pmd_chromosome* | Chromosome to compute MethylSeekR alpha values for. |
| *min_cov* | Minimum average coverage for methylation calling |
| *min_cpg* | Minimum number of CpGs in a DMR to be called |
| *min_diff* | Minimum average difference between the two groups for DMR calling |
| *output_dir* | Directory containing all files created by the pipeline |
| *promoter_tss_distances* | Distance interval around TSS's to be recognized as promoters in DMR annotation. |
| *qualimap_memory_gb* | Gigabytes of RAM to allocate for Qualimap. If samples are too large, this must be increased to prevent crashes. |
| *rawdir* | Directory containing .fastq files |
| *rawsuffixregex* | The regular expressions to match for paired reads. By default, Illumina naming conventions are accepted. |
| *ref* | .fasta reference file |
| *repeat_masker_annotation_file* | File containing repeat masker annotation |
| *repeat_masker_links* | Repeat masker files are relatively big and are only downloaded on demand from the links specified here. |
| *sample_fastq_csv* | Optional CSV file containing association between samples and read files. The CSV must contain a header with column names `sample`, `forward` and `reverse`. When this option is set, parameters *rawdir* and *rawsuffixregex* are ignored. |
| *samples* | All samples (usually concatenation of group1 and group2) |
| *target_files* | Files to be generated by the Snakemake workflow |
| *transcript_start_site_file* | File to read transcription start sites from (for DMR annotation) |

## Reporting errors / Requesting features
If anything goes wrong using `wg-blimp` or any features are missing, feel free to open an issue or to contact Marius Wöste ( mar.w@wwu.de )
