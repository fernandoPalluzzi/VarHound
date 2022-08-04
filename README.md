# VarHound (VH)

VH is a suite for genomic variant storage, retrieval, and validation through coverage analysis. It includes Python/R tools
for the storage, management, and retrieval of large VCF files. In addition, VH offers a set of sample-level and region-level
validation tools, through coverage analysis. VH is originally designed for truSight oncology 500 (TSO500) assays, but 
can be applied to any other DNA/RNA sequencing panel.

## Installation

The VH suite can be used by cloning or copying its repository to a given directory. Once downloaded, go to the VarHound directory and open the paths.ini file with any text editor. It looks kithe the following:

```python
[path]
# VarHound home absolute path
VarHound = ~/VarHound
# Coverage file suffix for high-level analysis
suffix = thresholds.bed
# Coverage file suffix for per-base analysis
pbcov = per-base.bed
# Relative path of TSO500 panel genes (BED file)
manifest = conf/TST500C_manifest.bed
# Relative path of reference genes (BED file)
reference = conf/UCSC_KnownGenes_hg19_TSO500_FPG.bed
```

The VH suite compute two kinds of coverage diagnostics: high-level and per-base. In the former, sequencing coverage is analyzed at exon level, while in the latter, the sequencing coverage profile is analyzed at single base level. Environmental variables are set to find the VH executable code, process thi input data at different levels, and use the required annotations (available in the VH repository).

Be sure that the **VarHound** variable is set to the absolute path of the VH directory inside your filesystem (by default it is set to the home directory: **~/VarHound**). The other environmental variables are defined as follows:

- **suffix**. The suffix used for the input coverage files. By default, filenames in the form xxx.thresholds.bed (or xxx.thresholds.bed.gz) are processed.
- **pbcov**. The suffix used for the input per-base coverage files. By default, filenames in the form xxx.per-base.bed (or xxx.per-base.bed.gz) are processed.
- **manifest**. The relative path of the BED format file containing the genomic coordinates (hg19 assembly) of all the regions covered by the TSO500 panel.
- **reference**. Reference gene coordinates (hg19 assembly). By default, the [**UCSC Genome Browser KnownGenes**](https://genome.ucsc.edu/) dataset is used.

## Coverage diagnostics for the TSO500 assay.

...

## VCF file indexing.

Under construction.

## Links to other repos

To be done.

## Disclaimer

VarHound is an ongoing project for TSO500 validation. Although I am happy to share the code with the community, it is still in its embryo stage!
