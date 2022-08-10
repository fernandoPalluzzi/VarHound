# VarHound (VH)

VH is a suite for genomic variant storage, retrieval, and validation through coverage analysis. It includes Python/R tools
for the storage, management, and retrieval of large VCF files. In addition, VH offers a set of sample-level and region-level
validation tools, through coverage analysis. VH is originally designed for truSight oncology 500 (TSO500) assays, but 
can be applied to any other DNA/RNA sequencing panel.

## 1. Installation

The VH suite can be used by cloning or copying its repository to a given directory. Once downloaded, go to the **VarHound** directory and open the **paths.ini** file with any text editor. It looks kithe the following:

```bash
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

The VH suite compute two kinds of coverage diagnostics: **high-level** and **per-base**. In the former, sequencing coverage is analyzed at **exon level**, while in the latter, the sequencing coverage profile is analyzed at **single base level**.

The **VarHound** variable must be set to the absolute path of the VH directory inside your filesystem (by default, it is set to the home directory: **~/VarHound**). The other environmental variables are defined as follows:

- **suffix**. The suffix used for the input coverage files. By default, filenames in the form xxx.thresholds.bed (or xxx.thresholds.bed.gz) are processed.
- **pbcov**. The suffix used for the input per-base coverage files. By default, filenames in the form xxx.per-base.bed (or xxx.per-base.bed.gz) are processed.
- **manifest**. The relative path of the BED format file containing the genomic coordinates (hg19 assembly) of all the regions covered by the TSO500 panel.
- **reference**. Reference gene coordinates (hg19 assembly). By default, a restricted panel of 64 genes from the [**UCSC Genome Browser**](https://genome.ucsc.edu/) KnownGenes database is used. The full KnownGenes dataset (conf/UCSC_KnownGenes_hg19_canonical.bed) is also available.

Once the VH variables are set, the user should add the VH absolute path to the PATH environmental variable:

```bash
# Extended PATH variables set
export PATH=$PATH:~/chosen_directory/VarHound
```

To work properly, the VH directory needs execution rights. If not (or not sure), they can be changed with:

```bash
chmod -R 700 ~/chosen_directory/VarHound
```

Now, from a Bash terminal, the `vhLaunch.py` command should give the following message:

```
          ## VarHound software ##
# Coverage Analysis (VarHound-CA) - v0.1.0 #

BASIC USAGE:
  vhLaunch.py COVERAGE_DIRECTORY [snv|cnv|rna|fusion]
Example:
  vhLaunch.py ~/TSO500 cnv

PER-BASE ANALYSIS:
  vhLaunch.py COVERAGE_DIRECTORY x<COVERAGE_THRESHOLD>
Example:
  vhLaunch.py ~/TSO500 x80
```

## 2. Coverage diagnostics for the TSO500 assay.

The VH coverage analysis tool takes in input a coverage directory and recursively search for (possibly gzipped) coverage input files, generating a number of useful reports. The input data can be placed anywhere in the directory tree, and can be easily generated with open-source tools, such as [**mosdepth**](https://github.com/brentp/mosdepth). The analysis scheme is shown below:

```
# INPUT SAMPLES #                      # ANALYSIS OUTPUT #
 
SNV --+                                +--- High-level analysis
      |                                |     +-- Whole-panel coverage
      +-- DNA --+                      |     +-- Sample-level coverage
      |         |                      |     +-- Gene-level coverage
      |         |                      |           +-- blacklisted regions (BED)
CNV --+         +--> RUN DIRECTORY --> + 
                |                      | 
                |                      +--- Base-level analysis
                |                            +-- Coverage drops BED file
                |                            +-- Gene drops BED file
          RNA --+                            +-- Coverage drops UCSC Genome Browser track file
```

### 2.1. High-level coverage diagnostics

High-level coverage diagnostics can be enabled using the `vhLauch.py` command specifying one of the available variant types (`snv`, `cnv`, `rna`, `fusion`):

```
vhLaunch.py ~/TSO500 snv
```

Input coverage files will be recursively searched within the given input direcory and subdirectories, recruiting all those files with the suffix specified by the environmental variables (by default, "thresholds.bed" and "thresholds.bed.gz"). The input folder is considered as a single *run*, where each coverage file comes from a subject. The thresholds.bed file should have the following mandatory fields: chromosome, region (exon) start, region (exon) end, region ID, 5x depth, 10x depth, 50x depth, 100x depth, 250x depth, and 500x depth. Coverages at a given depth can be obtained using SAM/BAM analysis tools, such as [**samtools**](https://github.com/samtools), [**deepTools**](https://deeptools.readthedocs.io/en/develop), and [**mosdepth**](https://github.com/brentp/mosdepth).

For each region, for a given depth, VarHound computes the number of bases covered at that depth (count data). Coverage count data is then converted into frequency data (i.e., the percent coverage, *p*), dividing by region length. VarHound coverage analysis starts from a *covdata* file, including the input fields and median values of *p* across samples (i.e., subjects), referred to as **median percent coverage** (MPC). As a good quality principle, a run should maximize the number of depths with MPC ≥ 75%. Similarly, at a given depth *D*, for each sample *x*, the first quartile of the exon percent coverage distribution should be ≥ 75% (i.e., at least 75% of the exons of *x* should be covered by at least 75%, at depth *D*).

The outputs include:
- run-level parplots showing the MPC distribution for each depth,
- subject-level boxplot reporting the distribution of exon MPC (y axis) for each subject (x axis),
- boxplot of exon MPC (y axis) for each gene (x axis),
- a blacklisted region (exon) BED file for each depth, including the following fields: chromosome, region start, region end, region name, gene symbol, genomic element, region id (RefSeq ID, in case of genes), entry ID (internal usage), time (internal usage), Q1 (first MPC quartile; it should be < 75), depth of coverage.

### 2.2. Base-level coverage diagnostics

...

## 3. VCF file indexing.

Under construction.

## 4. Links to other repos

**Data indexing**. [**Kew**](https://github.com/fernandoPalluzzi/KewTools) is a simple command line tool for [**SQLite**](https://sqlite.org) database creation and querying. It enables the creation and manipulation of big files and directories, also for non SQL or Python experts, by indexing them with the [**sqlite3**](https://docs.python.org/3.8/library/sqlite3.html) Python library.
