# VarHound (VH)

VH is a suite for genomic variant storage, retrieval, and validation through coverage analysis. It includes Python/R tools for the storage, management, and retrieval of large VCF files. In addition, VH offers a set of sample-level and region-level validation tools, through coverage analysis. VH is originally designed for truSight oncology 500 (TSO500) assays, but can be applied to any other DNA/RNA sequencing panel.

## 1. Installation

The VH suite can be used by cloning or copying its repository to a given directory. Once downloaded, go to the **VarHound** directory and open the **paths.ini** file with any text editor. It looks like the following:

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

The **VarHound** variable must be set to the absolute path of the VH directory inside your filesystem (by default, it is set to the home directory: **~/VarHound**). The VH suite compute two kinds of coverage diagnostics: **high-level** and **per-base**. In the former, sequencing coverage is analyzed at **exon level**, while in the latter, the sequencing coverage profile is analyzed at **single base level**. The environmental variables required for these analyses are defined as follows:


- **suffix**. The suffix used for the input coverage files. By default, filenames in the form xxx.thresholds.bed (or xxx.thresholds.bed.gz) are processed.
- **pbcov**. The suffix used for the input per-base coverage files. By default, filenames in the form xxx.per-base.bed (or xxx.per-base.bed.gz) are processed.
- **manifest**. The relative path of the BED format file containing the genomic coordinates (hg19 assembly) of all the regions covered by the TSO500 panel.
- **reference**. Reference gene coordinates (hg19 assembly). By default, a restricted panel of 64 genes from the [**UCSC Genome Browser**](https://genome.ucsc.edu/) KnownGenes database is used. The full KnownGenes dataset (conf/UCSC_KnownGenes_hg19_canonical.bed) is also available.

Once the VH variables are set, the user should add the VH absolute path to the PATH environmental variable:
To be done.

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

High-level coverage diagnostics can be enabled using the `vhLauch.py` command, followed by the input directory, and the variant type (`snv`, `cnv`, `rna`, `fusion`):

```
vhLaunch.py ~/TSO500 snv
```

Input coverage files will be recursively searched within the given input direcory and subdirectories, recruiting all those files with the suffix specified by the environmental variables (by default, "thresholds.bed" and "thresholds.bed.gz"). The input folder is considered as a single *run*, where each coverage file comes from a subject. The thresholds.bed file should have the following mandatory fields: chromosome, region (exon) start, region (exon) end, region ID, 5x depth, 10x depth, 50x depth, 100x depth, 250x depth, and 500x depth. Coverages at a given depth can be obtained using SAM/BAM analysis tools, such as [**samtools**](https://github.com/samtools), [**deepTools**](https://deeptools.readthedocs.io/en/develop), and [**mosdepth**](https://github.com/brentp/mosdepth).

For each region, for a given depth, VarHound computes the number of bases covered at that depth (count data). Coverage count data is then converted into frequency data (i.e., the percent coverage, *p*), dividing by region length. VarHound coverage analysis starts from a *covdata* file, including the input fields and median values of *p* across samples (i.e., subjects), referred to as **median percent coverage** (MPC). As a good quality principle, a run should maximize the number of depths with MPC ≥ 75%. Similarly, at a given depth *D*, for each sample *x*, the first quartile of the exon percent coverage distribution should be ≥ 75% (i.e., at least 75% of the exons of *x* should be covered by at least 75%, at depth *D*).

The outputs include:
- run-level barplots showing the MPC distribution for each depth,
- subject-level boxplot reporting the distribution of exon MPC (y axis) for each subject (x axis),
- boxplot of exon MPC (y axis) for each gene (x axis),
- a blacklisted region (exon) BED file for each depth, including the following fields: chromosome, region start, region end, region name, gene symbol, genomic element, region id (RefSeq ID, in case of genes), entry ID (internal usage), time (internal usage), Q1 (first MPC quartile; it should be < 75), depth of coverage.

### 2.2. Base-level coverage diagnostics

High-level coverage analysis cannot find small coverage drops (from 1 bp to roughly 100 bp) within highly covered regions. However, narrow coverage drops might hide false negative variant calls in apparently covered regions, potentially leading to undiagnosed pathogenic or malignant traits.

High-level coverage diagnostics can be enabled using the `vhLauch.py` command, followed by the input directory, and the desired coverage threshold preceeded by an "x" (in the example below, all the bases below 80 reads of coverage will be marked as drops):

```
vhLaunch.py ~/TSO500 x80
```

Here, the algorithm will recursively search for base-coverage files, as indicated by the *pbcov* environmental variable (by default, these file should have a "per-base.bed" suffix). These files must have the following columns: chromosome, start, end, coverage value.

This analysis generates three output files:
- A coverage drops BED file (fields: chromosome, drop start, drop end, region name, median coverage of the drop)
- A coverage drops UCSC Genome Browser BED track (fields: chromosome, start, end, name, score as the median coverage of the region, strand, thickStart, thickEnd, itemRgb). This file can be directly uploaded and viewed using the UCSC Genome Browser track viewer.
- A gene drops BED file (fields: chromosome, gene start, gene end, gene symbol, gene exon count, strand, number of exon drops at the given coverage threshold).

## 3. VCF file descriptive plots and tables.

The VH suite includes the `vhReadAnnotations.py` script to create a number of descriptive plots and tables for annotated [**VCF files**](https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/).
For a basic usage, the script only requires a folder with annotated VCF files:

```
vhReadAnnotations.py /path_to_annotated_VCF_files
```

The chosen path must contain VCF files without headers (i.e., plain tab-separated text). Headers can be recursively removed with the `cleanAnnotations.py` script. Annotated VCF files must have the following [**GDC-compliant**](https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format) mandatory fields:

- `--identifier`  A vector of length 2 specifying the gene name field and the generic label name used for genewise data aggregation. Deafult: ['SYMBOL', 'labels']
- `--vclass`  Variant class. Default: 'VARIANT_CLASS'
- `--impact`  Variant impact. Default: 'IMPACT'
- `--order`  Impact in degreasing order. Default: ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

Two additional annotation fields are requred: `CLIN_SIG` and `ClinVar_CLNSIG`. An annotated VCF containing these fields can be obtained using the [**variant effect predictor (VEP)**](https://grch37.ensembl.org/info/docs/tools/vep) tool. A full pipeline for genomic variant calling and annotation is available at the [**nfcore Sarek website**](https://nf-co.re/sarek). Additional options can be viewed by launching `vhReadAnnotations.py -h`

Descriptive plots and tables are generated both for each VCF file and for the entire cohort (represented by the files within the input directory). The whole output is organized within the input directory, according to the following scheme:

```
            +-- Input_annotated_VCF1                   
            +-- Input_annotated_VCF2                   
            +-- ...
            |
Input_dir --+-- Annotated_VCF1_outdir
            |     |
            |     +-- bySample
            |     |    +-- overall_variants_barplot
            |     |    +-- overall_variants_stackedbars
            |     |    +-- topk_high_impact_genes_barplot
            |     |    +-- topk_mutated_genes_barplot
            |     |    +-- ClinSig_impact_piechart
            |     |    +-- ClinVar_impact_piechart
            |     |    +-- Counts_directory
            |     |         |
            |     |         +-- Variant_impact_table_by_gene
            |     |         +-- Variant_impact_table_by_class
            |     |         +-- ClinSig_impact_freqs
            |     |         +-- ClinVar_impact_freqs
            |     |
            |     +-- byClass
            |          |
            |          +-- SNV
            |          |    +-- Variants_impact_barplot
            |          |    +-- topk_high_impact_genes_barplot
            |          |    +-- topk_mutated_genes_barplot
            |          |    +-- Counts_directory
            |          |         |
            |          |         +-- Variant_impact_table_by_gene
            |          |         +-- Variant_impact_freqs
            |          |
            |          +-- deletion
            |          |    +-- ...
            |          |
            |          +-- insertion
            |          |    +-- ...
            |          |
            |          +-- ...
            |
            +-- Annotated_VCF2_outdir
            |     +-- ...
            |
            +-- Annotated_VCFn_outdir
            |     +-- ...
            |
            +-- Overall
            |    |
            |    +-- SNV
            |    |    +-- Impact_piechart
            |    |    +-- topk_high_impact_genes_barplot
            |    |    +-- topk_mutated_genes_barplot
            |    |    +-- Genes_annotation_table
            |    |
            |    +-- deletion
            |    |    +-- ...
            |    |
            |    +-- insertion
            |    |    +-- ...
            |    |
            |    +-- ...
            |
            +-- trash (intermediate files)
```

## 4. VCF file indexing.

Under construction.

## 5. Links to other repos

**Data indexing**. [**Kew**](https://github.com/fernandoPalluzzi/KewTools) is a simple command line tool for [**SQLite**](https://sqlite.org) database creation and querying. It enables the creation and manipulation of big files and directories, also for non SQL or Python experts, by indexing them with the [**sqlite3**](https://docs.python.org/3.8/library/sqlite3.html) Python library.
