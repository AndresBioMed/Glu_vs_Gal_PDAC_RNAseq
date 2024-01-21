# Transcriptomic Analysis: Adherent Media vs Sphere Media in PDAC

## Overview

This repository contains the code and documentation for a comprehensive transcriptomic analysis comparing the effects of Adherent Media vs Sphere Media in Pancreatic Ductal Adenocarcinoma (PDAC). The analysis includes a fully crafted Transcriptomics Protocol, from read mapping and filtering to differential gene expression, ontology enrichment, and qualitative isoform usage analysis.

## Author

- **Andrés Gordo Ortiz**

## Date

- Report generated on: `current date`

## Repository Structure

The repository is organized into folders reflecting different stages and aspects of the analysis.

### 1. Mapping Folder

- **Kallisto**: Output from the pseudoaligner Kallisto for reads to the reference genome.
- **FastQC**: Reports on the quality and depth of reads for each sample.
- **MultiQC**: Summarized report of the quality, depth, and alignment success for all samples.

### 2. Report Folder

#### Data

- **log2_cpm.csv**: Counts per million of all genes across samples, non-filtered, non-normalized.
- **log2_cpm_filtered.csv**: Counts per million of all genes across samples, filtered, non-normalized.
- **log2_cpm_filtered_norm.csv**: Counts per million of all genes across samples, filtered, normalized.
- **pca_results.csv**: Results and loadings of the PCA analysis.

#### IsoformSwitchAnalyzer

- **common_switch_consequences.pdf**: Summary of all isoform switching events.
- Folder with plots from Isoform Switch-positive genes.

#### Plots

- **Volcano Plot**: Expression and statistical significance of every gene.
- **Bubble Plot C2 and H**: Top 20 enriched terms on the C2 collection or Hallmark collection.
- **Manhattan Plot**: Enriched terms on Gene Ontology.
- **GSEA**: Gene Set Enrichment Plots from the top enriched terms on the C2 and H collections.

#### Heatmaps

- **Heatmap_Full**: Gene expression across all samples and modules.
- **Heatmap_Up**: Gene expression across all samples for upregulated genes.
- **Heatmap_Down**: Gene expression across all samples for downregulated genes.

#### Multivariate Analysis

- **PCA**: Principal Component Analysis plot.
- **PCA SM**: Small Multiples Plot of loadings for each principal component.
- **All PCA**: Combined plot of the above two.
- **Scatter**: Scatter plot showing the average expression of each gene between groups.

#### Quality of Reads

- Various violin plots representing Counts per million at different preprocessing stages.

#### Tables

- **differentialexpression_results.csv**: Results of the differential expression analysis.
- Other tables related to gene expression, isoform switching, and functional enrichment.

## Alignment Protocol

### 1. Data Acquisition

- Raw sequencing reads obtained from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB10204).

### 2. Read Mapping

- Reads mapped to the human reference transcriptome [GRCh38.cdna](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz) using [Kallisto](https://pachterlab.github.io/kallisto/) version 0.48.

### 3. Quality Analysis

- Quality assessed using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiqc](https://multiqc.info/).

### 4. Experimental Design

- 20 single-end samples for two experimental conditions: with adherent media and with sphere media.

### 5. Pseudo-alignment and Automated Script

- Pseudo-alignment performed using Kallisto. Custom automated script available [here](https://github.com/AndresBioMed/Glu_vs_Gal_PDAC_RNAseq/blob/main/AutomatedKallistoGeneAlignment.sh).

### 6. Pre-processing

- Quality evaluation with *fastqc* for each sample.

### 7. Data Integration

- [TxImport](https://bioconductor.org/packages/release/bioc/html/tximport.html) used to import Kallisto outputs into the R environment.

### 8. Data Summarization

- Annotation data from Biomart used to summarize data from transcript-level to gene-level.

## Usage

- Clone the repository and follow the provided analysis protocol for a detailed transcriptomic investigation.

## License

This project is licensed under the [MIT License](LICENSE).

Feel free to reach out for any questions or collaboration opportunities.

Andrés 
