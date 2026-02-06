# Habit-free-OC_Analysis_APJCP
Codes used for the manuscript titled "Genomic landscape of oral squamous cell carcinoma in never smokers and never drinkers" for Asian Pacific Journal of Cancer Prevention, 2026, is deposited here
##TCGA processed dataset and NCBI RAW datset was used for the analysis.
## Objectives
- Perform **molecular subtyping** of oral cancer samples using transcriptomic data.
- Conduct **RNA-seq analysis** for differential expression and pathway enrichment.
- somatic variant calling pipelines** (Mutect2) for mutation profiling.

## WES data Analysis steps for the replication study
- Quality control (FastQC)
- Alignment (BWA-MEM, GRCh38)
- Somatic variant calling (GATK Mutect2, Tumor only with PoN)
- Variant annotation (ANNOVAR )
- Downstream analyses (TMB, COSMIC signatures, Survival analysis)
### mutect2_somatic_pipeline/
├── README.md
├── config/
│   └── params.sh
├── scripts/
│   └── mutect2_somatic.sh
├── resources/
│   ├── pon.vcf.gz
│   └── germline_resource.vcf.gz
└──bash/
    └── run.sh

## RNA-seq Analysis Pipeline 

- Normalization and differential expression (DESeq2)
- Functional enrichment (GO, GSEA, pathway analysis)
- Visualization (PCA, heatmaps, volcano plots)



## Repository Structure
- **Molecular_subtype.R**  
  Script for molecular subtyping based on gene expression profiles.
  
- **Rnaseq.R**  
  Workflow for RNA sequencing analysis, including preprocessing, normalization, and differential expression.
  
- **Somatic Variant Calling Pipeline (Mutect2)**  
  End-to-end pipeline for somatic mutation detection using GATK Mutect2.
  
- **README.md**  
  Project overview and usage instructions.
  
- **SECURITY.md**  
  Guidelines for secure data handling and compliance.




## Getting Started
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/Habit-free-OC_Analysis_APJCP.git

Citation
If you use this pipeline, please cite: Sinumol George et al., [Asian Pacific Journal of Cancer Prevention, 2026]
