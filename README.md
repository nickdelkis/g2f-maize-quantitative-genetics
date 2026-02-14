# Quantitative Genetics Analysis of Maize Using G2F Data

## Overview
This repository demonstrates complementary approaches to quantitative genetics using 
the Genomes to Fields maize dataset: genetic discovery through GWAS and phenotype 
prediction through machine learning.

## Dataset
Genomes to Fields (G2F) Initiative: 5,899 hybrid genotypes, 2,425 SNPs, 
phenotypic data from 272 environments (2014-2023).

## Projects

### 1. GWAS: Genetic Architecture of Flowering Time
- **Goal**: Identify loci associated with flowering time variation
- **Methods**: GWAS using GAPIT, 2,425 SNPs, 2022-2023 environments
- **Results**: [Summary findings]
- **[View detailed analysis →](01-gwas-flowering-time/)**

### 2. Machine Learning: Yield Prediction
- **Goal**: Predict grain yield from genomic and environmental data
- **Methods**: Random Forest and Neural Network, 673 environmental covariates
- **Results**: [Model performance]
- **[View detailed analysis →](02-ml-yield-prediction/)**

## Why Both Approaches?
In modern breeding programs, GWAS identifies causal biology while ML models 
optimize selection. This portfolio demonstrates capability in both discovery 
and application.

## Setup Instructions

### Prerequisites
- R (>= 4.3.0)
- Python (>= 3.12.0)
- Git

### Installation

#### 1. Clone Repository
```bash
git clone https://github.com/YOUR_USERNAME/g2f-maize-quantitative-genetics.git
cd g2f-maize-quantitative-genetics
```

#### 2. Download Data
See [data/README.md](data/README.md) for data download instructions.

#### 3. R Environment Setup
Open RStudio in this directory:
```r
# Install renv if needed
install.packages("renv")

# Restore R packages from lockfile
renv::restore()
```

This installs all R packages listed in `renv.lock`.

#### 4. Python Environment Setup
```bash
# Create virtual environment
python -m venv venv

# Activate environment
source venv/Scripts/activate  # Git Bash on Windows
# OR
venv\Scripts\activate  # PowerShell on Windows

# Install packages
pip install -r requirements.txt
```