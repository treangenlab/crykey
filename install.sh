#!/bin/bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -c conda-forge biopython">=1.79" dask tqdm pandas"=1.0.5" numpy"=1.19.0"
conda install -c bioconda pysam pyvcf">=0.6.8" samtools">=1.13" snpeff
