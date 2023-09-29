# Crykey: cryptic lineage detection for SARS-CoV-2 wastewater surveillance 

## Summary

Wastewater monitoring is an important tool that can complement clinical testing for COVID-19 surveillance. This repository contains an implementation of a novel bioinformatics tool (Crykey), for identifying SARS-CoV-2 mutations that are rare or entirely missing in existing databases (crytpic mutations) from wastewater. The benefits of QuaID are three-fold: 

1. identifies co-occurring mutations in wastewater samples that potentially being cryptic
2. enables fast queries of mutation combinations against all publicly available SARS-CoV-2 genomes
3. improving understanding of SARS-CoV-2 intrahost evolution and transmission events at a large scale

It is highly recommanded that the wastewater samples are processed with [QuaID](https://gitlab.com/treangenlab/quaid), a novel bioinformatics tool (QuaID) for VoC detection based on quasiunique mutations that are being developed by [Treangenlab](https://gitlab.com/treangenlab).

## 3rd party software requirements 

Below is the list of 3rd party software requirements. All required sofwtare can be installed via Miniconda after adding `bioconda` to the list of channels. Version specified in the parantheses is the version currently tested.

* vdb (2.7)
* samtools (1.7)
* SnpEff

## Python requirements

* numpy (1.19.0)
* pandas (1.0.5)
* pyvcf 
* pysam
* tqdm
* dask
* biopython

## Operating system requirements

The code has been tested on MacOS Monterey 12.4 (Intel MacBook) and Ubuntu 18.04.6 (kernel version 4.15.0-171-generic).

## Usage

### Download pre-build Crykey database

A pre-build database based on publicly available SARS-CoV-2 genomes till Jan, 10, 2023 is provided and can be [downloaded from this link](https://rice.box.com/s/1rdad8b6q0f4uccmqohqczcvkcf65mgb). If you are using the pre-build Crykey database, you can now skip to [Searching co-occurring mutations in wastewater samples](https://github.com/treangenlab/crykey/tree/main#searching-co-occurring-mutations-in-wastewater-samples).

### MSA preprocessing and vdb precomputation work

The first step of building a Crykey database is the download MSA and metadata provided in [GISAID EpiCov](https://gisaid.org/). After the MSA and metadata are downloaded. Some pre-computation work has to be done. If you are running [QuaID](https://gitlab.com/treangenlab/quaid), this step should be already included while building the QuaID database. If you wish to use your own BAM/VCF files generated with custom-made bioinformatics pipelines directly, please refer to [msa-preprocessing-and-database-construction](https://gitlab.com/treangenlab/quaid#msa-preprocessing-and-database-construction) for details.

### Crykey database construction

After the pre-processing step, you are now ready to construct the Crykey database. To do so, run the following python script:

```
usage: crykey_build_db.py [-h] -m METADATA -l LINEAGES -v VDB -d DATE
                          [-o OUTPUT] [--exclude-recombinant]

Building Crykey Database based on GISAID MSA.

optional arguments:
  -h, --help            show this help message and exit
  -m METADATA, --metadata METADATA
                        Path to the metadata of the GISAID MSA in tsv format.
  -l LINEAGES, --lineages LINEAGES
                        Path to the vdb lineage dataframe.
  -v VDB, --vdb VDB     Path to the trimmed vdb file.
  -d DATE, --date DATE  Date of the lastest record in the GISAID MSA, with the
                        format [YYYY-MM-DD].
  -o OUTPUT, --output OUTPUT
                        Output path of the Crykey Database.
                        Default:[crykey_dbs]
  --exclude-recombinant
                        Exclude recombinant lineages. Default:[False]
```

The required files are:
* metadata.tsv (downloaded from GISAID EpiCoV)
* vdb_lineage_df_week.csv (generated in the previous step)
* vdb_trimmed_nucl.txt (generated in the previous step)

You would also need to input the date (usually also serves as version number for MSA you have downloaded). Any records after the date will be ignored. 

The `crykey_build_db.py` script will generate three files:
* mutation0_lookup.pkl
* mutation50_lookup.pkl
* quarc_db.pkl

### Searching co-occurring mutations in wastewater samples

Once the database construction is complete, you can process your wastewater samples with the following script:

```
usage: crykey_wastewater.py [-h] -i METADATA -r REFERENCE [-d DATABASE]
                            [-o OUTPUT]

Search for cryptic lineages in wastewater samples.

optional arguments:
  -h, --help            show this help message and exit
  -i METADATA, --metadata METADATA
                        Path to the input metadata, tab seperated dataframe
                        with the following field on each row:
                        Sample_Collection_Date[MMDDYYYY], WWTP[str],
                        Sorted_BAM[path-to-sorted-bam-file], VCF[path-to-VCF].
  -r REFERENCE, --reference REFERENCE
                        Path to the Reference Genome.
  -d DATABASE, --database DATABASE
                        Path to the Crykey Database Directory.
                        Default:[crykey_dbs]
  -o OUTPUT, --output OUTPUT
                        Output directory. Default:[crykey_output]
```

An input metadata file is required for the software to associated BAM files and VCF files in the correct path. The input metadata file is formatted as following as a tab seperated dataframe with corresponding headers. If you are using QuaID, you should have all those files ready to use.

| Sample_Collection_Date | WWTP | Sorted_BAM | VCF |
|------------------------|------|------------|-----|
| MMDDYYYY | user-defined WWTP identifier | path to the sorted BAM file | path to the VCF file |

The output of this step should have the following file structure:

```
crykey_output/
|
 -- cryptic_alignment/ 
    |
     -- cryptic_reads_[date]_[WWTP].bam
     -- cryptic_reads_[date]_[WWTP].bam.bai
     -- ...
 -- cryptic_dataframe/ 
    |
     -- [date]_[WWTP].csv
     -- ...
 -- cryptic_reads/ 
    |
     -- read_ids_[date]_[WWTP].txt
     -- ...
 -- cryptic_vcfs/ 
    |
     -- [date]_[WWTP]_mutations.ann.vcf
     -- [date]_[WWTP]_mutations.vcf
     -- ...
```

The output contains information of the alignment, allele frequency, read ids of co-occuring mutations in your data. They can be useful if you are planning to do downstream analysis on certain cryptic lineages.

### Querying Crykey Database to determine the rarety of the co-occuring mutations

In the final step, run the following script to query all co-occuring mutation sets that are found in your data against Crykey database, to calculate the occurrence of each co-occuring mutation sets in the known SARS-CoV-2 genomes. 

```
usage: crykey_query.py [-h] [-d DATABASE] [-o OUTPUT]

Search for cryptic lineages in wastewater samples.

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Path to the Crykey Database Directory.
                        Default:[crykey_dbs]
  -o OUTPUT, --output OUTPUT
                        Output directory. Default:[crykey_output]
```

The script merges all information from all of the wastewater samples, and queries each unique co-occuring mutation sets against the database. The step will creates `final_result.csv` under the output directory. `final_result.csv` file will be in comma-separated values (CSV) format and will contain the following columns:

| Date | Site | Nt Mutations | AA Mutations | Support DP | Total DP | Combined Freq | Overall Occurrence | Lineage Occurrence |
|------|------|--------------|--------------|------------|----------|---------------|--------------------|--------------------|
| Sample collection date | user-defined WWTP identifier | multiple SNVs represented as [REF1][location1][ALT1];[REF2][location2][ALT2];... | the amino acid changes of the corresponding SNVs | number of reads that support all SNVs at the same time in the sample | total number of reads that span all SNV locations in the sample | allele frequency of the co-occurring SNVs | the overall occurrence of the co-occurring SNVs in the database | detailed information of the lineages that contain the co-occurring SNVs, in the format of [lineage1]:[occurrence1];[lineage2]:[occurrence2];... |

Based on such information, you could determine which of the co-occurring SNVs qualify as cryptic lineage. We recommand that the qualified cryptic lineage should:

* occurred in multiple WWTP samples,
* have sufficiant number of supporting reads,
* the occurence in the database should be low. in other words, the cryptic lineage should be novo or at least rare in the database. 

## Manuscript

You can find the manuscript describing QuaID and corresponding results at [doi.org/10.1101/2023.06.16.23291524](https://www.medrxiv.org/content/10.1101/2023.06.16.23291524v1).