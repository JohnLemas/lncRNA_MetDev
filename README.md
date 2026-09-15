# IWGC lncRNA Prediction
## Expanding annotations for non-model plants
### John Lemas
### 5-6-2025
---
## Citation
 This pipeline was originally developed by Xue-Chan Tian *et al*: \
 Original publication: [DOI](https://doi.org/10.1093/hr/uhae041) \
 Author GitHub: [Plant-LncPipe-v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2/tree/master)
 
---
### Table of contents
- [General overview](#general-overview)
- [Working Directory](#working-directory)
- [Gathering inputs](#gathering-inputs)
- [Set up](#set-up)
- [Running the pipeline](#running-the-pipeline)

---
## General overview
1. Obtain 'truth' dataset from [PLncDB](https://www.tobaccodb.org/plncdb/) for a given speices
2. Generate a `PMID` list associated with RNA-seq projects under the 'species' tab in PLncDB
3. Use `EDirect` and `SRA-Toolkit` to get raw data from each of these projects
4. Download a reference genome and annotation file
5. Trim raw `fastq` files and align trimmed reads to the reference genome
6. Assemble a reference-guided *de novo* transcriptome using `StringTie`
7. Obtain candidate lncRNAs using `grep` and `gffread`
8. Run these candidates through `FEELnc`, `lncBoost`, `lncFinder`, `cpat`, and `DIAMOND`
9. Intersect IDs that passed every tool to obtain high-confidence lncRNAs
10. Use sequence alignment to compare high-confidence *ab initio* candidates to 'truth' dataset
 ---
## Working Directory
- def_files
  - definition files to create necessary images for each tool
- scripts
  - sbatch script templates for hpcc
  - custom `R` and `python` scripts
- species
  - 01_inputs
    - PMID metadata
    - list of sra IDs saved as `species_sra_IDs.txt`
    - reference genome (`species.chromosomes.fa`) and annotation (`species_mRNA.gtf`)
  - 02_scripts
    - scripts copied from parent script directory (see `scripts/README.md`)
  - 03_outputs
    - pipeline output directories
- truth_data
  - lncRNA datasets gathered from [PLncDB](https://www.tobaccodb.org/plncdb/)
  - species directories
    - individual database lncRNA `fasta` and `gff3` files
    - full concatonated lncRNA files
    - condensed lncRNA `fasta` files (see `scripts/README.md` `cd-hit` step)
 ---
## Gathering inputs
- Navigate to [PLncDB](https://www.tobaccodb.org/plncdb/)
- Gather PMID information from the RNA-Seq table for a given species
  - Use EDirect to gather SRA accession numbers from PMIDs; save as list `species_sra_IDs.txt`
- Gather truth sequences from the downloads>species page on PLncDB
  - Download the lncRNA fasta files and gff3 files
- Create working directory
  - Clone this repository. Clone [Plant-lncPipe](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2) and [FEELnc](https://github.com/tderrien/feelnc) into the `scripts` directory.
  - Run `scripts/build_images.sb` array job to build each necessary image
  - Create the working directory for the species
    - `mkdir 01_inputs 02_scripts 03_outputs`
  - Copy scripts into `species/02_scripts` and run the `sed`  command (see `scripts/README.md`)
---
## Set up
### Prepare the reference genome assembly
The NCBI tool `datasets-cli` is an easy way to download a reference genome via command line; however the assembly usually contains messy fasta headers and gff files. The following steps can be used to obtain the latest reference assembly for a species and clean it for use in this pipeline.

```bash
# Create an environment with datasets-cli
ml purge
ml load Miniforge3
conda create -n datasets -c conda-forge ncbi-datasets-cli
conda activate datasets

# General datasets command
datasets download genome taxon "[Scientific Name]" --reference --include genome,gff3
unzip ncbi_dataset.zip # This will result in a directory, an md5sum file, and a README.md file
```

The script `ncbi_datasets_cleanup.py` pulls full chromosomes from the genome, filters the gff file, renames chromosome entries for both files so that they match, then outputs the two required downstream files: `species.chromosomes.fa` & `species_mRNA.gtf`.

```bash
# Clean up the raw download
cd ncbi_dataset/ # move into the unzipped directory
python3 /path/to/ncbi_datasets_cleanup.py [species] # run the cleanup script with the chosen species designation as the first argument
```

Move the new `chromosomes.fa` and `_mRNA.gtf` files into the `species/01_inputs` directory for indexing, and remove the ncbi_dataset directory, README file, and md5sum file from the datasets download. Optional: keep the compressed datasets download for later reference.

### Create the UniProt database for DIAMOND:
```bash
# Download the latest UniProt database release
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# Expand the file
gunzip uniprot_sprot.fasta.gz

# Use DIAMOND to make the reference database
ml purge
ml load DIAMOND
diamond makedb --in uniprot_sprot.fasta -d uniprot_out
```
---
## Preparing the truth datasets
Truth datasets are individually packaged from each respective source (ie NCBI, PLncBD, EVlncRNA, RNAcentral, etcetera). The goal of this step is to concatonate all of these input datasets into one master set, and condense this set using `CD-HIT` to filter out any lncRNAs that are present in more than one dataset. 

```bash
# Combine all fasta and gff files into one master set
cat IWGC_lncRNA_pipeline/truth_data/species/*.fa >> IWGC_lncRNA_pipeline/truth_data/species/species_PLncDB_lncRNA.fa

# Run the script to condense the master set
sbatch IWGC_lncRNA_pipeline/species/02_scripts/cd-hit.sb

```

This master set should be exhaustive of the known experimentally validated and high confidence predicted lncRNAs already published.

---
## Running the pipeline
**Detailed notes on running each step can be found in `scripts/README.md`**

Once all the inputs have been gathered and the species working directory is created, the order of operations is as follows:

### Steps 1 - 7: Preparing model inputs
1. `star_genomegenerate.sb`
2. `prefetch_array.sb`
3. `fastp.sb`
4. `alignment_array.sb`
5. `bam_indexing.sb`
6. `stringtie_assemble.sb`

### Step 7: Merging assemblies and gathering candidates

7. `stringtie_merge.sb`

Once the merge step is complete, the output `gtf` file contains all of the aligned sequence locations and any information about an overlapping mRNA in that region. This script parses out transcripts that did not overlap with mRNA locations in the master output gtf file and outputs `candidates.gtf` and `candidates.fasta`. These unfiltered candidate locations, sequences, and stringtie transcript IDs are used in all of the following *ab initio* prediction models.

### Steps 8 - 12: Running the models
*Note that each of the models can run simultaneously*

8. `FEELnc_filter.sb`
9. `cpat.sb`
10. `lncFinder.sb`
11. `lncBoost.sb`
12. `DIAMOND.sb`

### Step 13: Intersect predictions and classify high-confidence lncRNAs

13. `intersect_classify.sb`

This step merges output from all ML tools producing a `.txt` file with `StringTie` gene IDs for candidates that passed all filters and a table of concensus information. It then uses `FEELnc_classify.pl` to classify each predicted lncRNA. `intersect_classify.sb` eliminates the need for multiple R based images and handles all of the data parsing and file generation in one ~5 minute step using `intersect_ids.R` and `classifications.R`. 

### Step 14: Compare pipeline results to truth dataset

14. `minimap2.sb`

---
