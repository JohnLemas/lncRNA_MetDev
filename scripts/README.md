# Template script directory
### John Lemas
### 5-8-2026
---
This directory outlines the scripts specific to lncRNA prediction using public resources. These scripts should be copied into `${BASE}/species/02_scripts/` for editing.
```bash
# use a wild card to correct file paths and file names quickly for each sbatch script. Example:
cp ${BASE}/lncRNA_met_dev/scripts/*sb ${BASE}/A.tha/02_scripts
sed -i 's/species/A.tha/g' ${BASE}/A.tha/02_scripts/*sb
# Make sure to update the base working directory for each script
```
---
### Contents
- [CD-HIT](#cd-hit)
  - removes duplicate sequences from truth dataset
- [NCBI datasets](#ncbi-datasets)
  - clean up `datasets` reference assembly and annotation files for use in the pipeline
- [Minimap2](#minimap2)
  - mapping predictions against truth sequences for pipeline analysis
- [Prefetch array](#prefetch-array)
  - uses sra accessions list to run `prefetch` and `fasterq-dump` as an array
---
### CD-HIT
INPUTS
- the combined list of all sequences in the 'truth' dataset

COMMAND
```bash
cd-hit -i ${INPUT} -o ${OUTPUT} -c 0.9 -M 800 -T 8 -l 200
```
FLAGS
- -c
  - sequence identity threshold, default 0.9
- -M
  - memory limit (in MB) for the program, default 800; 0 for unlimitted
- -T
  - number of threads, default 1; with 0, all CPUs will be used
- -l
  - length of throw_away_sequences, default 10

NOTES
- very fast
  - 8 cpu, 1 G and 2 hours
---
### Minimap2
INPUTS
- Predicted lncRNA sequences
- Unfiltered lncRNA candidate sequences
- condensed 'truth' dataset sequences
- Optional: Annotated lncRNA sequences from reference assembly
- Output file path

COMMAND
```bash
minimap2 -x cdna -t 8 $TRUTH $PRED -o $OUT_FILE
```
FLAGS
- -x
  - sets the input sequence type

SLURM NOTES
32 G, 8 threads, 30 minutes
---
### NCBI datasets
INPUTS
- unzipped `ncbi_dataset` directory
- species designation (ie z.mays)

COMMAND
```bash
# Unzip the dataset
unzip ncbi_dataset.zip

# Run the clean up script
sbatch datasets_cleanup.sb
```
This script uses AGAT and SeqKit to output required input files:
1. Nuclear scaffolds from the reference assembly
  - Do not want mitochondrial or chloroplast scaffolds in the pipeline
2. An mRNA-only annotation file in `gtf` format for `StringTie`
  - Only want to filter out transcripts that belong to annotated proteins
3. Any lncRNA annotations and sequences from the reference assembly
  - These are optionally included in `Minimap2` as a pseudo-positive control during benchmarking
4. A species karyotype file
  - This file can be used along with the final annotation file to map lncRNA loci onto the genome assembly using `RIdeogram`

This script also gives the option of including features annotated as `antisense_RNA` in the set of lncRNAs. Note that both `AGAT` and `SeqKit` have available containers if modules are not available. 

*`AGAT 1.4.0` and `SeqKit 2.3.1` were used for pipeline development*

```bash
# Output a clean gff3 file from the nuclear-only annotation
agat_convert_sp_gxf2gxf.pl -g "$GFF_NUC" -o "$CLEAN"
```
---
### Prefetch array
INPUTS
- List of SRA accessions
- output directory

COMMAND
``` bash
# 1. Run Prefetch (only if not already downloaded)
prefetch $ID

# 2. Run fasterq-dump
echo "Starting fasterq-dump for $ID..."
fasterq-dump $ID \
    --outdir $OUTDIR \
    --split-files \
    --skip-technical \
    --threads $SLURM_CPUS_PER_TASK \
    --temp /tmp

# 3. Gzip the results
if [ $? -eq 0 ]; then
    echo "Conversion complete. Compressing files..."
    gzip $OUTDIR/${ID}*.fastq

else
    echo "fasterq-dump failed for $ID"
    exit 1
fi

```
FLAGS
- -skip-technical
  - only keep biological reads
- -split-files
  - 

NOTES
- `pigz` is a faster alternative to gzip if it is available.
- Define the `/tmp` space path ahead of time for NCBI downloads. If running on HPCC the default will be $USER/$HOME which may have a small quota.
---