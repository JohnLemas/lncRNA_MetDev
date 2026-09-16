# Template script directory
### John Lemas
### 5-8-2026
---
This directory contains all scripts used to predict and verify lncRNAs in a given species. Scripts can be copied into a species' sub-directory for editing.
```bash
# use a wild card to correct file paths and file names quickly for each sbatch script. Example:
cp IWGC_lncRNA_pipeline/scripts/* IWGC_lncRNA_pipeline/A.thaliana/02_scripts
sed -i 's/species/A.thaliana/g' IWGC_lncRNA_pipeline/A.thaliana/02_scripts/*
```
---
### Contents
- [Alignment array](#alignment-array)
  - Array `SLURM` job that uses an sra list to align trimmed `fastq` files
- [Bam indexing](#bam-indexing)
  - Array `SLURM` job to index `bam` files
- [CD-HIT](#cd-hit)
  - removes duplicate sequences from truth dataset
- [CPAT plant](#cpat-plant)
  - Logistic regression Ml tool for lncRNA prediction
- [DIAMOND](#diamond)
  - protein database homology search to eliminate candidates
- [fastp](#fastp)
  - template `sbatch` script called by `fastp_wrapper.sh`
  - wrapper script to submit all raw `fastq` files for trimming
- [FEELnc](#feelnc)
  - Filter: decision tree ensemble ML tool for lncRNA candidate identification
  - Classifier: edits gtf entries for predicted lncRNAs to include classifications
- [Intersect IDs](#intersect-ids)
  - Custom R script and sbatch script to combine outputs from up to four Ml tools and from DIAMOND
- [lncBoost](#lncboost)
  -  XGBoost model - iterative ML 'gradient boosting' model
- [lncFinder](#lncfinder)
  - support vector machine (SVM) geometric Ml tool - good for finding patterns in nonlinear data
  - image, r script, and sbatch script for this step
- [Minimap2](#minimap2)
  - mapping predictions against truth sequences for pipeline analysis
- [NCBI datasets](#ncbi-datasets)
  - clean up `datasets` reference assembly and annotation files for use in the pipeline
- [Prefetch array](#prefetch-array)
  - uses sra accessions list to run `prefetch` and `fasterq-dump` as an array
- [Genome generate](#genome-generate)
  - Indexes the reference genome
- [Stringtie](#stringtie)
  - Assemble: creates gtf for each fastq file
  - Merge: merges individual gtf files into one *de novo* transcriptome for the species
---
### Alignment array
INPUTS
- STAR indexed genome
- trimmed reads
- SRA accession list

COMMAND
```bash
STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX_DIR \
    --readFilesIn $READS \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUT_DIR}/${ID}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrIHstart 0 \
    --outSAMattributes NH HI AS nM NM \
    --outSAMstrandField intronMotif
```
FLAGS
 Flags
- outSAMattrIHstart 0
  - required for compatibility with certain older tools or specific pipelines that expect 0-based indexing for multi-mappers
- outSAMattributes
  - NH: Number of reported alignments for the query
  - HI: Query hit index (distinguishes multiple alignments of the same read)
  - AS: Alignment score
  - nM: Number of mismatches
  - NM: Edit distance to the reference (including insertions and deletions).
- outSAMstrandField
  - intronMotif is used for strand-specific RNA-seq data, which is common in eukaryotic transcriptomics
  - It helps STAR to correctly assign reads to the appropriate strand based on the presence of canonical splice motifs (GT/AG, GC/AG, AT/AC) at intron-exon boundaries
  - This is crucial for accurate gene expression quantification and downstream analyses.

SLURM NOTES
- cpus-per-task=8
- mem=64G
- time=07:00:00
  - based on 160 sra accessions
---
### Bam indexing
INPUTS
- STAR indexed genome
- trimmed reads
- SRA accession list

COMMAND
```bash
samtools index ${OUT_DIR}/${ID}_Aligned.sortedByCoord.out.bam
```

NOTES
- very light on resources
  - 1 cpu, 4 G, 1 hour
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
### CPAT plant
Model: Logistic regression
- maps features to a sigmoid curve to output a probability between 0 and 1
- simplest model while still considered ML

INPUTS:
1. Image
2. fasta file of candidate sequences
3. hexamer table from Plant-LncPipe github: 
  - `/base/dir/Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv`
4. logistic regression R data from Plant-LncPipe github
  - `/base/dir/Plant-LncRNA-pipeline-v2/Model/Plant.logit.RData`
5. output file

COMMAND:
```bash
singularity exec $IMAGE cpat.py -x $HEX -d $LOGIT -g $SEQS -o $OUT_ALL
```

FLAGS:
plain and simple. figure it out.

SLURM NOTES:
- 32 G
- 16 cpus-per-task
- 2 hours for soy
---
### DIAMOND
Double checks that any candidate transcripts do not match any known protein coding transcripts. See main `README.md` for database set-up.

Database path on hpcc: \
`/mnt/scratch/lemasjoh/BASF/uniprot_db/uniprot_out.dmnd`

INPUTS
1. candidate fasta sequences
2. output file path
3. database path

COMMAND
```bash
diamond blastx --threads $SLURM_CPUS_PER_TASK -d $DB \
	-q $INPUT -o $OUTPUT --outfmt 6 --sensitive \
	--max-target-seqs 1 --evalue 1e-5
```
FLAGS
- outfmt 6
  - creates a tsv file instead of a long human readable text file
- sensitive
  - diamond uses a 'seed' method to find matches
  - it is fast but might miss distant homologs
  - this flags tells it to use a more complex seeding strategy to find matches that are more divergent
  - if a transcript is acutally a degraded protein coding gene or other type of artefact this will hopefully catch it
- match-target-seqs 1
  - only reports the single best hit from matched in the db
	reduces output file size
- evalue 1e-5
  - e value: expectation value- is this hit just a random alignment or an actual hit? Lower expectation, lower value
  - currently set to a moderately strict threshold 

SLURM NOTES
- 32 cpus-per-task
- 64 G
- 2 hours for soy (works surprisingly fast- like, minutes)
---
### fastp
Cleans up raw `fastq` files from `prefetch_array.sb`. Trimmed reads are placed in the `species/03_outputs/01_fastp` directory. Double check that the array field in the SLURM parameters correctly reflects the number of accessions, and that the file paths for the species in question are correct.

INPUTS
- Raw `fastq` files from `prefetch_array.sb`
- SRA accession list
- Output directory for trimmed reads

COMMAND
```bash
sbatch fastp.sb
```
FLAGS

NOTES
Updated the original wrapper and template scripts to use an array instead. Easier on SLURM to organize the 25 jobs on its own without me submitting all of them at once. This script identified paired-end or singl end reads and selects fastp parameters accordingly. Works well when pulling lots of random SRA accessions from NCBI without having to filter for paired-end reads. This will be overkill in future runs where I am more selective of input data.

---
### FEELnc
Model: Decision tree ensemble
- builds hundreds of decision trees. each tree looks at a random subset of data and votes on whether a transcript is coding or noncoding
- the final result is the majority vote of all trees

INPUTS:
1. Image path (if using image)
2. Input gtf file: novel transcript gtfs from alignment/stringtie
3. reference annotation file
4. output file designation
5. gitcloned FEELnc directory

COMMAND:
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    env PERL5LIB="${FEELNC_DIR}/lib" \
    perl ${FEELNC_DIR}/scripts/FEELnc_filter.pl \
    -i $INPUT_GTF \
    -a $REF_ANN \
    --monoex=-1 -s 200 -p 20 > $OUTPUT_GTF
```
FLAGS:
- B
  - bind the working directory into the singularity environment; avoids file not found errors
- env
  - makes sure to call the correct perl library; essential for hpcc
- perl
  - the perl script from the gitcloned repo
- monoex=-1
  - monoexonic transcripts are usually noise in rnaseq data. most pipelines just throw them away
  - -1 means keep any monoexonic transcripts that are antisense to an mRNA
  - monoexonic antisense transcripts are much more likely to be a real noncoding transcript than a random monoexonic transcript out in intergenic space
- s 200
  - sets a transcripts size threashold to 200; a common value in lncRNA papers
- p 20
  - parallel threads to use

SLURM NOTES
- 64 G memory
- 20 cpus-per-task
- 4 hours for the 98,984 soy candidates
---
### Intersect IDs
Two scripts for this step: `intersect_ids.R` & `intersect_ids.sb`. The R script outputs a list of high confidence lncRNA IDs from `StringTie` and a table of concensus information for every submitted candidate.

INPUTS \
Designed to accept at least two outputs from any of the four tools. Fewer pipeline tools from the entire ensemble can be run if desired, or if a tool/model is not working. 

- FEELnc output
  - accepts either the raw `gtf` file or an indexed list of IDs
- `CPAT` output in `tsv` format
- `lncBoost` output in `tsv` format
- `lncFinder` output in `csv` format
- `DIAMOND` output in `tsv` format
- output file prefix
  - default is `concensus`


COMMAND
```bash
singularity run $IMAGE -f $FEELnc -c $CPAT -l $lncFinder -b $lncBoost -d $DIAMOND -o $OUTPFX
```
FLAGS

NOTES
4 G mem, 1 hour, 4 threads

---
### lncBoost
Model: XGBoost
- the iterative ML
- 'Gradient boosting' model
  - like random forest, but instead of making them all at once it makes them one by one
  - each new tree is designed to fix the 'errors' of the previous trees
  - currently one of the most powerful 'tabular' Ml tools in existence

INPUTS:
1. Image
2. path to script directory from Plant-LncPipe github
  - `/base/dir/Plant-LncRNA-pipeline-v2/Script`
3. path to model from Plant-LncPipe github
  - `/base/dir/Plant-LncRNA-pipeline-v2/Model/PlantLncBoost_model.cb`
4. fasta file of candidate sequences
5. output directory
  - makes two files. where do you want them
6. output features.csv file
7. output predictions.txt file
  - uses features to make predictions

COMMANDS:
1. feature extraction
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    python3 ${SCRIPT_DIR}/Feature_extraction.py -i $INPUT_FASTA -o $FEATURE_CSV
```

FLAGS:
`${SCRIPT_DIR}` is the path to script directory from Plant-LncPipe github

2. prediction
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    python3 ${SCRIPT_DIR}/PlantLncBoost_prediction.py \
    -i $FEATURE_CSV \
    -m $MODEL \
    -o $FINAL_RESULTS \
    -t 0.5
```
FLAGS:
`${SCRIPT_DIR}`
	path to script directory from Plant-LncPipe github
- i
  - input fasta
- m 
  - path to model from Plant-LncPipe github
- o
  - path to predictions.txt
- t 0.5
  - sets the threshold
  - lncBoost doesnt just 1 or 0 predictions, it assigns probability score between 0 and 1
	  - 0.5 means if 50% probability or higher, call it a lncRNA
  - standard, neutral decision for binary classification
	ok in this case because of the ensemble of Ml tools used in this pipeline
  - dont want one tool to filter out possibilities that other tools thought were lncRNAs

SLURM NOTES
- 20 cpus-per-task
- 64 G
- 4 hours for soy
---
### lncFinder
Model: support vector machine (SVM)
- the geometric ML
- plots transcripts as dots in high-dimensional space
- attempts to find the hyperplane (perfect boundary line) that creates the widest possible gap between the coding dots and noncoding dots
- good at finding patterns in nonlinear data

SIDE NOTE! different than other models. download from CRAN repo and run it through R
- The authors of Plant-LncPipe have an R script in the Scripts directory that runs lncFinder
- I didn't like it so Gemini and I rewrote it to work with the custom lncFinder.sif image I wanted to use
- much easier in hpcc not to have to deal with packages and paths and crap (as I quickly discovered)
- also, as I eventually discovered through much agony, the model that the Plant-LncPipe authors provide for lncFinder is broken/disfunctional
- luckily, lncFinder has a built in model that was trained on wheat data. use this model instead
	
INPUTS:
1. input fasta sequences
2. path to output csv file

PROCESS:
- edit the file paths in the lncFinder.R script before running the image command
  - maybe I should change this so that running the singularity command will feed R the fasta seqs and output path

COMMAND (in the R script):
```{r}
results <- lnc_finder(
  Seqs, 
  SS.features = FALSE, 
  format = "DNA", 
  frequencies.file = "wheat", 
  svm.model = "wheat", 
  parallel.cores = 32
)
```

FLAGS:
- Seqs
  - input fasta file
- SS.features = FALSE
  - lncFinder can predict lncRNAs after folding! only considers base pairs available for sequence interactions 
  - potentially more accurate, but way more computationally expensive
  - set to TRUE to have lncFinder pass the folding command onto another package
	folding uses the ViennaRNA package (RNAfold)
    - not in the image at this point! I'd have to edit the lncFinder.sif image to include ViennaRNA
    - it took ~24 hours just to fold all of the 156,364 corn candidates when I tried to use this
  - I need to rerun this step with the folded sequences that I generated to see if it works and compare results
- format = "DNA"
  - denotes fasta sequences
  - use "SS" if using already folded secondary structures
    - if SS.features = TRUE doesn't work well (i had some issues) run the Vienna package directly
    - use these as input instead of fasta seqs and set format to "SS"
- frequencies.file = "wheat"
  - the built-in wheat frequencies
- svm.model = "wheat"
  - the built-in wheat model
- parallel.cores = 32
  - uses parallel package (present in image) to utilize parallel threading
  - make sure this matches what is set in #SBATCH --cpus-per-task for SLURM

COMMAND (for the image itself)
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    Rscript $RSCRIPT
```
FLAGS:
- IMAGE
  - path to lncFinder.sif
- RSCRIPT
  - path to edited R script

SLURM NOTES
- 32 cpus-per-task
- 64 G
- 1 hour
  - HUGE NOTE: if having lncFinder fold seqs you'll need WAY more time
---
### Minimap2
INPUTS
- Predicted lncRNA sequences
- condensed 'truth' dataset sequences
- Output file pash

COMMAND
```bash
minimap2 -x cdna -t 8 $TRUTH $PRED -o $OUT_FILE
```
FLAGS
- -x
  - sets the input sequence type

NOTES
Super fast! 32 G, 8 threads, 30 minutes
---
### NCBI datasets
INPUTS
- unzipped `ncbi_dataset` directory
- species designation (ie z.mays)

COMMAND
```bash
# Unzip the dataset
unzip ncbi_dataset.zip

# Move into the uncompressed directory and run the cleanup script
cd ncbi_dataset
python3 IWGC_lncRNA_pipeline/scripts/ncbi_datasets_cleanup.py [species]

# Move the new clean assembly and annotation files into the inputs directory
mv species.chromosomes.fa IWGC_lncRNA_pipeline/species/01_inputs/
mv species_mRNA.gtf IWGC_lncRNA_pipeline/species/01_inputs/
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

---
### Genome generate
INPUTS
- Reference genome
  - `species/01_inputs/species.chromosomes.fa`
- Reference annotation 
  - `species/01_inputs/species_mRNA.gtf
- Index directory path

COMMAND
```bash
STAR \
     --runMode genomeGenerate \
     --runThreadN $SLURM_CPUS_PER_TASK \
     --genomeDir $INDEX \
     --genomeFastaFiles $GENOME \
     --sjdbGTFfile $ANNO \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 13
```
FLAGS
- -runMode genomeGenerate
- -sjdbGTFtagExonParentTranscript Parent
- -sjdbOverhang 99
- -genomeSAindexNbases 13
NOTES

---
### Stringtie
stringtie transcriptome assembly and lncRNA candidate discovery

two steps: assemble and merge.

assemble creates a gtf for each individual fastq file (also helpful to use an array job; see stringtie_assemble.sb) \
merge creates a single consolidated gtf file with input from each individual gtf. useful to make an array to loop over so that all the files get added to the path (see stringtie_merge.sb)

Once the merge step is complete the output file is a gtf file with all of the locations of aligned sequences and theyre associated annotations (assuming there is a match in the reference genome)

To obtain candidate lncRNAs, we grep for all of the aligned sequences/locations that did not overlap already annotated protein coding regions:
```bash
# grep for novel transcript alignments (maybe not necessarily novel, just not in the given reference)
grep 'transcript_id "MSTRG' merge.gtf > candidates.gtf
# write out the sequences for each novel transcript
gffread -w candidates.fasta -g genome.fasta candidates.gtf
# create a list with novel alignment IDs for use later
grep '>' candidates.fasta | sed 's/>//g' | sort -u > candidates.txt
```
These three steps set up the pipeline and result in the unfiltered files used as inputs for each model
---