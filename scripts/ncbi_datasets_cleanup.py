#!/usr/bin/env python3
import os
import glob
import sys
import re

# 1. Capture the species name from the command-line argument
if len(sys.argv) < 2:
    print("Error: Missing species name argument!")
    print("Usage: python3 ncbi_datasets_cleanup.py [species_name]")
    print("Example: python3 ncbi_datasets_cleanup.py g.max")
    sys.exit(1)

species_prefix = sys.argv[1]

print("Searching for NCBI dataset directory structure...")

# 2. Find the local GCF dataset directory
gcf_dirs = glob.glob("data/GCF*")

if not gcf_dirs:
    print("Error: Could not find any directory matching 'data/GCF*'!")
    print(f"Current working directory is: {os.getcwd()}")
    sys.exit(1)

target_dir = gcf_dirs[0]
print(f"Found active genome directory: {target_dir}")

# 3. Dynamically resolve file paths inside that specific GCF folder
fasta_matches = glob.glob(os.path.join(target_dir, "*_genomic.fna"))
input_fasta = fasta_matches[0] if fasta_matches else None
input_gff = os.path.join(target_dir, "genomic.gff")

# Safety check
missing_files = []
if not input_fasta: missing_files.append("*_genomic.fna")
if not os.path.exists(input_gff): missing_files.append("genomic.gff")

if missing_files:
    print(f"Error: Missing required files in {target_dir}: {', '.join(missing_files)}")
    sys.exit(1)

# Define ALL custom output paths
output_fasta = os.path.join(target_dir, f"{species_prefix}.chromosomes.fa")
output_full_gff3 = os.path.join(target_dir, f"{species_prefix}.chromosomes.gff3")
output_mrna_gff = os.path.join(target_dir, f"{species_prefix}_mRNA.gff")

print(f"-> Found FASTA:           {os.path.basename(input_fasta)}")
print(f"-> Found GFF:             {os.path.basename(input_gff)}")
print(f"-> Target FASTA Output:   {output_fasta}")
print(f"-> Target Full GFF3:      {output_full_gff3}")
print(f"-> Target mRNA GFF:       {output_mrna_gff}")

# --- Step 1: Extract Chromosome Map directly from FASTA headers ---
print("\nStep 1: Parsing FASTA headers to map primary chromosomes...")
chr_mapping = {}
chrom_regex = re.compile(r"chromosome\s+([A-Za-z0-9]+)")

with open(input_fasta, 'r') as infile:
    for line in infile:
        if line.startswith('>'):
            header = line.strip()
            if "unlocalized" in header.lower() or "unplaced" in header.lower() or "mitochondrion" in header.lower() or "chloroplast" in header.lower():
                continue
                
            match = chrom_regex.search(header)
            if match:
                accession = header.split()[0].replace('>', '')
                chrom_num = match.group(1)
                chr_mapping[accession] = f"chr{chrom_num}"

if not chr_mapping:
    print("Error: No primary chromosomes could be parsed from the FASTA headers!")
    sys.exit(1)
else:
    print(f"Success! Detected {len(chr_mapping)} main chromosomes.")

# --- Step 2: Clean FASTA ---
print("\nStep 2: Filtering and renaming FASTA file...")
with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
    write_sequence = False
    for line in infile:
        if line.startswith('>'):
            accession = line.split()[0].replace('>', '')
            if accession in chr_mapping:
                outfile.write(f">{chr_mapping[accession]}\n")
                write_sequence = True
            else:
                write_sequence = False
        else:
            if write_sequence:
                outfile.write(line)
print(f"Saved: {os.path.basename(output_fasta)}")

# --- Step 3: Clean GFF Files ---
print("\nStep 3: Filtering GFF annotations and writing both GFF3 and structural mRNA outputs...")

with open(input_gff, 'r') as infile, \
     open(output_full_gff3, 'w') as out_full, \
     open(output_mrna_gff, 'w') as out_mrna:
     
    # Ensure standard GFF3 header is written to the clean master file
    out_full.write("##gff-version 3\n")
    out_mrna.write("##gff-version 3\n")

    for line in infile:
        # Skip original comment headers since we explicitly wrote the crisp gff3 line above
        if line.startswith('#'):
            continue
            
        columns = line.split('\t')
        seqid = columns[0]
        feature_type = columns[2]
        
        # Check if feature belongs to our primary chromosomes
        if seqid in chr_mapping:
            # Update sequence ID naming scheme to 'chrX'
            columns[0] = chr_mapping[seqid]
            updated_line = '\t'.join(columns)
            
            # File 1: Write EVERYTHING (Full clean GFF3 reference for standard alignments)
            out_full.write(updated_line)
            
            # File 2: Write transcript tracking structures only (Specialized file for your downstream step)
            if feature_type in ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'lnc_RNA']:
                out_mrna.write(updated_line)

print(f"Saved: {os.path.basename(output_full_gff3)}")
print(f"Saved: {os.path.basename(output_mrna_gff)}")
print(f"\nProcessing complete! All files generated under prefix ID: '{species_prefix}'")