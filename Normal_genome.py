#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yongjie Wu
# @FileName: Normal_genome.py
# @QQ:1776262486 WeChat:RyougiSh1ki_0217
import sys
import re
import argparse
import subprocess
import os
from pathlib import Path


def read_fasta_simple(file_path):
    """Read a FASTA file and return a dictionary {id: sequence}, preserving the original IDs."""
    sequences = {}
    current_id = None
    current_seq = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]  # Use only the part before the first whitespace as the sequence ID
                    current_seq = []
                elif line:
                    current_seq.append(line)
            if current_id:
                sequences[current_id] = ''.join(current_seq)
    except Exception as e:
        print(f"Error: failed to read file {file_path} - {e}", file=sys.stderr)
        sys.exit(1)
    return sequences


def write_fasta_with_prefix(sequences, out_file, name_prefix=None):
    """Write a sequence dictionary to a FASTA file, with an optional ID prefix."""
    with open(out_file, 'w') as f:
        for seq_id, seq in sequences.items():
            if name_prefix:
                # Check whether the prefix is already present to avoid duplication
                if seq_id.startswith(name_prefix + '_'):
                    new_id = seq_id
                else:
                    new_id = f"{name_prefix}_{seq_id}"
            else:
                new_id = seq_id
            f.write(f">{new_id}\n")
            # Wrap sequences at 60 characters per line
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def remove_gff_comments(gff_in, gff_out):
    """Remove all comment lines starting with # from a GFF file and retain data lines."""
    with open(gff_in, 'r') as fin, open(gff_out, 'w') as fout:
        for line in fin:
            if not line.startswith('#'):
                fout.write(line)


def generate_bed(gff_in, bed_out, feature_type, note_attr, name_prefix):
    """Extract the specified feature type from a GFF file to generate a BED file, adding the name prefix to gene IDs."""
    with open(gff_in, 'r') as fin, open(bed_out, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != feature_type:
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # BED uses 0-based coordinates
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]

            # Parse the attributes in the ninth column
            attr_dict = {}
            for item in attributes.split(';'):
                item = item.strip()
                if '=' in item:
                    key, val = item.split('=', 1)
                    attr_dict[key] = val

            # Retrieve the specified note attribute value
            if note_attr in attr_dict:
                gene_id = attr_dict[note_attr]
                # Add the name prefix
                name = f"{name_prefix}_{gene_id}"
            else:
                name = "."

            # BED line: chrom, start, end, name, score, strand
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


def run_command(cmd, description):
    """Run a system command and check the return code."""
    print(f"  Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Error: {description} failed", file=sys.stderr)
        print(f"  stderr: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result


def generate_cds_with_gffread(gff_file, chr_file, output_file):
    """Generate a CDS file using gffread."""
    cmd = f"gffread -x {output_file} -g {chr_file} {gff_file}"
    run_command(cmd, "gffread")
    
    # Check whether the file was generated successfully
    if not Path(output_file).exists():
        print(f"  Error: gffread failed to generate {output_file}", file=sys.stderr)
        sys.exit(1)
    
    # Count the number of sequences
    count = 0
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def generate_protein_with_seqkit(cds_file, output_file, name_prefix):
    """Generate a protein file using seqkit translate and directly add the prefix."""
    temp_file = output_file + ".temp"
    
    # Step 1: seqkit translate
    cmd1 = f"seqkit translate {cds_file} -o {temp_file} -T 11 -F -x --clean"
    run_command(cmd1, "seqkit translate")
    
    # Step 2: remove the _frame= suffix and add the prefix
    with open(temp_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                # Extract the ID and remove the _frame= suffix
                old_id = line[1:].strip()
                # Remove the _frame= component
                clean_id = re.sub(r'_frame=\d+', '', old_id)
                # Add the prefix after checking whether it is already present
                if clean_id.startswith(name_prefix + '_'):
                    new_id = clean_id
                else:
                    new_id = f"{name_prefix}_{clean_id}"
                fout.write(f">{new_id}\n")
            else:
                fout.write(line)
    
    # Remove the temporary file
    os.remove(temp_file)
    
    # Count the number of sequences
    count = 0
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def show_help():
    """Display help information."""
    help_text = """
Usage:
  python Normal_genome.py \\
    --gff <input.gff> \\
    --chr <input.chr.fasta> \\
    [--cds <input.cds>] \\
    [--pep <input.protein>] \\
    --file <output_prefix> \\
    --feature <feature_type> \\
    --note <attribute> \\
    --name <name_prefix>

Required arguments:
  --gff               Input GFF file
  --chr               Input chromosome FASTA file
  --file              Output file prefix
  --feature           Feature type used for BED generation, e.g., mRNA or gene
  --note              Attribute in the ninth GFF column used as the gene ID, e.g., ID or Name
  --name              Prefix added before gene IDs

Optional arguments:
  --cds               Input CDS FASTA file; generated by gffread if not provided
  --pep               Input protein FASTA file; generated by seqkit if not provided
  --help              Display this help message

Required tools:
  gffread             Used to generate CDS files when --cds is not provided
  seqkit              Used to generate protein files when --pep is not provided

Example:
  python Normal_genome.py \\
    --gff DSlongest.gff \\
    --chr Dicot_Ros_Ros_Ros_Pyrus_bretschneideri_cv_Dangshansuli_AAAS_T2T.chr \\
    --file Pyr.br_cv.DS_AAAS \\
    --feature mRNA \\
    --note ID \\
    --name Magno_Ros_Ros_Ros_Pyr.br_cv.DS_AAAS
"""
    print(help_text)


def main():
    # Check whether any arguments were provided
    if len(sys.argv) == 1:
        print("Error: missing arguments!")
        show_help()
        sys.exit(1)

    # Check whether --help was specified
    if '--help' in sys.argv or '-h' in sys.argv:
        show_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(description='Process genome files in standard format', add_help=False)
    parser.add_argument('--gff', required=True, help='Input GFF file')
    parser.add_argument('--chr', required=True, help='Input chromosome FASTA file')
    parser.add_argument('--cds', help='Input CDS FASTA file; optional, automatically generated if not provided')
    parser.add_argument('--pep', help='Input protein FASTA file; optional, automatically generated if not provided')
    parser.add_argument('--file', required=True, help='Output file prefix')
    parser.add_argument('--feature', required=True, help='Feature type used for BED generation, e.g., gene or mRNA')
    parser.add_argument('--note', required=True, help='Attribute in the ninth GFF column used as the gene ID, e.g., ID or Name')
    parser.add_argument('--name', required=True, help='Prefix added before gene IDs')

    # Manually parse arguments to provide clearer error messages
    try:
        args = parser.parse_args()
    except SystemExit:
        # Display help when argument parsing fails
        show_help()
        sys.exit(1)

    # Process the chromosome file -> genome.fasta; chromosome IDs are not prefixed
    print("Processing chromosome file...")
    chr_seqs = read_fasta_simple(args.chr)
    genome_out = f"{args.file}.genome.fasta"
    write_fasta_with_prefix(chr_seqs, genome_out)  # Do not add a prefix

    # Generate the BED file; gene IDs are prefixed
    print("Generating BED file...")
    bed_out = f"{args.file}.bed"
    generate_bed(args.gff, bed_out, args.feature, args.note, args.name)

    # Remove GFF comment lines
    print("Cleaning GFF comments...")
    gff_out = f"{args.file}.gff"
    remove_gff_comments(args.gff, gff_out)

    # Prepare the output file list
    output_files = [
        f"  Genome sequence: {genome_out}",
        f"  Gene locations (BED): {bed_out}",
        f"  Cleaned GFF: {gff_out}"
    ]

    # Process the CDS file
    if args.cds:
        print(f"Using the provided CDS file: {args.cds}")
        cds_seqs = read_fasta_simple(args.cds)
        cds_out = f"{args.file}.cds"
        write_fasta_with_prefix(cds_seqs, cds_out, args.name)
    else:
        print("No CDS file provided; generating it with gffread...")
        cds_out = f"{args.file}.cds"
        # Generate the CDS file; IDs are not prefixed at this stage
        cds_count = generate_cds_with_gffread(args.gff, args.chr, cds_out)
        # Read the generated CDS file and add the prefix
        cds_seqs = read_fasta_simple(cds_out)
        write_fasta_with_prefix(cds_seqs, cds_out, args.name)
    output_files.append(f"  CDS sequence: {cds_out}")

    # Process the protein file
    if args.pep:
        print(f"Using the provided protein file: {args.pep}")
        pep_seqs = read_fasta_simple(args.pep)
        pep_out = f"{args.file}.protein"
        write_fasta_with_prefix(pep_seqs, pep_out, args.name)
    else:
        print("No protein file provided; generating it from CDS using seqkit...")
        pep_out = f"{args.file}.protein"
        # Directly generate a prefixed protein file
        generate_protein_with_seqkit(cds_out, pep_out, args.name)
    output_files.append(f"  Protein sequence: {pep_out}")

    # Print only the final results
    print("\n" + "=" * 60)
    print("Processing completed. Output files:")
    for f in output_files:
        print(f)
    print("=" * 60)


if __name__ == "__main__":
    main()
