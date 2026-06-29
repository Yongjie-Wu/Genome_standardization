#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yongjie Wu
# @FileName: GWH_genome.py
# @QQ:1776262486 WeChat:RyougiSh1ki_0217
import sys, re, argparse, os
import subprocess

class FastaReader:
    def __init__(self, file_path):
        self.file_path = file_path
        self.sequences = {}
        self.names = []
        self._read_fasta()

    def _read_fasta(self):
        current_name = ""
        current_sequence = []

        try:
            with open(self.file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_name and current_sequence:
                            self.sequences[current_name] = ''.join(current_sequence)
                            self.names.append(current_name)
                        current_name = line[1:].strip()
                        current_sequence = []
                    elif line:
                        current_sequence.append(line)

                if current_name and current_sequence:
                    self.sequences[current_name] = ''.join(current_sequence)
                    self.names.append(current_name)

        except FileNotFoundError:
            print(f"Error: file not found '{self.file_path}'")
        except Exception as e:
            print(f"Error occurred while reading the file: {str(e)}")

    @property
    def sequence_names(self):
        return self.names

    def get_sequence(self, name):
        return self.sequences.get(name)

    def __len__(self):
        return len(self.sequences)


def read_chr(chrfile):
    """Read the chromosome FASTA file and build a rawID-to-originalID mapping."""
    chrdict = {}
    seqdict = {}
    chrfasta = FastaReader(chrfile)

    for i in chrfasta.sequence_names:
        groups = re.search(r"(\S+?)\s+.*?OriSeqID=(\S+?)\s+.*", i)
        if groups:
            rawID = groups.group(1)
            orignalID = groups.group(2)
            chrdict[rawID] = orignalID
            seqdict[orignalID] = chrfasta.get_sequence(i)

    return chrdict, seqdict


def parse_gff_mappings(gfffile):
    """Parse the GFF file and return mappings of mRNA Accession to ID and protein Accession to ID."""
    mrna_acc_to_id = {}
    prot_acc_to_id = {}

    with open(gfffile, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            attributes = fields[8]

            if feature_type == 'mRNA':
                id_match = re.search(r'ID=([^;]+)', attributes)
                acc_match = re.search(r'Accession=([^;\s]+)', attributes)
                if id_match and acc_match:
                    mrna_id = id_match.group(1)
                    mrna_acc = acc_match.group(1)
                    mrna_acc_to_id[mrna_acc] = mrna_id

            elif feature_type == 'CDS':
                prot_acc_match = re.search(r'Protein_Accession=([^;\s]+)', attributes)
                parent_match = re.search(r'Parent=([^;]+)', attributes)
                if prot_acc_match and parent_match:
                    prot_acc = prot_acc_match.group(1)
                    parent_id = parent_match.group(1)  # mRNA ID
                    prot_acc_to_id[prot_acc] = parent_id

    return mrna_acc_to_id, prot_acc_to_id


def read_gff_to_bed(gfffile, output, chrdict, name_prefix):
    """Read the GFF file, extract mRNA features, and convert them to 6-column BED format using the mRNA ID field."""

    # Extract the Accession-to-Chromosome mapping from comment lines
    accession_to_chromosome = {}
    with open(gfffile, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if 'OriSeqID=' in line and 'Accession=' in line:
                    parts = line[1:].strip().split()
                    oriseqid = None
                    accession = None
                    for part in parts:
                        if part.startswith('OriSeqID='):
                            oriseqid = part.split('=', 1)[1]
                        elif part.startswith('Accession='):
                            accession = part.split('=', 1)[1]
                    if oriseqid and accession:
                        accession_to_chromosome[accession] = oriseqid

    with open(gfffile, 'r') as file, open(output, 'w') as f2:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            seq_id = fields[0]
            feature_type = fields[2]

            if feature_type == 'mRNA':
                attributes = fields[8]
                id_match = re.search(r'ID=([^;]+)', attributes)
                if not id_match:
                    continue
                mrna_id = id_match.group(1)

                start = fields[3]
                end = fields[4]
                strand = fields[6]

                # Determine the chromosome ID
                chromosome_id = None
                if seq_id in accession_to_chromosome:
                    chromosome_id = accession_to_chromosome[seq_id]
                elif seq_id in chrdict:
                    chromosome_id = chrdict[seq_id]
                else:
                    chromosome_id = seq_id

                prefixed_id = f"{name_prefix}_{mrna_id}"
                score = "0"
                bed_start = str(int(start) - 1)
                output_line = f"{chromosome_id}\t{bed_start}\t{end}\t{prefixed_id}\t{score}\t{strand}\n"
                f2.write(output_line)


def convert_gff_chromosome(gfffile, output, chrdict):
    """Convert the GFF file by updating only the chromosome ID in the first column, replacing Accession with OriSeqID."""
    accession_to_chromosome = {}
    with open(gfffile, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if 'OriSeqID=' in line and 'Accession=' in line:
                    parts = line[1:].strip().split()
                    oriseqid = None
                    accession = None
                    for part in parts:
                        if part.startswith('OriSeqID='):
                            oriseqid = part.split('=', 1)[1]
                        elif part.startswith('Accession='):
                            accession = part.split('=', 1)[1]
                    if oriseqid and accession:
                        accession_to_chromosome[accession] = oriseqid

    with open(gfffile, 'r') as file, open(output, 'w') as f2:
        for line in file:
            line = line.rstrip('\n')
            if line.startswith('#'):
                f2.write(line + '\n')
                continue
            fields = line.split('\t')
            if len(fields) >= 9:
                seq_id = fields[0]
                new_seq_id = None
                if seq_id in accession_to_chromosome:
                    new_seq_id = accession_to_chromosome[seq_id]
                elif seq_id in chrdict:
                    new_seq_id = chrdict[seq_id]
                if new_seq_id:
                    fields[0] = new_seq_id
                new_line = '\t'.join(fields)
                f2.write(new_line + '\n')
            else:
                f2.write(line + '\n')


def process_sequence_file(input_file, output_file, name_prefix, acc_to_id_map, is_pep=False):
    """
    Process a CDS or protein FASTA file: read FASTA records, replace each Accession with the corresponding ID, and add the specified prefix to the output ID.
    For protein files (is_pep=True), first replace 'GWHP' with 'GWHT' in the header, and then extract the Accession for ID mapping.
    """
    with open(output_file, 'w') as f_out:
        with open(input_file, 'r') as f_in:
            current_header = ""
            current_seq = []
            
            for line in f_in:
                line = line.strip()
                if line.startswith('>'):
                    # Process the previous sequence record
                    if current_header and current_seq:
                        # Process the header
                        if is_pep:
                            # Protein FASTA: first replace GWHP with GWHT
                            modified_header = current_header.replace('GWHP', 'GWHT')
                            # Extract the Accession, i.e., the first token of the modified header
                            acc = modified_header.split()[0].strip()
                        else:
                            # CDS FASTA: directly extract the Accession
                            acc = current_header.split()[0].strip()
                        
                        # Retrieve the mapped ID
                        gene_id = acc_to_id_map.get(acc, acc)
                        prefixed_id = f"{name_prefix}_{gene_id}"
                        
                        # Write the sequence record
                        f_out.write(f">{prefixed_id}\n")
                        f_out.write('\n'.join(current_seq) + '\n')
                    
                    # Start a new sequence record
                    current_header = line[1:].strip()
                    current_seq = []
                else:
                    if line:  # Add non-empty lines only
                        current_seq.append(line)
            
            # Process the last sequence record
            if current_header and current_seq:
                if is_pep:
                    modified_header = current_header.replace('GWHP', 'GWHT')
                    acc = modified_header.split()[0].strip()
                else:
                    acc = current_header.split()[0].strip()
                
                gene_id = acc_to_id_map.get(acc, acc)
                prefixed_id = f"{name_prefix}_{gene_id}"
                
                f_out.write(f">{prefixed_id}\n")
                f_out.write('\n'.join(current_seq) + '\n')


def generate_cds_pep_from_gff(gff_file, genome_file, file_prefix, name_prefix, mrna_acc_to_id):
    """
    Generate CDS and protein FASTA files from the GFF and genome FASTA files using gffread and seqkit.
    """
    import tempfile
    import os
    
    print("  Generating CDS file with gffread...")
    
    # Create temporary files
    temp_cds = tempfile.NamedTemporaryFile(delete=False, suffix='.cds', mode='w')
    temp_cds_name = temp_cds.name
    temp_cds.close()
    
    temp_prot_raw = tempfile.NamedTemporaryFile(delete=False, suffix='.prot', mode='w')
    temp_prot_raw_name = temp_prot_raw.name
    temp_prot_raw.close()
    
    temp_pep = tempfile.NamedTemporaryFile(delete=False, suffix='.pep', mode='w')
    temp_pep_name = temp_pep.name
    temp_pep.close()
    
    # Step 1: Generate the CDS file with gffread
    cmd_gffread = f"gffread -x {temp_cds_name} -g {genome_file} {gff_file}"
    print(f"  Running: {cmd_gffread}")
    
    result = subprocess.run(cmd_gffread, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Warning: gffread failed: {result.stderr}")
        return None, None
    
    # Step 2: Translate CDS sequences into proteins with seqkit
    print("  Translating CDS sequences into proteins with seqkit...")
    cmd_seqkit = f"seqkit translate {temp_cds_name} -o {temp_prot_raw_name} -T 11 -F -x --clean"
    print(f"  Running: {cmd_seqkit}")
    
    result = subprocess.run(cmd_seqkit, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Warning: seqkit translation failed: {result.stderr}")
        os.unlink(temp_cds_name)
        return None, None
    
    # Step 3: Process the protein file and remove frame information
    print("  Processing the protein file...")
    cmd_sed1 = f"sed 's/_frame=[0-9]\\+//' {temp_prot_raw_name} > {temp_pep_name}"
    subprocess.run(cmd_sed1, shell=True, check=True)
    
    # Step 4: Process stop codons
    print("  Processing stop codons...")
    cds_output = f"{file_prefix}.cds"
    pep_output = f"{file_prefix}.protein"
    
    # Process the CDS file by replacing IDs
    print("  Processing CDS IDs...")
    process_sequence_file(temp_cds_name, cds_output, name_prefix, mrna_acc_to_id, is_pep=False)
    
    # Process the protein file: first replace stop codons, then process IDs
    print("  Processing stop codons and IDs in the protein file...")
    
    # Create a temporary file for stop-codon processing
    temp_pep_processed = tempfile.NamedTemporaryFile(delete=False, suffix='.processed.pep', mode='w')
    temp_pep_processed_name = temp_pep_processed.name
    temp_pep_processed.close()
    
    # Replace stop codons
    with open(temp_pep_name, 'r') as f_in, open(temp_pep_processed_name, 'w') as f_out:
        for line in f_in:
            if not line.startswith('>'):
                # Replace . and U with *
                line = line.replace('.', '*').replace('U', '*')
            f_out.write(line)
    
    # Process protein IDs
    process_sequence_file(temp_pep_processed_name, pep_output, name_prefix, mrna_acc_to_id, is_pep=True)
    
    # Remove temporary files
    os.unlink(temp_cds_name)
    os.unlink(temp_prot_raw_name)
    os.unlink(temp_pep_name)
    os.unlink(temp_pep_processed_name)
    
    return cds_output, pep_output


def main():
    parser = argparse.ArgumentParser(description='Process GFF and sequence files')
    parser.add_argument('--gff', required=True, help='Input GFF file')
    parser.add_argument('--chr', required=True, help='Input chromosome file')
    parser.add_argument('--cds', help='Input CDS file (optional)')
    parser.add_argument('--pep', help='Input peptide file (optional)')
    parser.add_argument('--file', required=True, help='Prefix for output file names')
    parser.add_argument('--name', required=True, help='Prefix for gene IDs in CDS, PEP and BED files')
    args = parser.parse_args()

    # 1. Process the chromosome FASTA file
    chrdict, seqdict = read_chr(args.chr)
    genome_output = f"{args.file}.genome.fasta"
    with open(genome_output, 'w') as f:
        for chr_id, seq in seqdict.items():
            f.write(f">{chr_id}\n{seq}\n")

    # 2. Parse the GFF file and obtain the Accession-to-ID mapping (mRNA Accession -> mRNA ID)
    mrna_acc_to_id, _ = parse_gff_mappings(args.gff)  # Protein mapping is not used here

    # 3. Generate the BED file using mRNA IDs
    bed_output = f"{args.file}.bed"
    read_gff_to_bed(args.gff, bed_output, chrdict, args.name)

    # 4. Convert chromosome IDs in the GFF file
    gff_output = f"{args.file}.gff"
    convert_gff_chromosome(args.gff, gff_output, chrdict)

    # Prepare the output file list
    output_files = [
        f"  Genome sequence: {genome_output}",
        f"  Gene positions (BED): {bed_output}",
        f"  Chromosome-converted GFF: {gff_output}"
    ]

    # 5. Process CDS and protein files
    if args.cds and args.pep:
        # The user provided both CDS and protein files
        print("\nUsing user-provided CDS and protein files...")
        cds_output = f"{args.file}.cds"
        process_sequence_file(args.cds, cds_output, args.name, mrna_acc_to_id, is_pep=False)
        output_files.append(f"  CDS sequences: {cds_output}")
        
        pep_output = f"{args.file}.protein"
        process_sequence_file(args.pep, pep_output, args.name, mrna_acc_to_id, is_pep=True)
        output_files.append(f"  Protein sequences: {pep_output}")
    
    elif args.cds and not args.pep:
        # Only a CDS file was provided; generate the protein file from CDS sequences
        print("\nA CDS file was provided; generating the protein file from CDS sequences...")
        cds_output = f"{args.file}.cds"
        process_sequence_file(args.cds, cds_output, args.name, mrna_acc_to_id, is_pep=False)
        output_files.append(f"  CDS sequences: {cds_output}")
        
        print("  Translating CDS sequences to generate the protein file...")
        # Translate CDS sequences with seqkit
        temp_prot = f"{args.file}.temp.protein"
        pep_output = f"{args.file}.protein"
        
        # Translate sequences
        cmd_seqkit = f"seqkit translate {cds_output} -o {temp_prot} -T 11 -F -x --clean"
        subprocess.run(cmd_seqkit, shell=True, check=True)
        
        # Remove frame information
        temp_pep = f"{args.file}.temp.pep"
        cmd_sed1 = f"sed 's/_frame=[0-9]\\+//' {temp_prot} > {temp_pep}"
        subprocess.run(cmd_sed1, shell=True, check=True)
        
        # Process stop codons
        with open(temp_pep, 'r') as f_in, open(pep_output, 'w') as f_out:
            for line in f_in:
                if not line.startswith('>'):
                    line = line.replace('.', '*').replace('U', '*')
                f_out.write(line)
        
        # Remove temporary files
        os.unlink(temp_prot)
        os.unlink(temp_pep)
        output_files.append(f"  Protein sequences: {pep_output}")
    
    elif not args.cds and args.pep:
        # Only a protein file was provided, but a CDS file is required; reverse generation is not supported
        print("\nError: a protein file was provided without a CDS file, which cannot be processed. Please provide both CDS and protein files, or provide neither.")
        sys.exit(1)
    
    else:
        # Neither CDS nor protein files were provided; generate them from the GFF and genome FASTA files
        print("\nNo CDS or protein files were provided; generating them from the GFF and genome FASTA files...")
        cds_output, pep_output = generate_cds_pep_from_gff(args.gff, args.chr, args.file, args.name, mrna_acc_to_id)
        
        if cds_output and pep_output:
            output_files.append(f"  CDS sequences: {cds_output}")
            output_files.append(f"  Protein sequences: {pep_output}")
        else:
            print("  Warning: failed to generate CDS and protein files. Please check whether gffread and seqkit are installed.")
            print("  You can install them manually with the following command:")
            print("    conda install -c bioconda gffread seqkit")

    # Print the final results
    print("\n" + "=" * 60)
    print("Processing completed. Output files:")
    for file_info in output_files:
        print(file_info)
    print("=" * 60)


if __name__ == "__main__":
    main()
