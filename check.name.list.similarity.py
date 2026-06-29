#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gene ID Consistency Checker

This script validates:
1. Gene ID consistency between .bed, .protein, and .cds files
2. No duplicate gene IDs in any file
3. Chromosome ID consistency between .bed and .genome.fasta files
4. File formats are correct
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Set, List, Tuple
from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser

# ============================================================================
# Constants
# ============================================================================

MAX_DIFF_DISPLAY = 20
EXIT_SUCCESS = 0
EXIT_FILE_NOT_FOUND = 1
EXIT_VALIDATION_FAILED = 2


class Colors:
    """ANSI color codes for terminal output."""
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"

    @classmethod
    def green(cls, text: str) -> str:
        return f"{cls.GREEN}{text}{cls.ENDC}"

    @classmethod
    def red(cls, text: str) -> str:
        return f"{cls.RED}{text}{cls.ENDC}"

    @classmethod
    def yellow(cls, text: str) -> str:
        return f"{cls.YELLOW}{text}{cls.ENDC}"


def write_diff_to_file(diff_set: Set[str], file1: str, file2: str, output_dir: Path, prefix: str) -> Path:
    """Write diff to file and return the file path."""
    output_file = output_dir / f"{prefix}_diff_{file1}_vs_{file2}.txt"
    with open(output_file, 'w') as f:
        f.write(f"# Gene IDs in {file1} but not in {file2}\n")
        for gene_id in sorted(diff_set):
            f.write(f"{gene_id}\n")
    return output_file


def write_duplicates_to_file(duplicates: List[Tuple[str, int]], file_type: str, output_dir: Path, prefix: str) -> Path:
    """Write duplicates to file and return the file path."""
    output_file = output_dir / f"{prefix}_duplicates_{file_type}.txt"
    with open(output_file, 'w') as f:
        f.write(f"# Duplicate gene IDs in {file_type} file (count)\n")
        for gene_id, count in sorted(duplicates):
            f.write(f"{gene_id}\t{count}\n")
    return output_file


class FileValidator:
    """Validates consistency between genomics data files."""

    def __init__(self, abbreviation: str, base_path: Optional[Path] = None):
        self.abbreviation = abbreviation
        self.base_path = base_path or Path.cwd()
        self.diff_files = []  # Track generated diff files

    def _get_file_paths(self) -> dict[str, Path]:
        """Get all required file paths."""
        return {
            "bed": self.base_path / f"{self.abbreviation}.bed",
            "protein": self.base_path / f"{self.abbreviation}.protein",
            "cds": self.base_path / f"{self.abbreviation}.cds",
            "genome": self.base_path / f"{self.abbreviation}.genome.fasta",
        }

    def check_files_exist(self) -> tuple[bool, list[str]]:
        """Check if all required files exist."""
        paths = self._get_file_paths()
        missing = [name for name, path in paths.items() if not path.exists()]
        return len(missing) == 0, missing

    def _read_bed_ids(self, path: Path) -> Tuple[List[str], List[str], Set[str], List[Tuple[str, int]]]:
        """Read BED file and return gene IDs, chromosome IDs, and duplicates."""
        gene_ids = []
        chr_ids = []
        with open(path, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith(('#', 'track', 'browser')):
                    parts = line.split("\t")
                    if len(parts) >= 4:
                        chr_ids.append(parts[0])
                        gene_ids.append(parts[3])
        
        # Find duplicates
        counter = Counter(gene_ids)
        duplicates = [(item, count) for item, count in counter.items() if count > 1]
        
        return gene_ids, chr_ids, set(gene_ids), duplicates

    def _read_fasta_ids(self, path: Path) -> Tuple[List[str], Set[str], List[Tuple[str, int]]]:
        """Read FASTA file and return IDs and duplicates."""
        ids = []
        with open(path, "r") as f:
            for title, _ in SimpleFastaParser(f):
                ids.append(title.split()[0])
        
        # Find duplicates
        counter = Counter(ids)
        duplicates = [(item, count) for item, count in counter.items() if count > 1]
        
        return ids, set(ids), duplicates

    def validate(self) -> tuple[bool, dict]:
        """
        Run all validation checks.
        Returns: (is_valid, stats_dict)
        """
        paths = self._get_file_paths()
        stats = {}

        # Read BED file
        bed_gene_ids, bed_chr_ids, bed_set, bed_duplicates = self._read_bed_ids(paths["bed"])
        stats["bed_count"] = len(bed_gene_ids)
        stats["bed_duplicates"] = bed_duplicates

        # Read Protein file
        protein_ids, protein_set, protein_duplicates = self._read_fasta_ids(paths["protein"])
        stats["protein_count"] = len(protein_ids)
        stats["protein_duplicates"] = protein_duplicates

        # Read CDS file
        cds_ids, cds_set, cds_duplicates = self._read_fasta_ids(paths["cds"])
        stats["cds_count"] = len(cds_ids)
        stats["cds_duplicates"] = cds_duplicates

        # Read Genome file
        genome_ids, genome_set, genome_duplicates = self._read_fasta_ids(paths["genome"])
        stats["genome_count"] = len(genome_ids)
        stats["genome_duplicates"] = genome_duplicates
        stats["bed_chr_unique"] = len(set(bed_chr_ids))
        stats["genome_chr_unique"] = len(genome_set)
        
        # Check 1: Gene count consistency
        stats["gene_counts_match"] = (stats["bed_count"] == stats["protein_count"] == stats["cds_count"])
        
        # Check 2: Gene ID consistency between files
        stats["bed_vs_protein"] = bed_set == protein_set
        stats["bed_vs_cds"] = bed_set == cds_set
        stats["protein_vs_cds"] = protein_set == cds_set
        stats["gene_ids_match"] = (stats["bed_vs_protein"] and stats["bed_vs_cds"])
        
        # Calculate diffs
        stats["bed_only"] = bed_set - protein_set
        stats["protein_only"] = protein_set - bed_set
        stats["bed_vs_cds_only"] = bed_set - cds_set
        stats["cds_only"] = cds_set - bed_set
        
        # Check 3: No duplicate gene IDs
        stats["no_duplicates"] = (
            len(bed_duplicates) == 0 and
            len(protein_duplicates) == 0 and
            len(cds_duplicates) == 0
        )
        
        # Check 4: Chromosome ID consistency
        stats["chr_ids_match"] = (set(bed_chr_ids) == genome_set)
        stats["bed_chr_only"] = set(bed_chr_ids) - genome_set
        stats["genome_only"] = genome_set - set(bed_chr_ids)
        
        # Check 5: BED format (at least 4 columns)
        stats["bed_format_ok"] = True
        with open(paths["bed"], "r") as f:
            for line in f:
                if line.strip() and not line.startswith(('#', 'track', 'browser')):
                    if len(line.strip().split("\t")) < 4:
                        stats["bed_format_ok"] = False
                        break
        
        # Check 6: Protein sequences (no '.' or 'U')
        stats["protein_ok"] = True
        with open(paths["protein"], "r") as f:
            for _, seq in SimpleFastaParser(f):
                if "." in seq or "U" in seq:
                    stats["protein_ok"] = False
                    break

        is_valid = all([
            stats["gene_counts_match"],
            stats["gene_ids_match"],
            stats["no_duplicates"],
            stats["chr_ids_match"],
            stats["bed_format_ok"],
            stats["protein_ok"]
        ])
        
        return is_valid, stats


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Check gene ID consistency across .bed, .protein, .cds, and .genome.fasta files"
    )
    parser.add_argument("abbreviations", nargs="+", help="File prefix(es) to check")
    parser.add_argument("-d", "--directory", type=Path, default=Path.cwd(), help="Input directory")
    parser.add_argument("-o", "--output", type=Path, default=Path.cwd(), help="Output directory for diff files")
    args = parser.parse_args()

    exit_code = EXIT_SUCCESS

    for abbr in args.abbreviations:
        print(f"\n{'='*60}")
        print(f"Validating: {abbr}")
        print(f"{'='*60}")

        validator = FileValidator(abbr, args.directory)
        
        # Check files exist
        files_exist, missing = validator.check_files_exist()
        if not files_exist:
            print(Colors.red(f"✗ Missing files: {', '.join(missing)}"))
            exit_code = EXIT_FILE_NOT_FOUND
            continue

        # Run validation
        try:
            is_valid, stats = validator.validate()
            
            # Check gene counts
            if stats["gene_counts_match"]:
                print(Colors.green(f"✓ Line counts match: BED={stats['bed_count']}, Protein={stats['protein_count']}, CDS={stats['cds_count']}"))
            else:
                print(Colors.red(f"✗ Count mismatch: BED={stats['bed_count']}, Protein={stats['protein_count']}, CDS={stats['cds_count']}"))
                exit_code = EXIT_VALIDATION_FAILED
            
            # Check gene ID consistency between files
            if stats["gene_ids_match"]:
                print(Colors.green("✓ Gene IDs match between all files"))
            else:
                print(Colors.red("✗ Gene IDs don't match between files"))
                
                # Show specific mismatches
                if not stats["bed_vs_protein"]:
                    print(Colors.yellow(f"  BED vs Protein: BED has {len(stats['bed_only'])} unique, Protein has {len(stats['protein_only'])} unique"))
                    # Handle bed-only IDs
                    if len(stats["bed_only"]) <= MAX_DIFF_DISPLAY:
                        for gene_id in sorted(stats["bed_only"])[:MAX_DIFF_DISPLAY]:
                            print(Colors.yellow(f"    BED only: {gene_id}"))
                    else:
                        diff_file = write_diff_to_file(stats["bed_only"], "BED", "Protein", args.output, abbr)
                        print(Colors.yellow(f"    >20 unique IDs in BED, written to: {diff_file}"))
                    
                    # Handle protein-only IDs
                    if len(stats["protein_only"]) <= MAX_DIFF_DISPLAY:
                        for gene_id in sorted(stats["protein_only"])[:MAX_DIFF_DISPLAY]:
                            print(Colors.yellow(f"    Protein only: {gene_id}"))
                    else:
                        diff_file = write_diff_to_file(stats["protein_only"], "Protein", "BED", args.output, abbr)
                        print(Colors.yellow(f"    >20 unique IDs in Protein, written to: {diff_file}"))
                
                if not stats["bed_vs_cds"]:
                    print(Colors.yellow(f"  BED vs CDS: BED has {len(stats['bed_vs_cds_only'])} unique, CDS has {len(stats['cds_only'])} unique"))
                    # Handle BED vs CDS differences
                    if len(stats["bed_vs_cds_only"]) <= MAX_DIFF_DISPLAY:
                        for gene_id in sorted(stats["bed_vs_cds_only"])[:MAX_DIFF_DISPLAY]:
                            print(Colors.yellow(f"    BED only (vs CDS): {gene_id}"))
                    else:
                        diff_file = write_diff_to_file(stats["bed_vs_cds_only"], "BED", "CDS", args.output, abbr)
                        print(Colors.yellow(f"    >20 unique IDs in BED (vs CDS), written to: {diff_file}"))
                    
                    if len(stats["cds_only"]) <= MAX_DIFF_DISPLAY:
                        for gene_id in sorted(stats["cds_only"])[:MAX_DIFF_DISPLAY]:
                            print(Colors.yellow(f"    CDS only: {gene_id}"))
                    else:
                        diff_file = write_diff_to_file(stats["cds_only"], "CDS", "BED", args.output, abbr)
                        print(Colors.yellow(f"    >20 unique IDs in CDS, written to: {diff_file}"))
            
            # Check duplicates
            if stats["no_duplicates"]:
                print(Colors.green("✓ No duplicate gene IDs found"))
            else:
                print(Colors.red("✗ Duplicate gene IDs found"))
                
                # Check BED duplicates
                if stats["bed_duplicates"]:
                    if len(stats["bed_duplicates"]) <= MAX_DIFF_DISPLAY:
                        for gene_id, count in stats["bed_duplicates"]:
                            print(Colors.yellow(f"  BED duplicate: {gene_id} (x{count})"))
                    else:
                        diff_file = write_duplicates_to_file(stats["bed_duplicates"], "BED", args.output, abbr)
                        print(Colors.yellow(f"  >20 duplicate IDs in BED, written to: {diff_file}"))
                
                # Check Protein duplicates
                if stats["protein_duplicates"]:
                    if len(stats["protein_duplicates"]) <= MAX_DIFF_DISPLAY:
                        for gene_id, count in stats["protein_duplicates"]:
                            print(Colors.yellow(f"  Protein duplicate: {gene_id} (x{count})"))
                    else:
                        diff_file = write_duplicates_to_file(stats["protein_duplicates"], "Protein", args.output, abbr)
                        print(Colors.yellow(f"  >20 duplicate IDs in Protein, written to: {diff_file}"))
                
                # Check CDS duplicates
                if stats["cds_duplicates"]:
                    if len(stats["cds_duplicates"]) <= MAX_DIFF_DISPLAY:
                        for gene_id, count in stats["cds_duplicates"]:
                            print(Colors.yellow(f"  CDS duplicate: {gene_id} (x{count})"))
                    else:
                        diff_file = write_duplicates_to_file(stats["cds_duplicates"], "CDS", args.output, abbr)
                        print(Colors.yellow(f"  >20 duplicate IDs in CDS, written to: {diff_file}"))
            
            # Check chromosome IDs
            if stats["chr_ids_match"]:
                print(Colors.green(f"✓ Chromosome IDs match: BED has {stats['bed_chr_unique']}, Genome has {stats['genome_chr_unique']}"))
            else:
                print(Colors.red("✗ Chromosome mismatch"))
                if stats["bed_chr_only"]:
                    if len(stats["bed_chr_only"]) <= MAX_DIFF_DISPLAY:
                        for chr_id in sorted(stats["bed_chr_only"]):
                            print(Colors.yellow(f"  BED only: {chr_id}"))
                    else:
                        diff_file = write_diff_to_file(stats["bed_chr_only"], "BED_chr", "Genome", args.output, abbr)
                        print(Colors.yellow(f"  >20 unique chromosomes in BED, written to: {diff_file}"))
                
                if stats["genome_only"]:
                    if len(stats["genome_only"]) <= MAX_DIFF_DISPLAY:
                        for chr_id in sorted(stats["genome_only"]):
                            print(Colors.yellow(f"  Genome only: {chr_id}"))
                    else:
                        diff_file = write_diff_to_file(stats["genome_only"], "Genome", "BED_chr", args.output, abbr)
                        print(Colors.yellow(f"  >20 unique chromosomes in Genome, written to: {diff_file}"))
            
            # Final verdict
            if is_valid:
                print(Colors.green(f"\n✓ All checks passed for '{abbr}'!"))
            else:
                print(Colors.red(f"\n✗ Validation failed for '{abbr}'"))
                exit_code = EXIT_VALIDATION_FAILED
                
        except Exception as e:
            print(Colors.red(f"✗ Error: {e}"))
            exit_code = EXIT_VALIDATION_FAILED

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
