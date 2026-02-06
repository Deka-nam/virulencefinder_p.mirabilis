#!/usr/bin/env python3
"""
Proteus mirabilis Virulence Factor Detection Tool
=================================================

A standalone script for detecting virulence factors specific to Proteus mirabilis
in genome assemblies. Uses a curated database combining VFDB and NCBI reference genes.

Usage:
    python pmirabilis_vf_detect.py assembly.fasta [options]

Requirements:
    - NCBI BLAST+ (makeblastdb, blastn)
    - Python 3.6+ with BioPython (optional, for parsing)

Author: Generated for P. mirabilis uropathogen analysis
Database: VFDB setB + NCBI P. mirabilis HI4320 reference genes (61 sequences)
"""

import subprocess
import sys
import os
import argparse
import tempfile
import csv
import re
from pathlib import Path
from collections import defaultdict

# Database path - set to the directory containing this script
SCRIPT_DIR = Path(__file__).parent.absolute()
DEFAULT_DB = SCRIPT_DIR / "pmirabilis_vfdb.fasta"

# Virulence factor categories for P. mirabilis
VF_CATEGORIES = {
    'Adherence': ['PMF', 'MR/P', 'UCA', 'fimbrial', 'adhesin', 'pili'],
    'Hemolysin': ['hpmA', 'hpmB', 'hemolysin'],
    'Urease': ['ureA', 'ureB', 'ureC', 'ureD', 'ureE', 'ureF', 'ureG', 'urease'],
    'Autotransporter': ['pta', 'taaP', 'aipA', 'autotransporter', 'AipA', 'TaaP', 'Pta'],
    'Flagella/Motility': ['flhD', 'flhC', 'flagell', 'flg', 'fli'],
    'Iron Acquisition': ['hmu', 'irp', 'yersiniabactin', 'siderophore', 'iron'],
    'Protease': ['zapA', 'protease', 'ZapA'],
    'Other': []
}

def check_blast_available():
    """Check if BLAST+ is installed and available."""
    try:
        result = subprocess.run(['blastn', '-version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def create_blast_db(fasta_file, db_name=None):
    """Create a BLAST database from the assembly FASTA file."""
    if db_name is None:
        db_name = fasta_file
    cmd = ['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl', '-out', db_name]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error creating BLAST database: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return db_name

def run_blast(query_db, subject_assembly, min_identity=80, min_coverage=70, evalue=1e-10):
    """Run BLAST search of virulence genes against the assembly."""

    # BLAST output format 6 with custom columns
    outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'

    cmd = [
        'blastn',
        '-query', str(query_db),
        '-subject', str(subject_assembly),
        '-outfmt', outfmt,
        '-evalue', str(evalue),
        '-perc_identity', str(min_identity),
        '-max_target_seqs', '10',
        '-word_size', '11',
        '-dust', 'no'  # Don't mask low-complexity regions
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"BLAST error: {result.stderr}", file=sys.stderr)
        return []

    return result.stdout.strip().split('\n') if result.stdout.strip() else []

def parse_vf_header(header):
    """Parse virulence factor header to extract gene name, product, VF class."""
    header = header.lstrip('>')

    # Pattern for VFG style headers
    vfg_pattern = r'^(VFG\d+)\(([^)]+)\)\s+\(([^)]+)\)\s+(.+?)\s+\[([^\]]+)\]\s+\[([^\]]+)\]$'

    # Pattern for NCBI style headers with coordinates
    ncbi_pattern = r'^([^:]+):(\d+)\.\.(\d+)\s+(\S+)\s+\[([^\]]+)\]$'

    match = re.match(vfg_pattern, header)
    if match:
        vfg_id, accession, gene_name, product, vf_class, organism = match.groups()
        return {
            'id': vfg_id,
            'accession': accession.replace('gb|', '').replace('|', '_'),
            'gene': gene_name,
            'product': product,
            'vf_class': vf_class,
            'organism': organism
        }

    match = re.match(ncbi_pattern, header)
    if match:
        ref_id, start, end, gene_name, organism = match.groups()
        return {
            'id': f"{ref_id}_{start}_{end}",
            'accession': ref_id,
            'gene': gene_name,
            'product': gene_name,
            'vf_class': 'Reference gene',
            'organism': organism
        }

    # Fallback parsing
    parts = header.split()
    return {
        'id': parts[0] if parts else 'unknown',
        'accession': 'unknown',
        'gene': parts[1] if len(parts) > 1 else 'unknown',
        'product': ' '.join(parts[1:]) if len(parts) > 1 else header,
        'vf_class': 'unknown',
        'organism': 'Proteus mirabilis'
    }

def load_database_info(db_fasta):
    """Load information about all virulence factors in the database."""
    vf_info = {}
    current_header = None

    with open(db_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                current_header = line
                info = parse_vf_header(line)
                # Create a unique key that matches the BLAST qseqid
                # BLAST uses the first word of the header as qseqid
                qseqid = line[1:].split()[0]
                vf_info[qseqid] = info

    return vf_info

def categorize_vf(gene, product, vf_class):
    """Categorize a virulence factor based on its annotation."""
    combined = f"{gene} {product} {vf_class}".lower()

    for category, keywords in VF_CATEGORIES.items():
        if category == 'Other':
            continue
        for keyword in keywords:
            if keyword.lower() in combined:
                return category
    return 'Other'

def format_coordinates(sstart, send, strand):
    """Format genomic coordinates nicely."""
    if strand == '+':
        return f"{sstart}..{send} (+)"
    else:
        return f"{send}..{sstart} (-)"

def parse_blast_results(blast_lines, vf_info, min_coverage=70):
    """Parse BLAST results and format output."""
    results = []

    for line in blast_lines:
        if not line.strip():
            continue

        fields = line.split('\t')
        if len(fields) < 15:
            continue

        qseqid = fields[0]
        sseqid = fields[1]
        pident = float(fields[2])
        length = int(fields[3])
        qstart = int(fields[6])
        qend = int(fields[7])
        sstart = int(fields[8])
        send = int(fields[9])
        evalue = float(fields[10])
        bitscore = float(fields[11])
        qlen = int(fields[12])
        slen = int(fields[13])
        qcovs = float(fields[14])

        # Calculate coverage
        coverage = (length / qlen) * 100

        if coverage < min_coverage:
            continue

        # Determine strand
        strand = '+' if sstart < send else '-'

        # Get gene info from database
        info = vf_info.get(qseqid, {
            'gene': qseqid.split('(')[1].split(')')[0] if '(' in qseqid else qseqid,
            'product': qseqid,
            'vf_class': 'unknown'
        })

        category = categorize_vf(info.get('gene', ''), 
                                  info.get('product', ''), 
                                  info.get('vf_class', ''))

        result = {
            'gene': info.get('gene', 'unknown'),
            'product': info.get('product', 'unknown'),
            'vf_class': info.get('vf_class', 'unknown'),
            'category': category,
            'contig': sseqid,
            'start': min(sstart, send),
            'end': max(sstart, send),
            'strand': strand,
            'identity': pident,
            'coverage': coverage,
            'evalue': evalue,
            'bitscore': bitscore,
            'query_id': qseqid
        }
        results.append(result)

    return results

def print_results(results, output_format='table'):
    """Print results in specified format."""
    if not results:
        print("\nNo P. mirabilis virulence factors detected.", file=sys.stderr)
        return

    # Remove duplicates (keep best hit per gene per contig)
    unique_results = {}
    for r in results:
        key = (r['gene'], r['contig'], r['start'])
        if key not in unique_results or r['bitscore'] > unique_results[key]['bitscore']:
            unique_results[key] = r

    results = sorted(unique_results.values(), key=lambda x: (x['contig'], x['start']))

    if output_format == 'tsv':
        # TSV format compatible with ABRicate output
        print("\t".join(['#FILE', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE', 
                         'COVERAGE', '%IDENTITY', 'DATABASE', 'ACCESSION', 
                         'PRODUCT', 'CATEGORY']))
        for r in results:
            print("\t".join([
                'input',
                r['contig'],
                str(r['start']),
                str(r['end']),
                r['strand'],
                r['gene'],
                f"{r['coverage']:.1f}",
                f"{r['identity']:.1f}",
                'pmirabilis_vfdb',
                r['query_id'].split('(')[0] if '(' in r['query_id'] else r['query_id'],
                r['product'],
                r['category']
            ]))

    elif output_format == 'table':
        # Human-readable table format
        print("\n" + "="*100)
        print("PROTEUS MIRABILIS VIRULENCE FACTOR DETECTION RESULTS")
        print("="*100)
        print(f"\nTotal virulence genes detected: {len(results)}")

        # Group by category
        by_category = defaultdict(list)
        for r in results:
            by_category[r['category']].append(r)

        print(f"Categories found: {len(by_category)}")
        for cat in sorted(by_category.keys()):
            print(f"  - {cat}: {len(by_category[cat])} genes")

        print("\n" + "-"*100)
        print(f"{'Gene':<12} {'Contig':<20} {'Start':>10} {'End':>10} {'Strand':<6} {'%ID':>6} {'%Cov':>6} {'Category':<15}")
        print("-"*100)

        for r in results:
            contig = r['contig'][:18] + '..' if len(r['contig']) > 20 else r['contig']
            print(f"{r['gene']:<12} {contig:<20} {r['start']:>10} {r['end']:>10} {r['strand']:<6} {r['identity']:>5.1f}% {r['coverage']:>5.1f}% {r['category']:<15}")

        print("-"*100)
        print()

    elif output_format == 'gff':
        # GFF3 format for genome visualization
        print("##gff-version 3")
        for r in results:
            attributes = f"ID={r['gene']};gene={r['gene']};product={r['product']};category={r['category']};identity={r['identity']:.1f};coverage={r['coverage']:.1f}"
            print(f"{r['contig']}\tpmirabilis_vfdb\tvirulence_factor\t{r['start']}\t{r['end']}\t{r['bitscore']:.1f}\t{r['strand']}\t.\t{attributes}")

def main():
    parser = argparse.ArgumentParser(
        description='Detect P. mirabilis virulence factors in genome assemblies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s assembly.fasta                    # Basic detection with table output
  %(prog)s assembly.fasta -o results.tsv     # Save results to TSV file
  %(prog)s assembly.fasta -f gff             # Output in GFF3 format
  %(prog)s assembly.fasta --min-id 90        # Require 90%% identity
  %(prog)s *.fasta -f tsv > batch_results.tsv  # Batch analysis

Database:
  Uses curated P. mirabilis virulence factors from VFDB and NCBI (61 genes)
  Categories: Adherence, Hemolysin, Urease, Autotransporter, Flagella/Motility,
              Iron Acquisition, Protease, Other
        """
    )

    parser.add_argument('assembly', nargs='+', help='Assembly FASTA file(s)')
    parser.add_argument('-d', '--database', default=None, 
                        help='Path to virulence factor database FASTA (default: pmirabilis_vfdb.fasta in script directory)')
    parser.add_argument('-o', '--output', default=None, help='Output file (default: stdout)')
    parser.add_argument('-f', '--format', choices=['table', 'tsv', 'gff'], default='table',
                        help='Output format (default: table)')
    parser.add_argument('--min-id', type=float, default=80, 
                        help='Minimum percent identity (default: 80)')
    parser.add_argument('--min-cov', type=float, default=70,
                        help='Minimum percent coverage (default: 70)')
    parser.add_argument('--evalue', type=float, default=1e-10,
                        help='E-value threshold (default: 1e-10)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')

    args = parser.parse_args()

    # Check BLAST availability
    if not check_blast_available():
        print("Error: BLAST+ (blastn) not found. Please install NCBI BLAST+.", file=sys.stderr)
        print("  Ubuntu/Debian: sudo apt install ncbi-blast+", file=sys.stderr)
        print("  macOS: brew install blast", file=sys.stderr)
        print("  Conda: conda install -c bioconda blast", file=sys.stderr)
        sys.exit(1)

    # Determine database path
    if args.database:
        db_path = Path(args.database)
    else:
        db_path = DEFAULT_DB

    if not db_path.exists():
        print(f"Error: Database file not found: {db_path}", file=sys.stderr)
        print("Please ensure pmirabilis_vfdb.fasta is in the same directory as this script", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Using database: {db_path}", file=sys.stderr)

    # Load database information
    vf_info = load_database_info(db_path)
    if args.verbose:
        print(f"Loaded {len(vf_info)} virulence factor sequences", file=sys.stderr)

    # Redirect output if specified
    output_file = open(args.output, 'w') if args.output else sys.stdout
    original_stdout = sys.stdout
    sys.stdout = output_file

    try:
        # Process each assembly
        all_results = []
        for assembly in args.assembly:
            if not Path(assembly).exists():
                print(f"Warning: Assembly file not found: {assembly}", file=sys.stderr)
                continue

            if args.verbose:
                print(f"Processing: {assembly}", file=sys.stderr)

            # Run BLAST
            blast_results = run_blast(db_path, assembly, args.min_id, args.min_cov, args.evalue)

            if args.verbose:
                print(f"  BLAST hits: {len(blast_results)}", file=sys.stderr)

            # Parse results
            results = parse_blast_results(blast_results, vf_info, args.min_cov)

            # Add filename to results
            for r in results:
                r['file'] = assembly

            all_results.extend(results)

        # Print results
        print_results(all_results, args.format)

    finally:
        if args.output:
            output_file.close()
        sys.stdout = original_stdout

if __name__ == '__main__':
    main()
