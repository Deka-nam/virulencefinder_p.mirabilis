#!/bin/bash
#
# pmirabilis_vf.sh - Wrapper script for P. mirabilis virulence factor detection
#
# Usage: 
#   ./pmirabilis_vf.sh assembly.fasta
#   ./pmirabilis_vf.sh assembly.fasta tsv output.tsv
#   ./pmirabilis_vf.sh assembly.fasta gff output.gff
#

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/pmirabilis_vf_detect.py"
DATABASE="$SCRIPT_DIR/pmirabilis_vfdb.fasta"

# Check if python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: pmirabilis_vf_detect.py not found in $SCRIPT_DIR"
    exit 1
fi

# Check if database exists
if [ ! -f "$DATABASE" ]; then
    echo "Error: pmirabilis_vfdb.fasta not found in $SCRIPT_DIR"
    exit 1
fi

# Check arguments
if [ $# -lt 1 ]; then
    echo "P. mirabilis Virulence Factor Detection"
    echo ""
    echo "Usage:"
    echo "  $0 <assembly.fasta>                  # Table output to stdout"
    echo "  $0 <assembly.fasta> tsv              # TSV output to stdout"
    echo "  $0 <assembly.fasta> tsv output.tsv   # TSV output to file"
    echo "  $0 <assembly.fasta> gff output.gff   # GFF output to file"
    echo ""
    echo "Options:"
    echo "  table  - Human-readable table format"
    echo "  tsv    - Tab-separated values (ABRicate-like)"
    echo "  gff    - GFF3 format for genome browsers"
    exit 0
fi

ASSEMBLY="$1"
FORMAT="${2:-table}"
OUTPUT="${3:-}"

if [ ! -f "$ASSEMBLY" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY"
    exit 1
fi

# Run detection
if [ -n "$OUTPUT" ]; then
    python3 "$PYTHON_SCRIPT" "$ASSEMBLY" -f "$FORMAT" -o "$OUTPUT"
    echo "Results written to: $OUTPUT"
else
    python3 "$PYTHON_SCRIPT" "$ASSEMBLY" -f "$FORMAT"
fi
