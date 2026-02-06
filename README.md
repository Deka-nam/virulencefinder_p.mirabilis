# virulencefinder_p.mirabilis
Ready-to-use package for detecting Proteus mirabilis virulence factors from genome assemblies. 

# Proteus mirabilis Virulence Factor Detection Tool

## Overview

This tool provides a standalone solution for detecting **Proteus mirabilis** virulence factors in genome assemblies. Unlike the generic VFDB database used by ABRicate, this curated database includes P. mirabilis-specific sequences from both VFDB and NCBI reference genomes.

## Why This Tool?

When running ABRicate with the standard VFDB database, P. mirabilis-specific virulence factors are often missed because:
1. VFDB focuses on well-characterized pathogens (E. coli, Salmonella, etc.)
2. P. mirabilis virulence genes have lower sequence similarity to related species
3. Important uropathogenic factors (urease, MR/P fimbriae) need species-specific sequences

This curated database includes **61 virulence genes** specifically from P. mirabilis, ensuring comprehensive detection.

## Files

- `pmirabilis_vf_detect.py` - Main detection script
- `pmirabilis_vfdb.fasta` - Curated virulence factor database (61 sequences)

## Requirements

- Python 3.6+
- NCBI BLAST+ (blastn, makeblastdb)

### Installing BLAST+

```bash
# Ubuntu/Debian
sudo apt install ncbi-blast+

# macOS
brew install blast

# Conda
conda install -c bioconda blast
```

## Usage

### Basic Detection

```bash
python pmirabilis_vf_detect.py assembly.fasta
```

### Output Formats

```bash
# Table format (default, human-readable)
python pmirabilis_vf_detect.py assembly.fasta

# TSV format (machine-readable, similar to ABRicate)
python pmirabilis_vf_detect.py assembly.fasta -f tsv

# GFF format (for genome browsers)
python pmirabilis_vf_detect.py assembly.fasta -f gff
```

### Batch Analysis

```bash
# Process multiple assemblies
python pmirabilis_vf_detect.py *.fasta -f tsv > all_results.tsv

# Or use a loop for individual files
for f in assemblies/*.fasta; do
    python pmirabilis_vf_detect.py "$f" -f tsv >> batch_results.tsv
done
```

### Adjusting Parameters

```bash
# Higher stringency (90% identity, 80% coverage)
python pmirabilis_vf_detect.py assembly.fasta --min-id 90 --min-cov 80

# Lower stringency for divergent strains
python pmirabilis_vf_detect.py assembly.fasta --min-id 70 --min-cov 60

# Save to file
python pmirabilis_vf_detect.py assembly.fasta -o results.tsv -f tsv
```

## Output Interpretation

### Table Output Columns

| Column | Description |
|--------|-------------|
| Gene | Gene name/symbol |
| Contig | Contig/chromosome name |
| Start | Start position in genome |
| End | End position in genome |
| Strand | + (forward) or - (reverse) |
| %ID | Percent sequence identity |
| %Cov | Percent of gene covered |
| Category | Functional category |

### TSV Output Columns

| Column | Description |
|--------|-------------|
| #FILE | Input filename |
| SEQUENCE | Contig name |
| START | Start coordinate |
| END | End coordinate |
| STRAND | Orientation |
| GENE | Gene name |
| COVERAGE | Percent coverage |
| %IDENTITY | Percent identity |
| DATABASE | pmirabilis_vfdb |
| ACCESSION | Database accession |
| PRODUCT | Gene product description |
| CATEGORY | Functional category |

## Virulence Factor Categories

| Category | Description | Key Genes |
|----------|-------------|-----------|
| **Adherence** | Fimbriae and adhesins for attachment | mrpA-J, ucaA, PMF genes |
| **Urease** | Struvite stone formation, urinary pH | ureA-G (complete operon) |
| **Hemolysin** | Tissue damage and lysis | hpmA, hpmB |
| **Autotransporter** | Secreted toxins | pta, aipA, taaP |
| **Flagella/Motility** | Swarming motility | flhD, flhC |
| **Iron Acquisition** | Siderophore systems | HmuR2 system genes |
| **Protease** | IgA degradation | zapA |

## Database Contents (61 genes)

### Adherence (20 genes)
- MR/P fimbriae: mrpA, mrpB, mrpC, mrpD, mrpH, mrpI, mrpJ
- UCA fimbriae: ucaA
- PMF fimbriae: PMI_RS09265-09285
- Other fimbriae: PMI_RS01295-01355, PMI_RS02630-02645

### Urease (6 genes)
- Complete operon: ureA, ureB, ureC, ureE, ureF, ureG

### Hemolysin (2 genes)
- HpmA-HpmB system: hpmA, hpmB

### Autotransporter (3 genes)
- pta (Proteus toxic agglutinin)
- aipA (Autotransporter iron-regulated)
- taaP (Trimeric autotransporter adhesin)

### Flagella/Motility (2 genes)
- Master regulators: flhD, flhC

### Iron Acquisition (17 genes)
- HmuR2 system: PMI_RS06900-06920
- Yersiniabactin-like: PMI_RS12820-12865

### Protease (1 gene)
- IgA protease: zapA

## Comparison with ABRicate

| Feature | ABRicate + VFDB | This Tool |
|---------|-----------------|-----------|
| P. mirabilis genes | ~10-15 | 61 |
| Urease detection | Poor | Complete operon |
| MR/P fimbriae | Partial | All 7 genes |
| Coordinates | Yes | Yes |
| Categories | Generic | P. mirabilis specific |

## Clinical Significance

For uropathogenic P. mirabilis, key virulence factors to look for:

1. **Urease (ureA-G)**: Essential for struvite stone formation and ascending infection
2. **MR/P fimbriae (mrp genes)**: Key adhesins for urinary tract colonization
3. **Hemolysin (hpmAB)**: Causes tissue damage and hemolysis
4. **Flagella (flhDC)**: Required for swarming motility and biofilm formation
5. **ZapA protease**: IgA degradation for immune evasion

## Citation

If you use this database, please cite:
- Chen, L., Yang, J., Yu, J., Yao, Z., Sun, L., Shen, Y., & Jin, Q. (2005). VFDB: a reference database for bacterial virulence factors. Nucleic acids research, 33(Database issue), D325–D328. https://doi.org/10.1093/nar/gki008
- Pearson, M. M., Sebaihia, M., Churcher, C., Quail, M. A., Seshasayee, A. S., Luscombe, N. M., Abdellah, Z., Arrosmith, C., Atkin, B., Chillingworth, T., Hauser, H., Jagels, K., Moule, S., Mungall, K., Norbertczak, H., Rabbinowitsch, E., Walker, D., Whithead, S., Thomson, N. R., Rather, P. N., … Mobley, H. L. (2008). Complete genome sequence of uropathogenic Proteus mirabilis, a master of both adherence and motility. Journal of bacteriology, 190(11), 4027–4037. https://doi.org/10.1128/JB.01981-07
- Deka, N., Brauer, A. L., Connerton, K., Hanson, B., Walker, J. N., & Armbruster, C. E. (2025). Pangenome Analysis of Proteus mirabilis Reveals Lineage-Specific Antimicrobial Resistance Profiles and Discordant Genotype-Phenotype Correlations. bioRxiv : the preprint server for biology, 2025.11.21.689858. https://doi.org/10.1101/2025.11.21.689858

## Version

- Database version: 1.0
- Last updated: 2025
- Based on: VFDB setB + NCBI AM942759 (P. mirabilis HI4320)
