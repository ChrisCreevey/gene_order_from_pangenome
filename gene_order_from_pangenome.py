#!/usr/bin/env python3
"""
Generate ordered gene-family lists from a Roary presence-absence CSV and GFF annotation files.

For each genome defined in the presence-absence file, outputs a comma-separated list of
gene-family names in the order genes appear in the corresponding GFF file. Paralogs (multiple
copies of a gene family in one genome) are each listed at the position they occupy in the GFF.

Usage:
    python gene_order_from_pangenome.py \\
        --presence-absence gene_presence_absence.csv \\
        --gff-dir /path/to/gff/files \\
        [--output results.txt] \\
        [--unknown skip|TOKEN]
"""

import csv
import argparse
import sys
from pathlib import Path


def parse_attributes(attr_string):
    """Parse a GFF3 attribute string into a dictionary."""
    attrs = {}
    for part in attr_string.strip().split(';'):
        part = part.strip()
        if '=' in part:
            key, _, value = part.partition('=')
            attrs[key.strip()] = value.strip()
    return attrs


def parse_presence_absence(pa_path, n_meta_cols=14):
    """
    Parse a Roary-style presence-absence CSV file.

    Returns:
        genome_names (list[str]): Genome identifiers taken from the header row.
        gene_lookup (dict[str, str]): Maps each locus_tag to its gene-family name.
    """
    gene_lookup = {}
    genome_names = []

    with open(pa_path, newline='', encoding='utf-8') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        genome_names = [col.strip() for col in header[n_meta_cols:]]

        for row in reader:
            if not row:
                continue
            gene_family = row[0].strip()
            for cell in row[n_meta_cols:]:
                cell = cell.strip()
                if not cell:
                    continue
                # Roary separates paralogous copies within a cell with tab characters
                for locus_tag in cell.split('\t'):
                    locus_tag = locus_tag.strip()
                    if locus_tag:
                        gene_lookup[locus_tag] = gene_family

    return genome_names, gene_lookup


def parse_gff_order(gff_path):
    """
    Parse a GFF3 file and return locus_tags of CDS features in file order.

    Uses the 'locus_tag' attribute as the primary identifier, falling back to 'ID'.
    Skips duplicate locus_tags to handle genes split across multiple GFF rows.
    Stops reading if a FASTA section (lines starting with '>') is encountered.
    """
    seen = set()
    order = []

    with open(gff_path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            if line.startswith('>'):  # embedded FASTA section
                break
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue
            attrs = parse_attributes(parts[8])
            locus_tag = attrs.get('locus_tag') or attrs.get('ID', '').split('.')[0]
            if locus_tag and locus_tag not in seen:
                seen.add(locus_tag)
                order.append(locus_tag)

    return order


def find_gff(gff_dir, genome_name):
    """Return the GFF file path for a genome, or None if not found."""
    for ext in ('.gff', '.gff3'):
        candidate = gff_dir / f'{genome_name}{ext}'
        if candidate.is_file():
            return candidate
    return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Generate ordered gene-family lists from a Roary presence-absence file '
            'and per-genome GFF annotation files.'
        )
    )
    parser.add_argument(
        '--presence-absence', '-p',
        required=True,
        metavar='FILE',
        help='Roary gene_presence_absence.csv file'
    )
    parser.add_argument(
        '--gff-dir', '-g',
        required=True,
        metavar='DIR',
        help='Directory containing GFF files, named <genome_name>.gff or <genome_name>.gff3'
    )
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        help='Output file (default: stdout)'
    )
    parser.add_argument(
        '--unknown', '-u',
        default='*',
        metavar='TOKEN',
        help=(
            'Token for CDS features not found in the presence-absence file '
            '(default: *). Use "skip" to omit them entirely.'
        )
    )
    parser.add_argument(
        '--meta-cols', '-m',
        type=int,
        default=14,
        help='Number of metadata columns before genome columns (default: 14)'
    )
    args = parser.parse_args()

    pa_path = Path(args.presence_absence)
    gff_dir = Path(args.gff_dir)

    if not pa_path.is_file():
        sys.exit(f'Error: presence-absence file not found: {pa_path}')
    if not gff_dir.is_dir():
        sys.exit(f'Error: GFF directory not found: {gff_dir}')

    print(f'Parsing presence-absence file: {pa_path}', file=sys.stderr)
    genome_names, gene_lookup = parse_presence_absence(pa_path, args.meta_cols)
    print(
        f'  {len(genome_names)} genomes, {len(gene_lookup)} locus_tag mappings',
        file=sys.stderr
    )

    out_fh = open(args.output, 'w', encoding='utf-8') if args.output else sys.stdout
    try:
        for genome in genome_names:
            gff_path = find_gff(gff_dir, genome)
            if gff_path is None:
                print(f'Warning: no GFF file found for genome "{genome}"', file=sys.stderr)
                continue

            print(f'Processing {genome} ({gff_path.name})', file=sys.stderr)
            locus_tag_order = parse_gff_order(gff_path)

            result = []
            for lt in locus_tag_order:
                family = gene_lookup.get(lt)
                if family:
                    result.append(family)
                elif args.unknown != 'skip':
                    result.append(args.unknown)

            out_fh.write(f'>{genome}\n')
            out_fh.write(','.join(result) + '\n')

    finally:
        if args.output:
            out_fh.close()

    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()
