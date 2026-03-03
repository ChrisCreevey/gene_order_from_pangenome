#!/usr/bin/env python3
"""
Generate ordered gene-family lists from a Roary presence-absence CSV and GFF annotation files.

For each genome defined in the presence-absence file, outputs a comma-separated list of
gene-family names in the order genes appear in the corresponding GFF file. Paralogs (multiple
copies of a gene family in one genome) are each listed at the position they occupy in the GFF.

Optional features:
  --contig-sep SEP : insert a separator token wherever the contig/scaffold changes
  --strand mark    : append + or - to each gene-family name to show strand orientation
  --strand split   : write separate .plus and .minus output files instead of a combined file
  --reverse-minus  : reverse the gene order within each contig in the minus strand output,
                     so genes read 5'→3' along the minus strand (only meaningful with
                     --strand split or --strand mark)

Usage:
    python gene_order_from_pangenome.py \\
        --presence-absence gene_presence_absence.csv \\
        --gff-dir /path/to/gff/files \\
        [--output results.txt] \\
        [--unknown skip|TOKEN] \\
        [--contig-sep |] \\
        [--strand mark|split] \\
        [--reverse-minus]
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
    Parse a GFF3 file and return a list of (locus_tag, seqid, strand) tuples
    for CDS features in file order.

    Uses the 'locus_tag' attribute as the primary identifier, falling back to 'ID'.
    Skips duplicate locus_tags to handle genes split across multiple GFF rows.
    Stops reading if a FASTA section (lines starting with '>') is encountered.
    """
    seen = set()
    entries = []

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
            seqid = parts[0]
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            locus_tag = attrs.get('locus_tag') or attrs.get('ID', '').split('.')[0]
            if locus_tag and locus_tag not in seen:
                seen.add(locus_tag)
                entries.append((locus_tag, seqid, strand))

    return entries


def find_gff(gff_dir, genome_name):
    """Return the GFF file path for a genome, or None if not found."""
    for ext in ('.gff', '.gff3'):
        candidate = gff_dir / f'{genome_name}{ext}'
        if candidate.is_file():
            return candidate
    return None


def process_genome(entries, gene_lookup, unknown_token, mark_strand, contig_sep):
    """
    Convert (locus_tag, seqid, strand) entries into gene-family name lists.

    Returns:
        all_genes   (list[str]): All genes in GFF order, with optional strand marks
                                  and contig separators.
        plus_genes  (list[str]): Genes on the + strand only, in GFF order.
        minus_genes (list[str]): Genes on the - strand only, in GFF order.

    Contig separators (if requested) are inserted into all three lists whenever the
    sequence ID changes, regardless of strand.
    """
    all_genes = []
    plus_genes = []
    minus_genes = []
    prev_seqid = None

    for lt, seqid, strand in entries:
        # Insert contig boundary token when the sequence ID changes
        if contig_sep is not None and prev_seqid is not None and seqid != prev_seqid:
            all_genes.append(contig_sep)
            plus_genes.append(contig_sep)
            minus_genes.append(contig_sep)
        prev_seqid = seqid

        family = gene_lookup.get(lt)
        if family:
            label = f'{family}{strand}' if mark_strand else family
        elif unknown_token != 'skip':
            family = unknown_token
            label = f'{unknown_token}{strand}' if mark_strand else unknown_token
        else:
            continue  # omit unknown genes entirely

        all_genes.append(label)
        if strand == '+':
            plus_genes.append(label)
        else:
            minus_genes.append(label)

    return all_genes, plus_genes, minus_genes


def reverse_within_contigs(genes, contig_sep):
    """
    Reverse gene order within each contig segment independently.

    GFF files list all genes in reference (+ strand) coordinate order. For minus strand
    genes this means they appear 3'→5'. Reversing within each contig corrects this so
    that genes read 5'→3' along the minus strand, while keeping contigs in their
    original order.

    If no contig_sep is set the entire list is simply reversed.
    """
    if not contig_sep:
        return list(reversed(genes))

    # Split into contig segments, reverse each, then rejoin
    segments = []
    current = []
    for gene in genes:
        if gene == contig_sep:
            segments.append(current)
            current = []
        else:
            current.append(gene)
    segments.append(current)

    result = []
    for i, segment in enumerate(segments):
        result.extend(reversed(segment))
        if i < len(segments) - 1:
            result.append(contig_sep)
    return result


def strand_output_paths(base_output):
    """Derive + and - strand output file paths from the base output path."""
    if base_output is None:
        return Path('out.plus.txt'), Path('out.minus.txt')
    p = Path(base_output)
    suffix = p.suffix or '.txt'
    return (
        p.parent / f'{p.stem}.plus{suffix}',
        p.parent / f'{p.stem}.minus{suffix}',
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Generate ordered gene-family lists from a Roary presence-absence file '
            'and per-genome GFF annotation files.'
        )
    )
    parser.add_argument(
        '--presence-absence', '-p', required=True, metavar='FILE',
        help='Roary gene_presence_absence.csv file'
    )
    parser.add_argument(
        '--gff-dir', '-g', required=True, metavar='DIR',
        help='Directory containing GFF files, named <genome_name>.gff or <genome_name>.gff3'
    )
    parser.add_argument(
        '--output', '-o', metavar='FILE',
        help=(
            'Output file (default: stdout). '
            'With --strand split, used as the base name for .plus and .minus files '
            '(e.g. results.txt → results.plus.txt and results.minus.txt).'
        )
    )
    parser.add_argument(
        '--unknown', '-u', default='*', metavar='TOKEN',
        help=(
            'Token for CDS features not found in the presence-absence file '
            '(default: *). Use "skip" to omit them entirely.'
        )
    )
    parser.add_argument(
        '--meta-cols', '-m', type=int, default=14,
        help='Number of metadata columns before genome columns (default: 14)'
    )
    parser.add_argument(
        '--contig-sep', metavar='SEP',
        help=(
            'Token inserted into the output list at every contig/scaffold boundary '
            'to indicate that adjacent genes are on different genome fragments '
            '(e.g. "--contig-sep |"). Off by default.'
        )
    )
    parser.add_argument(
        '--strand', choices=['mark', 'split'],
        help=(
            'Strand handling. '
            '"mark": append + or - to each gene-family name (e.g. ppnP+, lptC-). '
            '"split": write separate output files for the + and - strands '
            'instead of a single combined file.'
        )
    )
    parser.add_argument(
        '--reverse-minus', action='store_true',
        help=(
            'Reverse the order of minus strand genes within each contig so they read '
            '5\'→3\' along the minus strand. In a GFF file all genes are listed in '
            'reference coordinate order, which means minus strand genes appear 3\'→5\'; '
            'this flag corrects that. Most useful with --strand split or --strand mark.'
        )
    )
    args = parser.parse_args()

    pa_path = Path(args.presence_absence)
    gff_dir = Path(args.gff_dir)

    if not pa_path.is_file():
        sys.exit(f'Error: presence-absence file not found: {pa_path}')
    if not gff_dir.is_dir():
        sys.exit(f'Error: GFF directory not found: {gff_dir}')

    mark_strand = args.strand == 'mark'
    split_strand = args.strand == 'split'

    print(f'Parsing presence-absence file: {pa_path}', file=sys.stderr)
    genome_names, gene_lookup = parse_presence_absence(pa_path, args.meta_cols)
    print(
        f'  {len(genome_names)} genomes, {len(gene_lookup)} locus_tag mappings',
        file=sys.stderr
    )

    # Open output file handles
    if split_strand:
        plus_path, minus_path = strand_output_paths(args.output)
        plus_fh = open(plus_path, 'w', encoding='utf-8')
        minus_fh = open(minus_path, 'w', encoding='utf-8')
        main_fh = None
        print(f'Strand split outputs: {plus_path}, {minus_path}', file=sys.stderr)
    else:
        main_fh = open(args.output, 'w', encoding='utf-8') if args.output else sys.stdout
        plus_fh = minus_fh = None

    try:
        for genome in genome_names:
            gff_path = find_gff(gff_dir, genome)
            if gff_path is None:
                print(f'Warning: no GFF file found for genome "{genome}"', file=sys.stderr)
                continue

            print(f'Processing {genome} ({gff_path.name})', file=sys.stderr)
            entries = parse_gff_order(gff_path)
            all_genes, plus_genes, minus_genes = process_genome(
                entries, gene_lookup, args.unknown, mark_strand, args.contig_sep
            )

            if args.reverse_minus:
                minus_genes = reverse_within_contigs(minus_genes, args.contig_sep)

            if split_strand:
                plus_fh.write(f'>{genome}\n')
                plus_fh.write(','.join(plus_genes) + '\n')
                minus_fh.write(f'>{genome}\n')
                minus_fh.write(','.join(minus_genes) + '\n')
            else:
                main_fh.write(f'>{genome}\n')
                main_fh.write(','.join(all_genes) + '\n')

    finally:
        if main_fh and args.output:
            main_fh.close()
        if plus_fh:
            plus_fh.close()
        if minus_fh:
            minus_fh.close()

    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()
