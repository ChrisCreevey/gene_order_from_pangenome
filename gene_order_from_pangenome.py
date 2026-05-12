#!/usr/bin/env python3
"""
Generate ordered gene-family lists from a Roary presence-absence CSV and GFF annotation files.

For each genome defined in the presence-absence file, outputs a comma-separated list of
gene-family names in the order genes appear in the corresponding GFF file. Paralogs (multiple
copies of a gene family in one genome) are each listed at the position they occupy in the GFF.

Optional features:
  --contig-sep SEP : write each contig/scaffold as a separate output sequence
  --strand mark    : append + or - to each gene-family name to show strand orientation
  --strand split   : write separate .plus and .minus output files instead of a combined file
  --reverse-minus  : reverse the gene order within each contig in the minus strand output,
                     so genes read 5'→3' along the minus strand (only meaningful with
                     --strand split or --strand mark)
  --unmatched FILE : write a TSV report of every CDS in the GFF that has no entry in
                     the presence-absence file (genome, locus_tag, product, contig, strand)

Usage:
    python gene_order_from_pangenome.py \\
        --presence-absence gene_presence_absence.csv \\
        --gff-dir /path/to/gff/files \\
        [--output results.txt] \\
        [--unknown skip|TOKEN] \\
        [--contig-sep |] \\
        [--strand mark|split] \\
        [--reverse-minus] \\
        [--unmatched unmatched.tsv]
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
    Parse a GFF3 file and return a list of (locus_tag, seqid, strand, product) tuples
    for CDS features in file order.

    Uses the 'locus_tag' attribute as the primary identifier, falling back to 'ID'.
    The 'product' attribute is used as the human-readable description, falling back to
    'Name', then to an empty string if neither is present.
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
            product = attrs.get('product') or attrs.get('Name', '')
            if locus_tag and locus_tag not in seen:
                seen.add(locus_tag)
                entries.append((locus_tag, seqid, strand, product))

    return entries


def find_gff(gff_dir, genome_name):
    """Return the GFF file path for a genome, or None if not found."""
    for ext in ('.gff', '.gff3'):
        candidate = gff_dir / f'{genome_name}{ext}'
        if candidate.is_file():
            return candidate
    return None


def process_genome(entries, gene_lookup, unknown_token, mark_strand):
    """
    Convert (locus_tag, seqid, strand, product) entries into gene-family name lists.

    Returns:
        all_fragments   (list[list[str]]) : All genes in GFF order, split by contig.
        plus_fragments  (list[list[str]]) : + strand genes in GFF order, split by contig.
        minus_fragments (list[list[str]]) : - strand genes in GFF order, split by contig.
        unmatched   (list[dict]) : CDS entries with no gene-family in the PA file.
                                    Each dict has keys: locus_tag, product, seqid, strand.
    """
    all_fragments = [[]]
    plus_fragments = [[]]
    minus_fragments = [[]]
    unmatched = []
    prev_seqid = None

    for lt, seqid, strand, product in entries:
        if prev_seqid is not None and seqid != prev_seqid:
            all_fragments.append([])
            plus_fragments.append([])
            minus_fragments.append([])
        prev_seqid = seqid

        family = gene_lookup.get(lt)
        if family:
            label = f'{family}{strand}' if mark_strand else family
        else:
            # Always record the unmatched gene for the QC report
            unmatched.append({
                'locus_tag': lt,
                'product':   product,
                'seqid':     seqid,
                'strand':    strand,
            })
            if unknown_token == 'skip':
                continue  # omit from gene lists but still recorded above
            label = f'{unknown_token}{strand}' if mark_strand else unknown_token

        all_fragments[-1].append(label)
        if strand == '+':
            plus_fragments[-1].append(label)
        else:
            minus_fragments[-1].append(label)

    return all_fragments, plus_fragments, minus_fragments, unmatched


def reverse_fragments(fragments):
    """
    Reverse gene order within each contig fragment independently.

    Empty fragments are preserved so fragment numbering remains aligned with GFF contig order.
    """
    return [list(reversed(fragment)) for fragment in fragments]


def write_output_sequences(fh, genome, fragments, split_by_fragment):
    """
    Write output records in FASTA-like format.

    If split_by_fragment is True, writes one record per fragment named genome.X.
    Otherwise writes one combined record named genome with all fragment genes concatenated.
    Empty fragments are written as blank sequences to preserve fragment indexing.
    """
    if split_by_fragment:
        for idx, genes in enumerate(fragments, start=1):
            fh.write(f'>{genome}.{idx}\n')
            fh.write(','.join(genes) + '\n')
        return

    all_genes = [gene for fragment in fragments for gene in fragment]
    fh.write(f'>{genome}\n')
    fh.write(','.join(all_genes) + '\n')


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
        '--contig-sep', nargs='?', const=True, metavar='SEP',
        help=(
            'Enable fragment-level output: each contig/scaffold is written as its own '
            'sequence entry named "genome.X" (X = fragment number in GFF order). '
            'For backward compatibility, an optional SEP value can still be provided '
            '(e.g. "--contig-sep |"), but it is ignored. Off by default.'
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
    parser.add_argument(
        '--unmatched', metavar='FILE',
        help=(
            'Write a TSV report of every CDS in the GFF files that has no corresponding '
            'entry in the presence-absence file. Columns: genome, locus_tag, product, '
            'contig, strand. Useful for checking what has been excluded from the output.'
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
    split_by_fragment = args.contig_sep is not None
    if isinstance(args.contig_sep, str):
        print(
            f'Warning: --contig-sep value "{args.contig_sep}" is ignored; '
            'fragment-level output is controlled by the flag itself.',
            file=sys.stderr
        )

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

    # Open unmatched report if requested
    unmatched_fh = None
    if args.unmatched:
        unmatched_fh = open(args.unmatched, 'w', encoding='utf-8')
        unmatched_fh.write('genome\tlocus_tag\tproduct\tcontig\tstrand\n')
        print(f'Unmatched gene report: {args.unmatched}', file=sys.stderr)

    total_unmatched = 0

    try:
        for genome in genome_names:
            gff_path = find_gff(gff_dir, genome)
            if gff_path is None:
                print(f'Warning: no GFF file found for genome "{genome}"', file=sys.stderr)
                continue

            print(f'Processing {genome} ({gff_path.name})', file=sys.stderr)
            entries = parse_gff_order(gff_path)
            all_fragments, plus_fragments, minus_fragments, unmatched = process_genome(
                entries, gene_lookup, args.unknown, mark_strand
            )

            if args.reverse_minus:
                minus_fragments = reverse_fragments(minus_fragments)

            if split_strand:
                write_output_sequences(plus_fh, genome, plus_fragments, split_by_fragment)
                write_output_sequences(minus_fh, genome, minus_fragments, split_by_fragment)
            else:
                write_output_sequences(main_fh, genome, all_fragments, split_by_fragment)

            # Write unmatched genes for this genome
            if unmatched_fh:
                for rec in unmatched:
                    unmatched_fh.write(
                        f'{genome}\t{rec["locus_tag"]}\t{rec["product"]}'
                        f'\t{rec["seqid"]}\t{rec["strand"]}\n'
                    )
            total_unmatched += len(unmatched)

    finally:
        if main_fh and args.output:
            main_fh.close()
        if plus_fh:
            plus_fh.close()
        if minus_fh:
            minus_fh.close()
        if unmatched_fh:
            unmatched_fh.close()

    if args.unmatched:
        print(f'  {total_unmatched} unmatched CDS entries written to {args.unmatched}',
              file=sys.stderr)
    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()
