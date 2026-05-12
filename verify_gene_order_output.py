#!/usr/bin/env python3
"""
Sample-based independent verification for gene_order_from_pangenome output.

This tool checks a user-specified sample of output records/positions and verifies:
  1) The sampled position maps to a real CDS locus in the genome GFF order.
  2) The output gene-family token matches the family assignment from the
     presence-absence pangenome table for that locus.
"""

import argparse
import csv
import random
import sys
from datetime import datetime, timezone
from pathlib import Path


def parse_presence_absence_by_genome(pa_path, n_meta_cols):
    """Return (genome_names, per_genome_lookup[genome][locus_tag] -> family)."""
    per_genome = {}
    with open(pa_path, newline='', encoding='utf-8') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        genome_names = [c.strip() for c in header[n_meta_cols:]]
        per_genome = {g: {} for g in genome_names}

        for row in reader:
            if not row:
                continue
            family = row[0].strip()
            for idx, genome in enumerate(genome_names, start=n_meta_cols):
                if idx >= len(row):
                    continue
                cell = row[idx].strip()
                if not cell:
                    continue
                for locus in cell.split('\t'):
                    locus = locus.strip()
                    if locus:
                        per_genome[genome][locus] = family
    return genome_names, per_genome


def parse_gff_loci_by_fragment(gff_path):
    """
    Return ordered locus tags grouped by contig fragment:
    [[(locus_tag, strand), ...], ...]
    """
    seen = set()
    fragments = [[]]
    prev_seqid = None

    with open(gff_path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            if line.startswith('>'):
                break
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] != 'CDS':
                continue

            seqid = cols[0]
            strand = cols[6]
            attrs = {}
            for part in cols[8].split(';'):
                if '=' in part:
                    k, _, v = part.partition('=')
                    attrs[k.strip()] = v.strip()

            locus = attrs.get('locus_tag') or attrs.get('ID', '').split('.')[0]
            if not locus or locus in seen:
                continue

            if prev_seqid is not None and seqid != prev_seqid:
                fragments.append([])
            prev_seqid = seqid

            seen.add(locus)
            fragments[-1].append((locus, strand))
    return fragments


def find_gff(gff_dir, genome):
    for ext in ('.gff', '.gff3'):
        p = gff_dir / f'{genome}{ext}'
        if p.is_file():
            return p
    return None


def parse_output_records(result_path):
    """Parse FASTA-like output records into {header: [tokens]}."""
    records = {}
    current = None
    chunks = []
    with open(result_path, encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current is not None:
                    joined = ''.join(chunks)
                    tokens = [t.strip() for t in joined.split(',')] if joined else []
                    records[current] = [t for t in tokens if t != '']
                current = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if current is not None:
            joined = ''.join(chunks)
            tokens = [t.strip() for t in joined.split(',')] if joined else []
            records[current] = [t for t in tokens if t != '']
    return records


def normalize_token(token, strip_strand_mark):
    token = token.strip()
    if strip_strand_mark and token.endswith(('+', '-')):
        return token[:-1]
    return token


def resolve_header(header, genome_set):
    """
    Resolve record header to (genome, fragment_index|None).
    fragment_index is 1-based when present.
    """
    if header in genome_set:
        return header, None
    if '.' in header:
        base, frag_str = header.rsplit('.', 1)
        if base in genome_set and frag_str.isdigit():
            return base, int(frag_str)
    return None, None


def sampled_positions(rng, n_items, requested):
    if n_items <= 0:
        return []
    k = min(n_items, requested)
    return sorted(rng.sample(range(n_items), k=k))


def main():
    ap = argparse.ArgumentParser(
        description=(
            'Independently verify a user-specified sample of gene-order output entries '
            'against genome GFF order and pangenome family assignments.'
        )
    )
    ap.add_argument('--presence-absence', '-p', required=True, metavar='FILE')
    ap.add_argument('--gff-dir', '-g', required=True, metavar='DIR')
    ap.add_argument('--result-file', '-r', required=True, metavar='FILE')
    ap.add_argument('--report', '-o', metavar='FILE',
                    help='Write report to this file (default: stdout)')
    ap.add_argument('--sample-records', type=int, default=5, metavar='N',
                    help='Number of output records (headers) to sample (default: 5)')
    ap.add_argument('--sample-genes', type=int, default=20, metavar='N',
                    help='Number of gene positions to sample per selected record (default: 20)')
    ap.add_argument('--seed', type=int, default=1,
                    help='Random seed for reproducible sampling (default: 1)')
    ap.add_argument('--meta-cols', '-m', type=int, default=14,
                    help='Metadata column count before genome columns in PA CSV (default: 14)')
    ap.add_argument('--unknown', '-u', default='*',
                    help='Unknown token expected in output for unmatched loci (default: *)')
    ap.add_argument('--strip-strand-mark', action='store_true',
                    help='Strip trailing + / - from output tokens before comparison')
    args = ap.parse_args()

    pa_path = Path(args.presence_absence)
    gff_dir = Path(args.gff_dir)
    result_path = Path(args.result_file)
    if not pa_path.is_file():
        sys.exit(f'Error: presence-absence file not found: {pa_path}')
    if not gff_dir.is_dir():
        sys.exit(f'Error: GFF directory not found: {gff_dir}')
    if not result_path.is_file():
        sys.exit(f'Error: result file not found: {result_path}')
    if args.sample_records <= 0:
        sys.exit('Error: --sample-records must be > 0')
    if args.sample_genes <= 0:
        sys.exit('Error: --sample-genes must be > 0')

    genome_names, per_genome_lookup = parse_presence_absence_by_genome(pa_path, args.meta_cols)
    genome_set = set(genome_names)
    records = parse_output_records(result_path)
    headers = sorted(records.keys())
    if not headers:
        sys.exit('Error: no output records found in result file')

    rng = random.Random(args.seed)
    selected_headers = rng.sample(headers, k=min(args.sample_records, len(headers)))

    gff_cache = {}
    total_checked = 0
    passed = 0
    failed = 0
    notes = []
    per_record = []

    for header in selected_headers:
        genome, frag_idx = resolve_header(header, genome_set)
        observed = records[header]
        if genome is None:
            notes.append(f'SKIP {header}: header not mappable to a genome in PA table')
            continue

        if genome not in gff_cache:
            gff_path = find_gff(gff_dir, genome)
            if gff_path is None:
                notes.append(f'SKIP {header}: GFF not found for genome {genome}')
                gff_cache[genome] = None
            else:
                gff_cache[genome] = parse_gff_loci_by_fragment(gff_path)

        fragments = gff_cache[genome]
        if fragments is None:
            continue

        if frag_idx is None:
            expected_loci = [item for fragment in fragments for item in fragment]
        else:
            if frag_idx < 1 or frag_idx > len(fragments):
                notes.append(
                    f'SKIP {header}: fragment index {frag_idx} out of range (1..{len(fragments)})'
                )
                continue
            expected_loci = fragments[frag_idx - 1]

        expected_tokens = []
        genome_map = per_genome_lookup.get(genome, {})
        for locus, _strand in expected_loci:
            expected_tokens.append(genome_map.get(locus, args.unknown))

        n_common = min(len(observed), len(expected_tokens))
        positions = sampled_positions(rng, n_common, args.sample_genes)

        mismatches = 0
        checks = []
        for pos in positions:
            locus_tag, strand = expected_loci[pos]
            observed_token = normalize_token(observed[pos], args.strip_strand_mark)
            expected_token = expected_tokens[pos]
            ok = observed_token == expected_token
            total_checked += 1
            if ok:
                passed += 1
            else:
                failed += 1
                mismatches += 1
            checks.append({
                'position': pos + 1,
                'locus_tag': locus_tag,
                'strand': strand,
                'observed': observed_token,
                'expected': expected_token,
                'match': ok,
            })

        length_match = len(observed) == len(expected_tokens)
        per_record.append({
            'header': header,
            'genome': genome,
            'fragment': frag_idx,
            'observed_len': len(observed),
            'expected_len': len(expected_tokens),
            'length_match': length_match,
            'sampled': len(positions),
            'mismatches': mismatches,
            'checks': checks,
        })

    out = open(args.report, 'w', encoding='utf-8') if args.report else sys.stdout
    try:
        out.write('Gene order verification report\n')
        out.write('==============================\n')
        out.write(f'Generated (UTC): {datetime.now(timezone.utc).isoformat()}\n')
        out.write(f'Presence-absence: {pa_path}\n')
        out.write(f'GFF directory: {gff_dir}\n')
        out.write(f'Result file: {result_path}\n')
        out.write(f'Sample records requested: {args.sample_records}\n')
        out.write(f'Sample genes per record requested: {args.sample_genes}\n')
        out.write(f'Random seed: {args.seed}\n')
        out.write(f'Unknown token: {args.unknown}\n')
        out.write(f'Strip strand mark: {args.strip_strand_mark}\n\n')

        out.write('Summary\n')
        out.write('-------\n')
        out.write(f'Sampled output records: {len(per_record)}\n')
        out.write(f'Sampled positions checked: {total_checked}\n')
        out.write(f'Matches: {passed}\n')
        out.write(f'Mismatches: {failed}\n')
        if total_checked:
            out.write(f'Pass rate: {(passed / total_checked) * 100:.2f}%\n')
        out.write('\n')

        if notes:
            out.write('Notes / skipped records\n')
            out.write('-----------------------\n')
            for note in notes:
                out.write(f'- {note}\n')
            out.write('\n')

        out.write('Per-record verification\n')
        out.write('-----------------------\n')
        for rec in per_record:
            frag_text = f'.{rec["fragment"]}' if rec['fragment'] is not None else ''
            out.write(f'Record: {rec["header"]} (resolved: {rec["genome"]}{frag_text})\n')
            out.write(
                f'  Length check: observed={rec["observed_len"]}, '
                f'expected={rec["expected_len"]}, '
                f'match={rec["length_match"]}\n'
            )
            out.write(
                f'  Sampled positions: {rec["sampled"]}, '
                f'mismatches in sample: {rec["mismatches"]}\n'
            )
            for c in rec['checks']:
                out.write(
                    '    '
                    f'pos={c["position"]} locus={c["locus_tag"]} strand={c["strand"]} '
                    f'observed={c["observed"]} expected={c["expected"]} '
                    f'match={c["match"]}\n'
                )
            out.write('\n')
    finally:
        if args.report:
            out.close()

    if failed > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
