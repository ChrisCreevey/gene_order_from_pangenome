# gene_order_from_pangenome

A command-line tool that reconstructs the chromosomal gene order of each genome in a
pangenome, using gene-family names rather than locus tags.

It combines a **Roary-style presence-absence table** with the **GFF annotation file for
each genome** to produce, per genome, a comma-separated list of gene-family names in the
exact order the genes appear in the genome.

---

## Background

Pangenome tools such as [Roary](https://sanger-pathogens.github.io/Roary/) cluster
protein-coding genes across many genomes into gene families and report which locus tag(s)
from each genome belong to each family. However, the presence-absence table does not
preserve genomic order. This tool restores that order by walking each genome's GFF file
and replacing every locus tag with its gene-family name.

---

## Requirements

- Python 3.7 or later
- No third-party packages — only the Python standard library is used

---

## Inputs

### 1. Presence-absence CSV (`--presence-absence`)

The `gene_presence_absence.csv` file produced by Roary (or a compatible tool).

- **Column 1:** Gene-family name
- **Columns 2–14:** Metadata fields (Non-unique gene name, Annotation, No. isolates, etc.)
- **Columns 15+:** One column per genome; each cell contains the locus tag(s) from that
  genome that belong to the gene family. Paralogs (multiple copies) are separated by a
  tab character within the cell, as per standard Roary output.

Example header:

```
"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences",
"Avg sequences per isolate","Genome Fragment","Order within Fragment",
"Accessory Fragment","Accessory Order with Fragment","QC",
"Min group size nuc","Max group size nuc","Avg group size nuc",
"1004153.3","1005563.3","1005566.3", ...
```

### 2. GFF annotation files (`--gff-dir`)

One GFF3 file per genome, located in a single directory. Files must be named after the
genome identifier exactly as it appears in the presence-absence table header:

```
<genome_name>.gff      e.g.  1004153.3.gff
<genome_name>.gff3     e.g.  1004153.3.gff3
```

Only `CDS` feature rows are read. The gene identifier is taken from the `locus_tag`
attribute, falling back to `ID` if `locus_tag` is absent. If the GFF contains an
embedded FASTA section (lines beginning with `>`), it is ignored.

---

## Usage

```bash
python gene_order_from_pangenome.py \
    --presence-absence gene_presence_absence.csv \
    --gff-dir /path/to/gff/files \
    --output results.txt
```

Output is written to stdout if `--output` is omitted.

### All options

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--presence-absence` | `-p` | *(required)* | Path to the Roary presence-absence CSV |
| `--gff-dir` | `-g` | *(required)* | Directory containing the per-genome GFF files |
| `--output` | `-o` | stdout | Output file path. With `--strand split`, used as the base name for `.plus` and `.minus` files |
| `--unknown` | `-u` | `*` | Token written for CDS features not found in the presence-absence table. Pass `skip` to omit them entirely |
| `--meta-cols` | `-m` | `14` | Number of metadata columns before the genome columns (default matches standard Roary output) |
| `--contig-sep` | | *(off)* | Token inserted at every contig/scaffold boundary (e.g. `\|`) to indicate genes are on different genome fragments |
| `--strand` | | *(off)* | `mark`: append `+` or `-` to each gene-family name. `split`: write separate `.plus` and `.minus` output files |
| `--reverse-minus` | | *(off)* | Reverse minus strand gene order within each contig so genes read 5'→3' along the minus strand (see below) |
| `--unmatched` | | *(off)* | Path for a TSV report of every CDS in the GFF with no entry in the presence-absence file |

---

## Output format

The output contains one block per genome:

```
>1004153.3
ppnP,lptC,ptsN,*,*,ribE,yhcB,fadJ,fadJ,*,mraZ,...
```

- The header line (`>genome_name`) identifies the genome.
- The following line is a comma-separated list of gene-family names in chromosomal order.
- Genes not assigned to any family in the presence-absence table are represented by the
  `--unknown` token (default `*`).
- Paralogs appear at each of their respective positions in the GFF; for example, two
  copies of `fadJ` will appear twice, at the positions of each copy.

### With `--contig-sep`

A separator token is inserted wherever the sequence ID changes, marking the boundary
between contigs, scaffolds, or chromosomes:

```
>1004153.3
ppnP,lptC,|,ptsN,ribE,yhcB,|,fadJ,fadJ,mraZ,...
```

The two genes either side of `|` are on different genome fragments and are not physically
adjacent. The separator appears in all output modes (combined, `mark`, and `split`).

### With `--strand mark`

A `+` or `-` character is appended to each gene-family name to indicate which strand the
gene is encoded on:

```
>1004153.3
ppnP+,lptC-,ptsN+,*+,ribE+,fadJ-,fadJ-,mraZ+,...
```

### With `--strand split`

Instead of a single combined file, two files are written — one per strand — each
containing only the genes from that strand in their GFF order. Given
`--output results.txt`, the output files are:

```
results.plus.txt   →  genes on the + strand
results.minus.txt  →  genes on the - strand
```

Each file uses the same `>genome_name` / comma-separated-list format. Contig separators
(`--contig-sep`) are included in split files if requested.

### With `--reverse-minus`

GFF files list all genes in reference coordinate order (lowest to highest start position),
regardless of strand. This means:

- **`+` strand genes** — GFF order already reads 5'→3' ✓
- **`−` strand genes** — GFF order reads 3'→5', because the gene with the lowest
  coordinate is at the 3' end of the minus strand ✗

Passing `--reverse-minus` corrects this by reversing the minus strand genes **within
each contig independently**, so they read 5'→3' along the minus strand while the contigs
themselves remain in their original order.

This flag affects the minus strand genes in all output modes — the combined file
(`--strand mark`), the `.minus.txt` split file, or the default combined output.

### Unmatched gene report (`--unmatched`)

Any CDS in a GFF file whose `locus_tag` has no entry in the presence-absence table is
written to a tab-separated report for later inspection. This is independent of the
`--unknown` setting — genes are recorded even if they are being skipped in the main
output. The report has five columns:

```
genome    locus_tag         product                           contig           strand
1004153.3 BFGHOIAM_00002    hypothetical protein              AFET01000001     +
1004153.3 BFGHOIAM_00003    hypothetical protein              AFET01000001     +
1004153.3 BFGHOIAM_00007    putative membrane protein         AFET01000001     -
```

A summary count of unmatched entries is printed to stderr when the run completes.

---

## Examples

**Basic run, output to file:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    -o gene_order.txt
```

**Skip genes not in the presence-absence table:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    --unknown skip \
    -o gene_order_core_only.txt
```

**Print to terminal (useful for a quick check):**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ | head -4
```

**Mark contig boundaries with `|`:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    --contig-sep '|' \
    -o gene_order.txt
```

**Append strand to every gene-family name:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    --strand mark \
    -o gene_order_stranded.txt
```

**Write separate files for each strand:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    --strand split \
    -o gene_order.txt
# produces gene_order.plus.txt and gene_order.minus.txt
```

**Combine contig boundaries and strand splitting:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    --contig-sep '|' \
    --strand split \
    -o gene_order.txt
```

**Generate an unmatched gene report for QC:**
```bash
python gene_order_from_pangenome.py \
    -p gene_presence_absence.csv \
    -g gff/ \
    -o gene_order.txt \
    --unmatched unmatched_genes.tsv
```

---

## Memory requirements

The presence-absence CSV is read line-by-line, so the raw file size is not the limiting
factor. Memory is dominated by the `gene_lookup` dictionary, which holds a
`locus_tag → gene_family` entry for every gene assignment in the table.

A rough estimate can be obtained from the dimensions of your dataset:

```
total entries  =  number of genomes  ×  average genes per genome
RAM (GB)       ≈  total entries  ×  115  /  1,000,000,000
```

Add ~0.5 GB for Python and OS overhead. For typical bacterial datasets (3,000–6,000
genes per genome), **16 GB is sufficient for most large runs**. If you are unsure, check
the number of columns and rows with:

```bash
head -1 gene_presence_absence.csv | tr ',' '\n' | wc -l   # genomes (subtract 14)
wc -l < gene_presence_absence.csv                          # gene families (subtract 1)
```

---

## Notes

- Progress and warnings are written to **stderr**, so they do not mix with the output
  when piping to a file.
- A warning is printed for any genome in the presence-absence table for which no GFF
  file is found; that genome is skipped.
- The tool uses the GFF file order as-is. Genes on different contigs/chromosomes appear
  in the order the contigs are listed in the GFF.
