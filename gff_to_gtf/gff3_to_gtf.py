#!/usr/bin/env python3
"""Convert a gene -> transcript -> exon GFF3 into an exon-only GTF for STAR.

Example (no extra Python packages are needed):
    python3 gff3_to_gtf.py annotation.gff3 annotation.gtf --fasta genome.fa

The output contains exon rows with gene_id, transcript_id, and gene_name.
STAR uses these rows for transcript models, splice junctions, and gene counting.
Exons include any annotated untranslated regions (UTRs), not just coding bases.
CDS, gene, transcript, and separate UTR rows are deliberately not copied: this
is a STAR-focused conversion, not a lossless conversion of every GFF3 feature.

Supported input: plain-text GFF3 with explicit genes, transcripts, and exons.
Each transcript must belong to exactly one gene. Shared exons are supported.
Exons directly attached to genes, missing parents, and missing exon annotations
cause an error rather than guessed transcript models. Input order does not matter.

Other limitations:
    - Exons must have a known strand (+ or -); . and ? are rejected.
    - Gene features must be named gene, ncRNA_gene (Ensembl), or pseudogene. Transcript features must
      be named transcript, end in RNA, or end in _transcript. Other naming
      conventions or parent hierarchies may require changes to this script.
    - Only gene_id, transcript_id, and gene_name are written as attributes.
      Product descriptions, database links, biotypes, and other metadata are
      omitted. If a gene has no Name, its ID is used as gene_name.
    - Spaces, semicolons, = and quotes in gene names stay percent-encoded as
      in GFF3 (e.g. CYCA2%3B3), because STAR would cut the name at them. IDs
      containing these characters are rejected instead.
    - This is not a complete GFF3 validator. For another annotation source,
      inspect its structure and independently validate the converted output.

Format references:
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    https://github.com/alexdobin/STAR/blob/master/source/parametersDefault
"""

import argparse
from pathlib import Path
from urllib.parse import unquote


def read_gff3(filename):
    """Read one feature at a time, returning its line number and nine columns."""
    with open(filename, encoding="utf-8") as gff3_file:
        for line_number, line in enumerate(gff3_file, start=1):
            # A GFF3 can contain DNA after this marker. It is not annotation.
            if line.strip() == "##FASTA":
                break
            if line.startswith("#") or not line.strip():
                continue

            columns = line.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"Line {line_number}: expected 9 tab-separated columns.")
            yield line_number, columns


def read_attributes(text):
    """Turn 'ID=tx1;Parent=gene1;' into a dictionary of attribute lists."""
    attributes = {}
    for entry in text.split(";"):
        if not entry or entry == ".":
            continue
        key, value = entry.split("=", 1)
        if key in attributes:
            raise ValueError(f"Repeated GFF3 attribute: {key}")
        # Split BEFORE decoding: %2C means a literal comma inside an ID.
        attributes[key] = [unquote(item) for item in value.split(",")]
    return attributes


def single_attribute(attributes, key, line_number):
    """Get one required value, with a useful error if it is absent or ambiguous."""
    values = attributes.get(key, [])
    if len(values) != 1 or not values[0]:
        raise ValueError(f"Line {line_number}: expected exactly one {key} value.")
    return values[0]


def read_gene_and_transcript_ids(filename):
    """First pass: learn which gene each transcript belongs to."""
    gene_names = {}
    transcript_to_gene = {}
    seen_ids = set()

    for line_number, columns in read_gff3(filename):
        feature_type = columns[2]
        is_gene = feature_type in ("gene", "ncRNA_gene", "pseudogene")
        is_transcript = (
            feature_type in ("transcript", "pseudogenic_transcript")
            or feature_type.endswith("RNA")
            or feature_type.endswith("_transcript")
        )
        if not is_gene and not is_transcript:
            continue

        attributes = read_attributes(columns[8])
        feature_id = single_attribute(attributes, "ID", line_number)
        if feature_id in seen_ids:
            raise ValueError(f"Line {line_number}: repeated gene/transcript ID {feature_id!r}.")
        seen_ids.add(feature_id)

        if is_gene:
            # A gene ID is a useful fallback when no readable gene name exists.
            gene_names[feature_id] = attributes.get("Name", [feature_id])[0]
        else:
            gene_id = single_attribute(attributes, "Parent", line_number)
            transcript_to_gene[feature_id] = gene_id

    for transcript_id, gene_id in transcript_to_gene.items():
        if gene_id not in gene_names:
            raise ValueError(f"Transcript {transcript_id!r} refers to missing gene {gene_id!r}.")

    return gene_names, transcript_to_gene


def read_fasta_lengths(filename):
    """Count bases without loading the whole genome into memory."""
    lengths = {}
    sequence_name = None
    with open(filename, encoding="utf-8") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                # STAR uses the first word of each FASTA header as the name.
                words = line[1:].split()
                if not words:
                    raise ValueError("FASTA contains an empty sequence header.")
                sequence_name = words[0]
                if sequence_name in lengths:
                    raise ValueError(f"Duplicate FASTA sequence name: {sequence_name!r}")
                lengths[sequence_name] = 0
            elif line.strip():
                if sequence_name is None:
                    raise ValueError("FASTA sequence appears before its header.")
                lengths[sequence_name] += len("".join(line.split()))
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError("FASTA is empty or contains an empty sequence.")
    return lengths


# STAR turns these into spaces and keeps only the first word of each attribute
# value (see source/GTF.cpp), so 'CYCA2;3' would be read as 'CYCA2'.
STAR_SEPARATORS = (" ", "\t", "\r", "\n", ";", "=", '"')


def quote_gtf(value):
    """Put a value in GTF quotes; reject characters STAR may misinterpret."""
    if not value or any(character in value for character in STAR_SEPARATORS + ("\\",)):
        raise ValueError(f"Cannot safely write this value as a STAR GTF attribute: {value!r}")
    return f'"{value}"'


def encode_gene_name(name):
    """Keep GFF3 percent-encoding for characters STAR splits on: 'CYCA2;3' -> 'CYCA2%3B3'."""
    unsafe = STAR_SEPARATORS + ("\\", "%")  # Encoding % too keeps this reversible.
    return "".join(f"%{ord(character):02X}" if character in unsafe else character for character in name)


def convert_gff3(input_file, output_file, fasta_file=None):
    """Convert explicit exons, keeping their coordinates and biological IDs."""
    gene_names, transcript_to_gene = read_gene_and_transcript_ids(input_file)
    fasta_lengths = None
    if fasta_file is not None:
        fasta_lengths = read_fasta_lengths(fasta_file)

    # Keep the output in memory until validation succeeds. For this annotation,
    # this is only tens of MB and prevents writing a partial GTF after an error.
    output_lines = []
    used_genes = set()
    used_transcripts = set()

    # Second pass: follow exon -> transcript -> gene to create GTF attributes.
    for line_number, columns in read_gff3(input_file):
        if columns[2] != "exon":
            continue
        attributes = read_attributes(columns[8])
        transcript_ids = attributes.get("Parent", [])
        if not transcript_ids:
            raise ValueError(f"Line {line_number}: exon has no Parent transcript.")

        sequence_name = unquote(columns[0])
        start = int(columns[3])
        end = int(columns[4])
        if start < 1 or end < start:
            raise ValueError(f"Line {line_number}: invalid exon coordinates {start}-{end}.")
        if columns[6] not in ("+", "-"):
            raise ValueError(f"Line {line_number}: exon needs a known strand (+ or -).")
        if any(character.isspace() for character in sequence_name):
            raise ValueError(f"Line {line_number}: sequence name contains whitespace.")
        if fasta_lengths is not None:
            if sequence_name not in fasta_lengths:
                raise ValueError(f"Line {line_number}: {sequence_name!r} is absent from FASTA.")
            if end > fasta_lengths[sequence_name]:
                raise ValueError(f"Line {line_number}: exon extends beyond {sequence_name!r}.")

        # Both formats use 1-based, inclusive coordinates. Do not subtract 1.
        columns[0] = sequence_name
        columns[7] = "."  # Reading-frame phase is not applicable to exon rows.
        for transcript_id in transcript_ids:
            if transcript_id not in transcript_to_gene:
                raise ValueError(
                    f"Line {line_number}: exon parent {transcript_id!r} is not a supported "
                    "transcript with a gene parent. No transcript will be invented."
                )
            gene_id = transcript_to_gene[transcript_id]
            columns[8] = (
                f"gene_id {quote_gtf(gene_id)}; "
                f"transcript_id {quote_gtf(transcript_id)}; "
                f"gene_name {quote_gtf(encode_gene_name(gene_names[gene_id]))};"
            )
            # A shared exon gets one row per parent transcript.
            output_lines.append("\t".join(columns) + "\n")
            used_genes.add(gene_id)
            used_transcripts.add(transcript_id)

    if not output_lines:
        raise ValueError("No exons found. CDS rows cannot substitute for full exons.")
    missing_transcripts = set(transcript_to_gene) - used_transcripts
    if missing_transcripts:
        example = sorted(missing_transcripts)[0]
        raise ValueError(
            f"{len(missing_transcripts)} transcripts have no explicit exons (e.g. {example!r}). "
            "Supply an annotation with exon features."
        )

    # Exclusive creation protects existing files, including the input files.
    with open(output_file, "x", encoding="utf-8") as gtf_file:
        gtf_file.write("# Exon-only GTF for STAR; generated by gff3_to_gtf.py\n")
        gtf_file.writelines(output_lines)

    print(f"Wrote {len(output_lines):,} exon rows to {output_file}")
    print(f"Included {len(used_genes):,} genes and {len(used_transcripts):,} transcripts.")
    if fasta_lengths is not None:
        print("All exon sequence names and coordinates match the FASTA.")
    unused_genes = set(gene_names) - used_genes
    if unused_genes:
        print(f"Note: {len(unused_genes):,} genes without exon-bearing transcripts were omitted.")


def main():
    """Read the filenames supplied on the command line, then run the conversion."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gff3", type=Path, help="Input annotation (plain-text GFF3)")
    parser.add_argument("gtf", type=Path, help="Output exon-only GTF (must not already exist)")
    parser.add_argument("--fasta", type=Path, help="Optional genome FASTA for checking exon names and coordinates")
    args = parser.parse_args()
    if args.gtf.exists():
        parser.error(f"Output already exists: {args.gtf}. Choose a new filename.")
    try:
        convert_gff3(args.gff3, args.gtf, args.fasta)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
