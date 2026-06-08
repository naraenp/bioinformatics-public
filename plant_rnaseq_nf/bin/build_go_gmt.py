#!/usr/bin/env python3
"""Build a GO gene-set GMT for the enrichment step from Ensembl Plants BioMart.

Pulls the gene -> GO term mapping for a species (default rice, osativa_eg_gene)
and writes it as a GMT (term_id, term_name, gene1, gene2, ...). Because the gene
identifiers come from the same Ensembl Plants release as the genome/GTF used for
alignment, they match the featureCounts gene_ids without any extra mapping.

Stdlib only. Example:
    build_go_gmt.py --out data/real/rice_go.gmt
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="{schema}" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="{dataset}" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="go_id"/>
    <Attribute name="name_1006"/>
    <Attribute name="namespace_1003"/>
  </Dataset>
</Query>"""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--dataset", default="osativa_eg_gene")
    p.add_argument("--schema", default="plants_mart")
    p.add_argument("--host", default="https://plants.ensembl.org/biomart/martservice")
    p.add_argument("--namespace", default="biological_process",
                   help="GO namespace to keep (or 'all')")
    p.add_argument("--min-size", type=int, default=5)
    p.add_argument("--max-size", type=int, default=1000)
    p.add_argument("--timeout", type=int, default=300)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    url = f"{args.host}?query={quote(XML.format(schema=args.schema, dataset=args.dataset))}"
    print(f"querying BioMart: {args.dataset} ...", file=sys.stderr)
    with urlopen(url, timeout=args.timeout) as resp:
        text = resp.read().decode("utf-8", "replace")

    names: dict[str, str] = {}
    members: dict[str, set[str]] = defaultdict(set)
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        gene, go_id, go_name, ns = (p.strip() for p in parts[:4])
        if not gene or not go_id:
            continue
        if args.namespace != "all" and ns != args.namespace:
            continue
        names[go_id] = go_name
        members[go_id].add(gene)

    kept = 0
    with open(args.out, "w") as fh:
        for go_id, genes in sorted(members.items()):
            if not (args.min_size <= len(genes) <= args.max_size):
                continue
            fh.write("\t".join([go_id, names.get(go_id, go_id), *sorted(genes)]) + "\n")
            kept += 1

    if kept == 0:
        raise SystemExit("no gene sets written — check dataset/namespace/network")
    print(f"wrote {args.out}: {kept} GO sets "
          f"(size {args.min_size}-{args.max_size}, namespace={args.namespace})")


if __name__ == "__main__":
    main()
