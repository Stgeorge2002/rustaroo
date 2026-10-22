#!/usr/bin/env python3
"""Look each reference's gene name up in UniProt and compare lengths.

The sequence is not sent and not translated. The FASTA header is the query.
UniProt returns an amino-acid length; the base-pair length is that number
times 3, plus 3 where the reference counts the stop codon. A name UniProt
does not know is left out of the match rate.
"""

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def read_lengths(path):
    name = None
    length = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        yield name, length


def request(url, retries=5):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "rustroo-uniprot-lengths/1.0"}
            )
            with urllib.request.urlopen(req, timeout=60) as response:
                return response.read().decode()
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504) or attempt + 1 == retries:
                detail = error.read().decode(errors="replace")[:300]
                raise RuntimeError("HTTP %s: %s" % (error.code, detail))
        except urllib.error.URLError:
            if attempt + 1 == retries:
                raise
        time.sleep(2 ** attempt)
    raise RuntimeError(url)


def lookup(name, taxon):
    quoted = name.replace('"', "")
    query = '(gene:"%s") AND taxonomy_id:%s' % (quoted, taxon)
    encoded = urllib.parse.urlencode({
        "query": query,
        "fields": "accession,length,gene_names",
        "format": "tsv",
        "size": "25",
    })
    lines = request(SEARCH + "?" + encoded).splitlines()[1:]
    lengths = []
    accessions = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2 or not parts[1].isdigit():
            continue
        accessions.append(parts[0])
        lengths.append(int(parts[1]))
    if not lengths:
        return {"status": "absent"}
    modal = Counter(lengths).most_common(1)[0][0]
    return {
        "status": "found",
        "uniprot": accessions[0],
        "uniprot_aa": modal,
        "n_entries": len(lengths),
        "n_lengths": len(set(lengths)),
    }


def agrees(reference_bp, amino_acids):
    return reference_bp in (amino_acids * 3, amino_acids * 3 + 3)


def safe_cache_name(name):
    slug = "".join(character if character.isalnum() else "_" for character in name)[:80]
    digest = hashlib.sha256(name.encode("utf-8")).hexdigest()[:12]
    return "%s_%s" % (slug, digest)


def main():
    root = Path(__file__).resolve().parents[1] / "macs_curve_fine2"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rustroo", type=Path, default=root / "n1500" / "rustroo" / "genes.fa")
    parser.add_argument(
        "--panaroo",
        type=Path,
        default=root / "n1500" / "panaroo" / "pan_genome_reference.fa",
    )
    parser.add_argument("--output", type=Path, default=root / "uniprot_lengths")
    parser.add_argument("--taxon", type=int, default=1313)
    parser.add_argument("--limit", type=int, default=0, help="stop after this many distinct names; 0 means all")
    args = parser.parse_args()

    records = []
    for tool, path in (("rustroo", args.rustroo), ("panaroo", args.panaroo)):
        for name, bp in read_lengths(path):
            records.append((tool, name, bp))
    names = []
    seen = set()
    for _tool, name, _bp in records:
        if name not in seen:
            seen.add(name)
            names.append(name)
    if args.limit:
        names = names[:args.limit]
    wanted = set(names)

    cache_dir = args.output / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    found = {}
    pending = 0
    for index, name in enumerate(names, 1):
        path = cache_dir / (safe_cache_name(name) + ".json")
        if path.is_file():
            found[name] = json.loads(path.read_text())
            continue
        found[name] = lookup(name, args.taxon)
        path.write_text(json.dumps(found[name]))
        pending += 1
        if pending % 50 == 0:
            print("looked up %s/%s" % (index, len(names)))
        time.sleep(0.2)
    print("names %s, new lookups %s" % (len(names), pending))

    report = args.output / "lengths.tsv"
    counts = {}
    with report.open("w") as handle:
        handle.write("\t".join([
            "tool", "gene", "reference_bp", "status", "uniprot", "uniprot_aa",
            "uniprot_bp", "n_entries", "n_lengths", "length_match",
        ]) + "\n")
        for tool, name, bp in records:
            if name not in wanted:
                continue
            row = found[name]
            amino = row.get("uniprot_aa")
            match = row.get("status") == "found" and agrees(bp, amino)
            handle.write("\t".join([
                tool,
                name,
                str(bp),
                row.get("status", ""),
                row.get("uniprot", ""),
                str(amino or ""),
                str(amino * 3 + 3 if amino else ""),
                str(row.get("n_entries", "")),
                str(row.get("n_lengths", "")),
                "yes" if match else "no",
            ]) + "\n")
            bucket = counts.setdefault(tool, Counter())
            bucket["references"] += 1
            bucket[row.get("status", "absent")] += 1
            bucket["match"] += int(match)

    print()
    for tool in ("rustroo", "panaroo"):
        bucket = counts.get(tool, Counter())
        known = bucket["found"]
        rate = 100.0 * bucket["match"] / known if known else 0.0
        print(
            "%s: %s references, %s names in UniProt, %s length matches (%.1f%%)"
            % (tool, bucket["references"], known, bucket["match"], rate)
        )
    print("per-gene table: %s" % report)


if __name__ == "__main__":
    main()
