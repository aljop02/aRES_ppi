#!/usr/bin/env python3
"""
interface_contacts.py -- inter-chain interface analysis for aRES / phage-protein complexes.

Reads one or more complex structures (PDB or mmCIF, so both HADDOCK and AlphaFold3
output work) and reports, per structure and aggregated over all of them:

  * inter-chain heavy-atom contacts below a distance cutoff
  * per-residue contact counts for each chain
  * contact distribution split across user-defined domains
  * the "P1 readout": what fraction of the interface on a given chain falls inside
    a region of interest (for this project: T7 gp5 insertion1, residues ~102-122)

Aggregating over many models is the point. A single top-ranked docking pose is weak
evidence; a contact-frequency map over every model in a run is much harder to fool.

Dependencies: numpy only (no Biopython needed).

Usage
-----
  # single model, chains auto-detected
  python3 interface_contacts.py model.cif

  # whole HADDOCK run, with domains defined and the P1 region flagged
  python3 interface_contacts.py haddock_run/*.pdb \
      --domain A:1-150:NTD1 --domain A:151-459:RES \
      --roi B:102-122:insertion1 \
      --tsv results/contacts.tsv

  # explicit chain pair, tighter cutoff
  python3 interface_contacts.py model.cif --chains A B --cutoff 4.5

Caveats
-------
Contact counting is geometry only. It says nothing about whether an interface is
real, and a high contact count in a low-confidence region (check pLDDT / HADDOCK
score separately) means little. Use --min-plddt to drop unreliable residues when
analysing predicted models.
"""

from __future__ import annotations

import argparse
import glob
import sys
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

# Residues are identified by (chain, resseq, icode) so insertion codes don't collide.
ResKey = tuple


@dataclass
class Atoms:
    """Flat arrays of heavy atoms from one structure."""

    chain: np.ndarray  # str
    resseq: np.ndarray  # int
    icode: np.ndarray  # str
    resname: np.ndarray  # str
    name: np.ndarray  # str
    bfac: np.ndarray  # float (pLDDT for AF models)
    xyz: np.ndarray  # (N, 3) float

    def __len__(self) -> int:
        return len(self.chain)

    def mask(self, m: np.ndarray) -> "Atoms":
        return Atoms(
            self.chain[m], self.resseq[m], self.icode[m], self.resname[m],
            self.name[m], self.bfac[m], self.xyz[m],
        )


# Anything that is clearly not protein/nucleic: water, ions, common buffer junk.
# DNA is dropped too -- 1T7P carries a primer/template that must not be docked.
SKIP_RESNAMES = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "MPD", "PEG", "ACT",
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "CD", "NI", "CU",
    "DA", "DT", "DG", "DC", "DU", "DI", "A", "U", "G", "C", "I",
}


def _finish(records: list, source: str) -> Atoms:
    if not records:
        raise ValueError(f"no usable heavy atoms parsed from {source}")
    cols = list(zip(*records))
    return Atoms(
        chain=np.array(cols[0], dtype=object),
        resseq=np.array(cols[1], dtype=int),
        icode=np.array(cols[2], dtype=object),
        resname=np.array(cols[3], dtype=object),
        name=np.array(cols[4], dtype=object),
        bfac=np.array(cols[5], dtype=float),
        xyz=np.array(cols[6], dtype=float),
    )


def parse_pdb(path: str) -> Atoms:
    records = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                if line.startswith("ENDMDL"):
                    break  # first model only
                continue
            element = line[76:78].strip().upper()
            name = line[12:16].strip()
            if element == "H" or (not element and name.startswith("H")):
                continue
            resname = line[17:20].strip().upper()
            if resname in SKIP_RESNAMES:
                continue
            chain = line[21].strip() or "_"
            try:
                resseq = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            try:
                bfac = float(line[60:66])
            except ValueError:
                bfac = 0.0
            records.append((chain, resseq, line[26].strip(), resname, name, bfac, xyz))
    return _finish(records, path)


def parse_cif(path: str) -> Atoms:
    """Minimal mmCIF atom_site reader -- enough for AF3 and PDB-downloaded CIFs."""
    records = []
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()

    i, n = 0, len(lines)
    while i < n:
        if not lines[i].strip().startswith("loop_"):
            i += 1
            continue
        # Collect the column names of this loop.
        j = i + 1
        tags = []
        while j < n and lines[j].strip().startswith("_"):
            tags.append(lines[j].strip().split()[0])
            j += 1
        if not tags or not tags[0].startswith("_atom_site."):
            i = j
            continue

        idx = {t.split(".", 1)[1]: k for k, t in enumerate(tags)}

        def col(*names, default=None):
            for nm in names:
                if nm in idx:
                    return idx[nm]
            return default

        c_group = col("group_PDB")
        c_sym = col("type_symbol")
        c_atom = col("label_atom_id", "auth_atom_id")
        c_comp = col("auth_comp_id", "label_comp_id")
        c_asym = col("auth_asym_id", "label_asym_id")
        c_seq = col("auth_seq_id", "label_seq_id")
        c_ins = col("pdbx_PDB_ins_code")
        c_b = col("B_iso_or_equiv")
        c_x, c_y, c_z = col("Cartn_x"), col("Cartn_y"), col("Cartn_z")
        c_model = col("pdbx_PDB_model_num")
        if None in (c_atom, c_comp, c_asym, c_seq, c_x, c_y, c_z):
            i = j
            continue

        first_model = None
        while j < n:
            raw = lines[j].strip()
            if not raw or raw.startswith(("#", "loop_", "_", "data_")):
                break
            f = raw.split()
            j += 1
            if len(f) < len(tags):
                continue
            if c_group is not None and f[c_group] not in ("ATOM", "HETATM"):
                continue
            if c_model is not None:
                if first_model is None:
                    first_model = f[c_model]
                elif f[c_model] != first_model:
                    break  # first model only
            name = f[c_atom].strip('"')
            if (c_sym is not None and f[c_sym].upper() == "H") or name.startswith("H"):
                continue
            resname = f[c_comp].upper()
            if resname in SKIP_RESNAMES:
                continue
            try:
                resseq = int(f[c_seq])
                xyz = (float(f[c_x]), float(f[c_y]), float(f[c_z]))
            except ValueError:
                continue
            try:
                bfac = float(f[c_b]) if c_b is not None else 0.0
            except ValueError:
                bfac = 0.0
            ins = "" if c_ins is None or f[c_ins] in ("?", ".") else f[c_ins]
            records.append((f[c_asym], resseq, ins, resname, name, bfac, xyz))
        i = j
    return _finish(records, path)


def load(path: str) -> Atoms:
    return parse_cif(path) if path.lower().endswith((".cif", ".mmcif")) else parse_pdb(path)


def contacts(a: Atoms, b: Atoms, cutoff: float) -> set:
    """Residue pairs with at least one heavy-atom pair closer than `cutoff`."""
    if len(a) == 0 or len(b) == 0:
        return set()
    cut2 = cutoff * cutoff
    found = set()
    # Chunked brute force: exact, and fast enough at these sizes.
    for s in range(0, len(a), 512):
        block = a.xyz[s : s + 512]
        d2 = ((block[:, None, :] - b.xyz[None, :, :]) ** 2).sum(-1)
        ia, ib = np.nonzero(d2 < cut2)
        for p, q in zip(ia, ib):
            p += s
            found.add((
                (a.chain[p], int(a.resseq[p]), a.icode[p], a.resname[p]),
                (b.chain[q], int(b.resseq[q]), b.icode[q], b.resname[q]),
            ))
    return found


def parse_region(spec: str, what: str) -> tuple:
    """'A:102-122:insertion1' -> ('A', 102, 122, 'insertion1')"""
    parts = spec.split(":")
    if len(parts) < 2:
        sys.exit(f"bad --{what} '{spec}'; expected CHAIN:START-END[:NAME]")
    chain = parts[0]
    try:
        lo, hi = (int(x) for x in parts[1].split("-"))
    except ValueError:
        sys.exit(f"bad range in --{what} '{spec}'; expected e.g. A:102-122")
    name = parts[2] if len(parts) > 2 else f"{chain}:{lo}-{hi}"
    return chain, lo, hi, name


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Inter-chain interface contacts, aggregated over models.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("structures", nargs="+", help="PDB/mmCIF files (globs ok)")
    ap.add_argument("--chains", nargs=2, metavar=("A", "B"),
                    help="chain pair to analyse (default: the two largest)")
    ap.add_argument("--cutoff", type=float, default=5.0,
                    help="heavy-atom contact cutoff in Angstrom (default 5.0)")
    ap.add_argument("--domain", action="append", default=[], metavar="CHAIN:LO-HI:NAME",
                    help="domain to break the interface down by; repeatable")
    ap.add_argument("--roi", action="append", default=[], metavar="CHAIN:LO-HI:NAME",
                    help="region of interest for the fraction readout; repeatable")
    ap.add_argument("--min-plddt", type=float, default=None,
                    help="drop residues whose mean B-factor/pLDDT is below this")
    ap.add_argument("--top", type=int, default=25, help="residues to list per chain")
    ap.add_argument("--tsv", help="write per-residue contact frequencies here")
    args = ap.parse_args()

    paths = []
    for pattern in args.structures:
        hits = sorted(glob.glob(pattern))
        paths.extend(hits if hits else [pattern])
    if not paths:
        sys.exit("no input structures")

    domains = [parse_region(d, "domain") for d in args.domain]
    rois = [parse_region(r, "roi") for r in args.roi]

    res_freq = defaultdict(int)  # (chain,resseq,icode,resname) -> models in contact
    pair_freq = defaultdict(int)
    n_ok = 0
    chain_a = chain_b = None

    for path in paths:
        try:
            atoms = load(path)
        except Exception as exc:  # noqa: BLE001 - report and carry on
            print(f"  ! skipping {path}: {exc}", file=sys.stderr)
            continue

        if args.min_plddt is not None:
            keep = np.ones(len(atoms), dtype=bool)
            by_res = defaultdict(list)
            for k in range(len(atoms)):
                by_res[(atoms.chain[k], int(atoms.resseq[k]), atoms.icode[k])].append(k)
            for ks in by_res.values():
                if float(np.mean(atoms.bfac[ks])) < args.min_plddt:
                    keep[ks] = False
            atoms = atoms.mask(keep)
            if len(atoms) == 0:
                print(f"  ! skipping {path}: nothing left after --min-plddt", file=sys.stderr)
                continue

        present = list(dict.fromkeys(atoms.chain.tolist()))
        if args.chains:
            ca, cb = args.chains
            if ca not in present or cb not in present:
                print(f"  ! skipping {path}: chains {ca}/{cb} not in {present}", file=sys.stderr)
                continue
        else:
            if len(present) < 2:
                print(f"  ! skipping {path}: only one chain ({present})", file=sys.stderr)
                continue
            sizes = sorted(present, key=lambda c: -int((atoms.chain == c).sum()))
            ca, cb = sizes[0], sizes[1]

        if chain_a is None:
            chain_a, chain_b = ca, cb
        elif (ca, cb) != (chain_a, chain_b):
            print(f"  ! {path}: chain pair {ca}/{cb} differs from {chain_a}/{chain_b}",
                  file=sys.stderr)

        pairs = contacts(atoms.mask(atoms.chain == ca), atoms.mask(atoms.chain == cb),
                         args.cutoff)
        n_ok += 1
        seen = set()
        for ra, rb in pairs:
            pair_freq[(ra, rb)] += 1
            for r in (ra, rb):
                if r not in seen:
                    seen.add(r)
                    res_freq[r] += 1

    if n_ok == 0:
        sys.exit("no structures could be analysed")

    print(f"\n{'=' * 72}")
    print(f"models analysed : {n_ok}")
    print(f"chains          : {chain_a} vs {chain_b}")
    print(f"cutoff          : {args.cutoff} A (heavy atoms)")
    if args.min_plddt is not None:
        print(f"pLDDT filter    : >= {args.min_plddt}")
    print(f"{'=' * 72}")

    for ch in (chain_a, chain_b):
        rows = sorted(((r, c) for r, c in res_freq.items() if r[0] == ch),
                      key=lambda x: (-x[1], x[0][1]))
        total = len(rows)
        print(f"\n--- chain {ch}: {total} interface residues over {n_ok} model(s) ---")
        if not rows:
            print("    (none)")
            continue
        print(f"    {'residue':>14}   freq   in models")
        for r, c in rows[: args.top]:
            label = f"{r[3]}{r[1]}{r[2]}"
            print(f"    {label:>14}   {c / n_ok:5.2f}   {c}/{n_ok}")
        if total > args.top:
            print(f"    ... {total - args.top} more")

        dom = [d for d in domains if d[0] == ch]
        if dom:
            print(f"\n    domain breakdown (by weighted contact frequency):")
            weighted = {r: c for r, c in rows}
            tot_w = sum(weighted.values()) or 1
            assigned = 0
            for _, lo, hi, name in dom:
                w = sum(c for r, c in weighted.items() if lo <= r[1] <= hi)
                assigned += w
                print(f"      {name:<20} {w / tot_w * 100:5.1f}%  "
                      f"({sum(1 for r in weighted if lo <= r[1] <= hi)} residues)")
            if assigned < tot_w:
                print(f"      {'(unassigned)':<20} {(tot_w - assigned) / tot_w * 100:5.1f}%")

        for _, lo, hi, name in (r for r in rois if r[0] == ch):
            weighted = {r: c for r, c in rows}
            tot_w = sum(weighted.values()) or 1
            w = sum(c for r, c in weighted.items() if lo <= r[1] <= hi)
            hit = sorted(r[1] for r in weighted if lo <= r[1] <= hi)
            print(f"\n    >>> ROI '{name}' ({ch}:{lo}-{hi}): "
                  f"{w / tot_w * 100:.1f}% of interface weight, "
                  f"{len(hit)}/{hi - lo + 1} residues contacted")
            if hit:
                print(f"        residues: {', '.join(map(str, hit))}")

    top_pairs = sorted(pair_freq.items(), key=lambda x: -x[1])[: args.top]
    if top_pairs:
        print(f"\n--- top residue-residue pairs ---")
        for (ra, rb), c in top_pairs:
            print(f"    {ra[3]}{ra[1]}{ra[2]} ({ra[0]}) -- {rb[3]}{rb[1]}{rb[2]} ({rb[0]})"
                  f"   {c / n_ok:5.2f}  {c}/{n_ok}")

    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("chain\tresseq\ticode\tresname\tn_models\tfrequency\n")
            for r, c in sorted(res_freq.items(), key=lambda x: (x[0][0], x[0][1])):
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{c}\t{c / n_ok:.4f}\n")
        print(f"\nwrote {args.tsv}")
    print()


if __name__ == "__main__":
    main()
