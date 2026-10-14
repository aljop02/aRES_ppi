#!/usr/bin/env python3
"""
af3_zn_site.py -- measure a modelled Zn site across an AlphaFold 3 job's models.

Companion to zn_site_check.py. That script infers a zinc site from cysteine
geometry alone, because the models it was written for carried no metal. This one
is for AF3 jobs where Zn *was* supplied as a ligand: the metal is present, so the
coordination shell can be measured directly instead of inferred.

What it reports, per model and aggregated over the job:

  * the coordination shell: every protein heavy atom within --cutoff of the Zn
  * Zn-ligand bond lengths against the 2.3 A reference for Zn-S(Cys)
  * L-Zn-L angles and their RMS deviation from tetrahedral (109.47 deg)
  * pLDDT of the coordinating residues and of the Zn itself
  * the job's confidence metrics from the summary_confidences JSONs

Two cautions on reading the output:

  ipTM with a single ion is not an interface score in the usual sense. It says the
  model places the ion confidently relative to the chain; it does not validate the
  chain. Read chain_ptm (the protein's own pTM) separately -- a job can show
  ipTM 0.97 and chain_ptm 0.61 at the same time, and both numbers are honest.

  AF3 will place a supplied ion somewhere regardless. A site is only credible if
  the geometry is right (tetrahedral, 2.3 A) AND the ligands are confident AND the
  same shell recurs across models. All three are reported so none is assumed.

Dependencies: numpy.

Usage
-----
  python3 af3_zn_site.py <dir-with-model-cifs> [...] \
      --regions 1-48:Zn-ribbon 49-120:NTD1 126-447:RES 448-459:tail \
      --residues 89 93 104 \
      --tsv results/zn_sites.tsv
"""

from __future__ import annotations

import argparse
import glob
import itertools
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

ZN_S_REF = 2.30  # Zn-S(Cys) bond length, A
ZN_N_REF = 2.05  # Zn-N(His) bond length, A
ZN_O_REF = 2.00  # Zn-O(Asp/Glu) bond length, A
TETRAHEDRAL = 109.4712

# Side-chain atoms that actually coordinate zinc. Restricting the shell to these
# matters for the angle statistics: a Cys CB sits ~2.8 A from the metal purely
# because its own SG is bonded to it, and counting CB as a fifth ligand inflates
# the tetrahedral RMSD from ~5 deg to ~28 deg on an otherwise perfect site.
# Non-donor atoms inside the cutoff are reported separately, not scored.
DONOR_ATOMS = {
    ("CYS", "SG"), ("HIS", "ND1"), ("HIS", "NE2"),
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
    ("MET", "SD"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"),
    ("HOH", "O"),
}


def ref_bond(element: str) -> float:
    return {"S": ZN_S_REF, "N": ZN_N_REF, "O": ZN_O_REF}.get(element, ZN_O_REF)


# ---------------------------------------------------------------- CIF parsing


def parse_cif_atoms(path: str):
    """Minimal mmCIF atom_site reader -> list of dicts. Handles AF3 server output.

    Reads the loop header to find column positions rather than assuming an order,
    so it survives AF3 changing the column layout.
    """
    cols: list[str] = []
    rows = []
    in_loop = False
    in_atoms = False
    with open(path, errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("loop_"):
                in_loop, in_atoms, cols = True, False, []
                continue
            if in_loop and s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1].split()[0])
                in_atoms = True
                continue
            if in_atoms:
                if s.startswith("#") or s.startswith("loop_") or s.startswith("_") or not s:
                    if s.startswith("loop_"):
                        in_loop, in_atoms, cols = True, False, []
                        continue
                    in_atoms = False
                    continue
                parts = shlex_split_cif(s)
                if len(parts) != len(cols):
                    continue
                rows.append(dict(zip(cols, parts)))
    if not rows:
        raise SystemExit(f"{path}: no _atom_site rows parsed")

    atoms = []
    for r in rows:
        if r.get("pdbx_PDB_model_num", "1") not in ("1", "."):
            continue
        try:
            xyz = np.array([float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])])
        except (KeyError, ValueError):
            continue
        try:
            resseq = int(r.get("auth_seq_id", r.get("label_seq_id", "0")))
        except ValueError:
            resseq = 0
        try:
            plddt = float(r.get("B_iso_or_equiv", "nan"))
        except ValueError:
            plddt = float("nan")
        atoms.append(
            {
                "group": r.get("group_PDB", "ATOM"),
                "element": r.get("type_symbol", "").strip().upper(),
                "atom": r.get("label_atom_id", "").strip().strip('"'),
                "resname": r.get("label_comp_id", "").strip().upper(),
                "chain": r.get("auth_asym_id", r.get("label_asym_id", "A")).strip(),
                "resseq": resseq,
                "xyz": xyz,
                "plddt": plddt,
            }
        )
    return atoms


def shlex_split_cif(line: str) -> list[str]:
    """Split a CIF data line, respecting single and double quotes."""
    out, cur, quote = [], "", None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
            else:
                cur += ch
        elif ch in "'\"":
            quote = ch
        elif ch.isspace():
            if cur:
                out.append(cur)
                cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


# ---------------------------------------------------------------- measurement


def residue_plddt(atoms) -> dict:
    acc = defaultdict(list)
    for a in atoms:
        if not np.isnan(a["plddt"]):
            acc[(a["chain"], a["resseq"])].append(a["plddt"])
    return {k: float(np.mean(v)) for k, v in acc.items()}


def zn_shell(atoms, cutoff: float):
    """Coordination shell of the first Zn.

    -> (zn_atom, donors, others, n_zn) where donors are side-chain donor atoms
    within the cutoff (these define the geometry) and others are any remaining
    close atoms, reported for transparency but excluded from the angle stats.
    """
    zns = [a for a in atoms if a["element"] == "ZN"]
    if not zns:
        return None, [], [], 0
    zn = zns[0]
    donors, others = [], []
    for a in atoms:
        if a["element"] in ("ZN", "H") or a is zn:
            continue
        d = float(np.linalg.norm(a["xyz"] - zn["xyz"]))
        if d <= cutoff:
            bucket = donors if (a["resname"], a["atom"]) in DONOR_ATOMS else others
            bucket.append((d, a))
    key = lambda t: t[0]
    return zn, sorted(donors, key=key), sorted(others, key=key), len(zns)


def tetrahedral_stats(zn, shell):
    """RMS deviation of L-Zn-L angles from 109.47 deg, plus the angle list."""
    if len(shell) < 3:
        return float("nan"), []
    vecs = [(a["xyz"] - zn["xyz"]) for _, a in shell]
    angles = []
    for (i, u), (j, v) in itertools.combinations(list(enumerate(vecs)), 2):
        cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
        angles.append((i, j, np.degrees(np.arccos(np.clip(cos, -1, 1)))))
    rms = float(np.sqrt(np.mean([(a - TETRAHEDRAL) ** 2 for _, _, a in angles])))
    return rms, angles


def region_means(plddt: dict, regions: list[str]) -> list[tuple]:
    out = []
    for spec in regions:
        rng, _, label = spec.partition(":")
        label = label or rng
        try:
            lo, hi = (int(x) for x in rng.split("-"))
        except ValueError:
            print(f"  bad region '{spec}', expected LO-HI[:NAME]", file=sys.stderr)
            continue
        vals = [v for (ch, r), v in plddt.items() if lo <= r <= hi]
        out.append((label, float(np.mean(vals)) if vals else float("nan"), len(vals)))
    return out


# ---------------------------------------------------------------- per job


def load_summaries(d: str) -> list[dict]:
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*summary_confidences_*.json"))):
        with open(p) as fh:
            s = json.load(fh)
        s["_file"] = os.path.basename(p)
        out.append(s)
    return out


def analyse_job(d: str, args) -> dict:
    models = sorted(glob.glob(os.path.join(d, "*model_*.cif")))
    if not models:
        print(f"{d}: no *model_*.cif found", file=sys.stderr)
        return {}
    name = os.path.basename(os.path.normpath(d))

    print("=" * 78)
    print(f"{name}   ({len(models)} models)")
    print("=" * 78)

    summaries = load_summaries(d)
    if summaries:
        print("\n  confidence metrics per model")
        print(f"    {'model':<7}{'rank':>7}{'pTM':>7}{'ipTM':>7}{'chainPTM':>10}"
              f"{'fracDis':>9}{'clash':>7}{'PAEmin':>8}")
        for i, s in enumerate(summaries):
            cp = s.get("chain_ptm") or [None]
            chain_ptm = cp[0] if cp and cp[0] is not None else float("nan")
            pae = s.get("chain_pair_pae_min") or [[None, None]]
            pae_ab = pae[0][1] if len(pae[0]) > 1 and pae[0][1] is not None else float("nan")
            print(f"    {i:<7}{s.get('ranking_score', float('nan')):>7.3f}"
                  f"{s.get('ptm', float('nan')):>7.2f}{s.get('iptm', float('nan')):>7.2f}"
                  f"{chain_ptm:>10.2f}{s.get('fraction_disordered', float('nan')):>9.2f}"
                  f"{s.get('has_clash', float('nan')):>7.1f}{pae_ab:>8.2f}")

    rows = []
    shells_seen = defaultdict(int)
    for mi, path in enumerate(models):
        atoms = parse_cif_atoms(path)
        plddt = residue_plddt(atoms)
        zn, shell, others, n_zn = zn_shell(atoms, args.cutoff)
        print(f"\n  --- model {mi}  ({os.path.basename(path)})")
        if zn is None:
            print("      no Zn in the model")
            continue
        if n_zn > 1:
            print(f"      note: {n_zn} Zn atoms present; measuring the first")
        rms, _ = tetrahedral_stats(zn, shell)
        print(f"      Zn pLDDT {zn['plddt']:.1f}   "
              f"donors within {args.cutoff:.1f} A: {len(shell)}   "
              f"tetrahedral RMSD {rms:.1f} deg")
        if shell:
            print(f"      {'ligand':<16}{'atom':>6}{'dist':>8}{'ref':>7}{'dev':>8}{'pLDDT':>8}")
        key = []
        for d_, a in shell:
            r = ref_bond(a["element"])
            print(f"      {a['resname']}{a['resseq']:<12}{a['atom']:>6}"
                  f"{d_:>8.2f}{r:>7.2f}{d_ - r:>+8.2f}"
                  f"{plddt.get((a['chain'], a['resseq']), float('nan')):>8.1f}")
            key.append(f"{a['resname']}{a['resseq']}")
            rows.append({
                "job": name, "model": mi, "lig_res": f"{a['resname']}{a['resseq']}",
                "lig_atom": a["atom"], "dist": round(d_, 3),
                "ref": r, "dev": round(d_ - r, 3),
                "lig_plddt": round(plddt.get((a["chain"], a["resseq"]), float("nan")), 1),
                "zn_plddt": round(zn["plddt"], 1), "tet_rmsd_deg": round(rms, 2),
            })
        shells_seen["+".join(sorted(set(key), key=lambda s: int(re.sub(r"\D", "", s))))] += 1

        if others:
            print("      non-donor atoms inside the cutoff (not scored): " + ", ".join(
                f"{a['resname']}{a['resseq']}:{a['atom']} {d_:.2f}" for d_, a in others))

        if args.residues:
            picks = [(r, plddt.get(("A", r))) for r in args.residues]
            seen = [(r, v) for r, v in picks if v is not None]
            if seen:
                res_by_num = {}
                for a in atoms:
                    res_by_num.setdefault(a["resseq"], a["resname"])
                print("      residues of interest: " + ", ".join(
                    f"{res_by_num.get(r, '?')}{r} pLDDT {v:.1f}" for r, v in seen))

        if args.regions:
            print("      region pLDDT: " + "   ".join(
                f"{lab} {m:.1f}" for lab, m, n in region_means(plddt, args.regions)))

    if shells_seen:
        print(f"\n  coordination shell reproducibility across {len(models)} models")
        for shell_key, n in sorted(shells_seen.items(), key=lambda t: -t[1]):
            print(f"    {n}/{len(models)}  {shell_key}")

    if rows:
        s_rows = [r for r in rows if r["lig_atom"] == "SG"]
        if s_rows:
            ds = np.array([r["dist"] for r in s_rows])
            print(f"\n  Zn-SG bonds over all models: n={len(ds)}  "
                  f"mean {ds.mean():.2f} A  range {ds.min():.2f}-{ds.max():.2f} A  "
                  f"(reference {ZN_S_REF:.2f} A)")
        tet = np.array(sorted({(r["model"], r["tet_rmsd_deg"]) for r in rows}))
        if len(tet):
            print(f"  tetrahedral RMSD over models: mean {tet[:, 1].mean():.1f} deg  "
                  f"range {tet[:, 1].min():.1f}-{tet[:, 1].max():.1f} deg")
    print()
    return {"name": name, "rows": rows, "summaries": summaries}


def main() -> None:
    ap = argparse.ArgumentParser(description="Measure a modelled Zn site in AF3 output.")
    ap.add_argument("dirs", nargs="+", help="directories containing *model_*.cif")
    ap.add_argument("--cutoff", type=float, default=3.0,
                    help="Zn-ligand distance cutoff in A (default 3.0)")
    ap.add_argument("--regions", nargs="*", default=[], metavar="LO-HI:NAME")
    ap.add_argument("--residues", nargs="*", type=int, default=[],
                    help="residue numbers to report pLDDT for")
    ap.add_argument("--tsv", help="write per-ligand rows to this TSV")
    args = ap.parse_args()

    results = [analyse_job(d, args) for d in args.dirs]

    if args.tsv:
        rows = [r for res in results for r in res.get("rows", [])]
        if rows:
            cols = list(rows[0].keys())
            os.makedirs(os.path.dirname(os.path.abspath(args.tsv)), exist_ok=True)
            with open(args.tsv, "w", newline="") as fh:
                fh.write("\t".join(cols) + "\n")
                for r in rows:
                    fh.write("\t".join(str(r[c]) for c in cols) + "\n")
            print(f"wrote {args.tsv}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
