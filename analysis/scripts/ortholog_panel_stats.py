#!/usr/bin/env python3
"""
ortholog_panel_stats.py -- does one region of a protein diverge faster than another,
across a panel of orthologs, once the panel's own phylogeny is accounted for?

The problem this solves: docs/10 section 2 compared per-region sequence identity
between aRES modules with a two-proportion z-test on alignment positions, which
assumes (a) positions are independent observations and (b) orthologs are
independent of each other. Neither is true. Position counts inflate n enormously
-- 71 positions in NTD-1 is not 71 independent draws -- and orthologs share
ancestry, so two close relatives agreeing about NTD-1 is close to one observation,
not two. The uncorrected p = 0.037 in docs/10 rested on effectively one ortholog.

So this script moves the unit of analysis from the position to the ortholog, and
then corrects for the orthologs not being independent either:

  1. Per ortholog, the paired difference  d = identity(region A) - identity(region B).
     One number per ortholog per region pair; the position count drops out of the
     test and survives only as a noise estimate that is reported alongside.

  2. Three tests of  median(d) = 0 / mean(d) = 0:
       * exact Wilcoxon signed-rank (enumerated, no normal approximation)
       * exact sign test (binomial)
       * paired t-test
     The first two are distribution-free, which matters because identity
     differences over a 48-residue region are not plausibly normal.

  3. Holm-Bonferroni across the region pairs being compared, because comparing
     three regions means three tests and the table in docs/10 quoted the best one.

  4. Phylogenetic generalised least squares (PGLS). A UPGMA tree is built from the
     panel's own all-against-all identities, turned into a Brownian-motion
     covariance matrix V (shared root-to-MRCA path length), and the intercept-only
     model  d = beta + e,  e ~ N(0, sigma^2 V)  is fitted by generalised least
     squares. Testing beta = 0 asks whether the ancestral difference is nonzero
     rather than whether 18 correlated observations happen to be negative. This is
     the "phylogenetically aware comparison" docs/10 section 8 asked for. It is a
     correction for shared ancestry, not a substitute for a real gene tree: the
     tree here is distance-based, from the same alignments being tested, and the
     trait is identity *to the reference*, which is not a per-species character in
     the strict comparative sense. It is reported as a check on the paired tests,
     and where the two disagree, say so.

  5. Redundancy checks: single-linkage clustering of the panel at a range of
     identity thresholds, re-testing on cluster means, plus leave-one-out and
     leave-one-group-out jackknives. These answer "is this driven by a handful of
     close relatives", which is the specific failure mode of docs/10 section 2.

No scipy: the exact Wilcoxon distribution is enumerated by dynamic programming
over rank sums, the binomial is summed directly, and the Student-t tail comes from
a continued-fraction incomplete beta.

Dependencies: numpy.

Usage
-----
  python3 ortholog_panel_stats.py results/panel18_region_identity.tsv \
      --panel data/sequences/pf18870_res_panel18.fasta \
      --focal NTD1 --against Zn-ribbon RES \
      --groups results/panel18_groups.tsv \
      --out-prefix results/panel18
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homolog_compare import align, read_fasta_records  # noqa: E402


# ------------------------------------------------------------------ small stats


def norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (modified Lentz)."""
    TINY, EPS, ITMAX = 1e-30, 3e-16, 300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < TINY:
        d = TINY
    d = 1.0 / d
    h = d
    for m in range(1, ITMAX + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < TINY:
            d = TINY
        c = 1.0 + aa / c
        if abs(c) < TINY:
            c = TINY
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < TINY:
            d = TINY
        c = 1.0 + aa / c
        if abs(c) < TINY:
            c = TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h


def betai(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_two_sided(t: float, df: int) -> float:
    """Two-sided p-value for Student's t."""
    if df <= 0:
        return float("nan")
    if not math.isfinite(t):
        return 0.0
    return betai(0.5 * df, 0.5, df / (df + t * t))


def paired_t(d: np.ndarray):
    n = len(d)
    if n < 2:
        return float("nan"), float("nan"), n - 1
    sd = d.std(ddof=1)
    if sd == 0:
        return (float("inf") if d.mean() else 0.0), (0.0 if d.mean() else 1.0), n - 1
    t = d.mean() / (sd / math.sqrt(n))
    return t, t_two_sided(t, n - 1), n - 1


def sign_test(d: np.ndarray):
    """Exact two-sided binomial test on the signs, zeros discarded."""
    nz = d[d != 0]
    n = len(nz)
    k = int((nz < 0).sum())
    if n == 0:
        return 0, 0, float("nan")
    # two-sided: sum of all outcomes no more likely than the observed one
    probs = np.array([math.comb(n, i) for i in range(n + 1)], dtype=float) / 2.0 ** n
    p = float(probs[probs <= probs[k] + 1e-15].sum())
    return k, n, min(p, 1.0)


def wilcoxon_signed_rank(d: np.ndarray):
    """Exact two-sided Wilcoxon signed-rank test. -> (W_minus, n, p).

    Zeros are dropped (Wilcoxon's own handling). Ranks are midranks, so with ties
    the enumerated null is the distribution of sums of the *observed* rank
    multiset, which is the conditional exact test rather than an approximation.
    """
    nz = d[d != 0]
    n = len(nz)
    if n == 0:
        return float("nan"), 0, float("nan"), 0.0
    order = np.argsort(np.abs(nz), kind="mergesort")
    a = np.abs(nz)[order]
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and a[j + 1] == a[i]:
            j += 1
        ranks[i:j + 1] = 0.5 * (i + 1 + j + 1)
        i = j + 1
    signs = np.sign(nz)[order]
    w_minus = float(ranks[signs < 0].sum())
    w_plus = float(ranks[signs > 0].sum())
    total = float(ranks.sum())

    # Null: every sign pattern equally likely. Enumerate the distribution of the
    # negative-rank sum by DP over doubled ranks (midranks are half-integers).
    scaled = np.rint(ranks * 2).astype(int)
    size = int(scaled.sum()) + 1
    dist = np.zeros(size)
    dist[0] = 1.0
    for r in scaled:
        nxt = np.zeros(size)
        nxt[:size - r] += dist[:size - r]      # rank assigned +
        nxt[r:] += dist[:size - r]             # rank assigned -
        dist = nxt
    dist /= dist.sum()
    stat = min(w_minus, w_plus)
    idx = int(round(stat * 2))
    p = float(dist[:idx + 1].sum() + dist[size - 1 - idx:].sum())
    if 2 * idx == size - 1:                     # perfectly symmetric point
        p = 1.0
    return w_minus, n, min(p, 1.0), total


def holm(pvals, labels):
    """Holm-Bonferroni. -> {label: (p_raw, p_adjusted)}."""
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    m = len(pvals)
    adj, running = {}, 0.0
    for rank, i in enumerate(order):
        a = min(1.0, (m - rank) * pvals[i])
        running = max(running, a)               # enforce monotonicity
        adj[labels[i]] = (pvals[i], running)
    return adj


# ------------------------------------------------------- distances and the tree


def pairwise_identity(seqs):
    """All-against-all percent identity over aligned columns. -> (names, matrix)."""
    names = [n for n, _ in seqs]
    n = len(seqs)
    ident = np.full((n, n), 100.0)
    for i in range(n):
        for j in range(i + 1, n):
            pairs = align(seqs[i][1], seqs[j][1])
            eq = [(a, b) for a, b in pairs if a is not None and b is not None]
            if not eq:
                ident[i, j] = ident[j, i] = 0.0
                continue
            same = sum(1 for a, b in eq if seqs[i][1][a] == seqs[j][1][b])
            ident[i, j] = ident[j, i] = 100.0 * same / len(eq)
    return names, ident


def upgma(dist):
    """UPGMA on a square distance matrix. -> cophenetic merge-distance matrix.

    Only the cophenetic distances are needed downstream, so no tree object is
    built: coph[i, j] is the distance at which i and j were first joined, and
    half of it is the height of their MRCA above the leaves.
    """
    n = dist.shape[0]
    d = dist.astype(float).copy()
    np.fill_diagonal(d, np.inf)
    members = [[i] for i in range(n)]
    sizes = np.ones(n)
    alive = np.ones(n, dtype=bool)
    coph = np.zeros((n, n))
    for _ in range(n - 1):
        sub = np.where(alive)[0]
        m = d[np.ix_(sub, sub)]
        k = int(np.argmin(m))
        a, b = sub[k // len(sub)], sub[k % len(sub)]
        h = d[a, b]
        for x in members[a]:
            for y in members[b]:
                coph[x, y] = coph[y, x] = h
        # average linkage, weighted by cluster size
        for c in np.where(alive)[0]:
            if c in (a, b):
                continue
            d[a, c] = d[c, a] = (
                (sizes[a] * d[a, c] + sizes[b] * d[b, c]) / (sizes[a] + sizes[b])
            )
        members[a] = members[a] + members[b]
        sizes[a] += sizes[b]
        alive[b] = False
        d[b, :] = d[:, b] = np.inf
        d[a, a] = np.inf
    return coph


def brownian_covariance(coph):
    """Ultrametric cophenetic distances -> Brownian-motion covariance.

    Under Brownian motion the covariance of two tips is the branch length they
    share, i.e. the root-to-MRCA path. With MRCA heights h = coph / 2 and root
    height H, that is V[i, j] = H - h[i, j], and V[i, i] = H.
    """
    h = coph / 2.0
    H = float(h.max()) if h.size else 1.0
    if H <= 0:
        H = 1.0
    V = H - h
    np.fill_diagonal(V, H)
    # Numerically nudge onto the PSD cone if the linkage left it marginal.
    w = np.linalg.eigvalsh(V)
    if w.min() <= 1e-9 * H:
        V = V + np.eye(V.shape[0]) * (1e-6 * H - min(0.0, w.min()))
    return V / H          # scale-free; sigma^2 absorbs the units


def pgls_intercept(y, V):
    """GLS fit of y = beta + e, e ~ N(0, sigma^2 V). -> (beta, se, t, p, df).

    Testing beta = 0 under a Brownian covariance is the phylogenetic analogue of
    a one-sample t-test: it asks whether the trait's root value differs from zero,
    weighting tips by how much independent information each contributes, instead
    of counting correlated tips as independent observations.
    """
    n = len(y)
    one = np.ones((n, 1))
    Vi = np.linalg.inv(V)
    xtx = (one.T @ Vi @ one).item()
    beta = (one.T @ Vi @ y.reshape(-1, 1)).item() / xtx
    resid = y - beta
    df = n - 1
    sigma2 = (resid.reshape(1, -1) @ Vi @ resid.reshape(-1, 1)).item() / df
    se = math.sqrt(max(sigma2, 0.0) / xtx)
    if se == 0:
        return beta, se, float("nan"), float("nan"), df
    t = beta / se
    return beta, se, t, t_two_sided(t, df), df


def single_linkage_clusters(ident, threshold):
    """Group sequences joined by a chain of >= threshold identity."""
    n = ident.shape[0]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if ident[i, j] >= threshold:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


# ------------------------------------------------------------------------ input


def read_region_tsv(path, metric, unit_col="ortholog", region_col="region"):
    """-> ({unit: {region: value}}, {unit: {region: n_observations}}).

    Shaped for homolog_compare.py's --identity-tsv, but any long-format table of
    one value per unit per region works -- charge_composition.py's --tsv is read
    the same way, so a composition claim gets the same corrected tests as an
    identity claim instead of a separate ad hoc calculation.
    """
    vals, counts = {}, {}
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for name in (unit_col, region_col, metric):
            if name not in cols:
                raise SystemExit(f"{path}: no column {name!r} (have {cols})")
        n_col = next((c for c in ("n_aligned", "len") if c in cols), None)
        for line in fh:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            if not r.get(metric, "").strip():
                continue
            o, reg = r[unit_col], r[region_col]
            vals.setdefault(o, {})[reg] = float(r[metric])
            counts.setdefault(o, {})[reg] = int(r[n_col]) if n_col and r[n_col] else 0
    # A unit missing one of the regions would silently shorten a paired vector.
    regs = {reg for v in vals.values() for reg in v}
    for o, v in vals.items():
        if set(v) != regs:
            raise SystemExit(f"{path}: {o} has {sorted(v)}, expected {sorted(regs)}")
    return vals, counts


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Ortholog-level, phylogeny-corrected test of per-region divergence.")
    ap.add_argument("identity_tsv", help="--identity-tsv output of homolog_compare.py, "
                                         "or any long-format per-region value table")
    ap.add_argument("--panel", help="panel FASTA, for the all-against-all tree")
    ap.add_argument("--focal", default="NTD1", help="region under test (default NTD1)")
    ap.add_argument("--against", nargs="*", default=[],
                    help="regions to compare the focal region against "
                         "(default: every other region in the TSV)")
    ap.add_argument("--metric", default="id_aligned",
                    help="value column to test (default id_aligned; id_allpos for "
                         "the indel-inclusive identity, net_per100 for a charge TSV)")
    ap.add_argument("--unit-col", default="ortholog",
                    help="column naming the ortholog (default ortholog; 'record' in "
                         "a charge_composition.py TSV)")
    ap.add_argument("--region-col", default="region")
    ap.add_argument("--groups", help="TSV with columns ortholog<TAB>group, for the "
                                     "leave-one-group-out jackknife")
    ap.add_argument("--thresholds", nargs="*", type=float,
                    default=[40.0, 50.0, 60.0, 70.0],
                    help="identity thresholds for the redundancy check")
    ap.add_argument("--out-prefix", help="write TSVs with this path prefix")
    ap.add_argument("--dist-cache", help="cache the all-against-all identity matrix here")
    args = ap.parse_args()

    vals, counts = read_region_tsv(args.identity_tsv, args.metric,
                                   args.unit_col, args.region_col)
    orthologs = sorted(vals)
    regions = [r for r in dict.fromkeys(
        reg for o in orthologs for reg in vals[o]) if r != "overall"]
    if args.focal not in regions:
        raise SystemExit(f"focal region {args.focal!r} not in {regions}")
    against = args.against or [r for r in regions if r != args.focal]
    for r in against:
        if r not in regions:
            raise SystemExit(f"region {r!r} not in {regions}")

    print("=" * 78)
    print(f"per-region comparison, ortholog-level tests   ({args.metric})")
    print("=" * 78)
    print(f"  panel          {len(orthologs)} orthologs")
    print(f"  focal region   {args.focal}")
    print(f"  compared with  {', '.join(against)}")

    # ---- the descriptive table, which is what the tests are about
    print(f"\nmean {args.metric} by region, over the panel")
    print(f"  {'region':<14}{'mean':>8}{'sd':>8}{'min':>8}{'max':>8}"
          f"{'median n':>13}")
    for r in regions:
        v = np.array([vals[o][r] for o in orthologs])
        n_al = np.array([counts[o][r] for o in orthologs])
        print(f"  {r:<14}{v.mean():>8.1f}{v.std(ddof=1):>8.1f}{v.min():>8.1f}"
              f"{v.max():>8.1f}{np.median(n_al):>13.0f}")
    print("  (n is the residues or alignment columns each ortholog's own value rests "
          "on --\n   a 48-residue region is a far noisier estimate than a "
          "322-residue one)")

    # ---- paired, ortholog-level tests
    print(f"\npaired differences  d = {args.metric}({args.focal}) - "
          f"{args.metric}(other), one per ortholog")
    print(f"  {'comparison':<26}{'mean d':>9}{'median d':>10}{'neg/n':>9}"
          f"{'Wilcoxon p':>12}{'sign p':>9}{'t p':>9}")
    raw = {}
    diffs = {}
    for other in against:
        d = np.array([vals[o][args.focal] - vals[o][other] for o in orthologs])
        diffs[other] = d
        w_minus, nw, p_w, _tot = wilcoxon_signed_rank(d)
        k, ns, p_s = sign_test(d)
        t, p_t, _df = paired_t(d)
        lab = f"{args.focal} - {other}"
        raw[lab] = dict(other=other, d=d, p_w=p_w, p_s=p_s, p_t=p_t, t=t,
                        k=k, ns=ns)
        print(f"  {lab:<26}{d.mean():>+9.2f}{np.median(d):>+10.2f}"
              f"{f'{k}/{ns}':>9}{p_w:>12.4f}{p_s:>9.4f}{p_t:>9.4f}")
    print("  (neg/n = orthologs in which the focal region has the lower value)")

    # ---- multiple-comparison correction
    labs = list(raw)
    print(f"\nHolm-Bonferroni across the {len(labs)} region pair(s) tested "
          f"(Wilcoxon p adjusted)")
    adj = holm([raw[l]["p_w"] for l in labs], labs)
    for l in labs:
        p0, p1 = adj[l]
        print(f"  {l:<26}p = {p0:.4f}  ->  p_adj = {p1:.4f}"
              f"{'   *' if p1 < 0.05 else ''}")

    # ---- phylogenetic correction
    pgls_rows = []
    ident = None
    names = None
    if args.panel:
        seqs = read_fasta_records(args.panel)
        seqs = [(n, s) for n, s in seqs if n in set(orthologs)]
        if len(seqs) != len(orthologs):
            missing = set(orthologs) - {n for n, _ in seqs}
            raise SystemExit(f"panel FASTA is missing {sorted(missing)}")
        cache = args.dist_cache
        if cache and os.path.exists(cache):
            arr = np.loadtxt(cache, delimiter="\t", skiprows=1, dtype=str)
            names = list(arr[:, 0])
            ident = arr[:, 1:].astype(float)
            print(f"\nall-against-all identity read from {cache}")
        else:
            print(f"\naligning {len(seqs)}x{len(seqs)} for the panel tree "
                  f"({len(seqs) * (len(seqs) - 1) // 2} pairwise alignments)")
            names, ident = pairwise_identity(seqs)
            if cache:
                os.makedirs(os.path.dirname(os.path.abspath(cache)), exist_ok=True)
                with open(cache, "w", newline="") as fh:
                    fh.write("name\t" + "\t".join(names) + "\n")
                    for i, nm in enumerate(names):
                        fh.write(nm + "\t" + "\t".join(f"{x:.4f}" for x in ident[i]) + "\n")
        order = [names.index(o) for o in orthologs]
        ident = ident[np.ix_(order, order)]
        names = orthologs

        off = ident[~np.eye(len(names), dtype=bool)]
        print(f"  panel internal identity: {off.min():.1f}-{off.max():.1f}% "
              f"(mean {off.mean():.1f}%)")

        for dist_label, D in (
            ("p-distance", 1.0 - ident / 100.0),
            ("Poisson-corrected", -np.log(np.clip(ident / 100.0, 1e-6, 1.0))),
        ):
            np.fill_diagonal(D, 0.0)
            V = brownian_covariance(upgma(D))
            print(f"\nPGLS, intercept only, Brownian covariance from a UPGMA tree "
                  f"({dist_label})")
            print(f"  {'comparison':<26}{'beta':>9}{'se':>8}{'t':>8}{'df':>5}"
                  f"{'p':>10}")
            ps = []
            for other in against:
                lab = f"{args.focal} - {other}"
                beta, se, t, p, df = pgls_intercept(diffs[other], V)
                ps.append(p)
                pgls_rows.append(dict(dist=dist_label, comparison=lab, beta=beta,
                                      se=se, t=t, df=df, p=p))
                print(f"  {lab:<26}{beta:>+9.2f}{se:>8.2f}{t:>8.2f}{df:>5}"
                      f"{p:>10.4f}")
            a2 = holm(ps, [f"{args.focal} - {o}" for o in against])
            for o in against:
                lab = f"{args.focal} - {o}"
                print(f"    Holm-adjusted  {lab:<24}p_adj = {a2[lab][1]:.4f}"
                      f"{'   *' if a2[lab][1] < 0.05 else ''}")
                for row in pgls_rows:
                    if row["dist"] == dist_label and row["comparison"] == lab:
                        row["p_adj"] = a2[lab][1]

        # ---- redundancy: are a few close relatives carrying the result?
        print(f"\nredundancy check: collapse sequences joined at >= T% identity, "
              f"re-test on cluster means")
        print(f"  {'T':>5}{'clusters':>10}  " +
              "".join(f"{args.focal + ' - ' + o:>26}" for o in against))
        for T in args.thresholds:
            cl = single_linkage_clusters(ident, T)
            cells = []
            for other in against:
                d = diffs[other]
                dm = np.array([d[c].mean() for c in cl])
                if len(dm) < 3:
                    cells.append("n<3")
                    continue
                _w, _n, p_w, _t = wilcoxon_signed_rank(dm)
                cells.append(f"mean {dm.mean():+.2f}, p = {p_w:.4f}")
            print(f"  {T:>5.0f}{len(cl):>10}  " +
                  "".join(f"{c:>26}" for c in cells))

    # ---- jackknives
    print(f"\nleave-one-ortholog-out jackknife (Wilcoxon p over the "
          f"{len(orthologs) - 1} remaining)")
    print(f"  {'comparison':<26}{'p min':>9}{'p max':>9}{'worst drop':>26}")
    for other in against:
        d = diffs[other]
        res = []
        for k in range(len(d)):
            dd = np.delete(d, k)
            _w, _n, p, _t = wilcoxon_signed_rank(dd)
            res.append((p, orthologs[k]))
        res.sort()
        print(f"  {args.focal + ' - ' + other:<26}{res[0][0]:>9.4f}"
              f"{res[-1][0]:>9.4f}{res[-1][1]:>26}")

    groups = {}
    if args.groups:
        with open(args.groups) as fh:
            first = fh.readline()
            if "ortholog" not in first.lower():
                fh.seek(0)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    groups[f[0]] = f[1]
    if groups:
        gs = sorted(set(groups.get(o, "?") for o in orthologs))
        print(f"\nleave-one-group-out jackknife ({len(gs)} groups: {', '.join(gs)})")
        print(f"  {'dropped group':<22}{'n left':>8}  " +
              "".join(f"{args.focal + ' - ' + o:>26}" for o in against))
        for g in gs:
            keep = [i for i, o in enumerate(orthologs) if groups.get(o, "?") != g]
            if len(keep) < 4:
                continue
            cells = []
            for other in against:
                dd = diffs[other][keep]
                _w, _n, p, _t = wilcoxon_signed_rank(dd)
                cells.append(f"mean {dd.mean():+.2f}, p = {p:.4f}")
            print(f"  {g:<22}{len(keep):>8}  " + "".join(f"{c:>26}" for c in cells))

    # ---- output
    if args.out_prefix:
        d = os.path.dirname(os.path.abspath(args.out_prefix))
        os.makedirs(d, exist_ok=True)

        p1 = f"{args.out_prefix}_paired_diffs.tsv"
        with open(p1, "w", newline="") as fh:
            fh.write("ortholog\t" + "\t".join(f"id_{r}" for r in regions) + "\t"
                     + "\t".join(f"d_{args.focal}_minus_{o}" for o in against) + "\n")
            for k, o in enumerate(orthologs):
                fh.write(o + "\t"
                         + "\t".join(f"{vals[o][r]:.4f}" for r in regions) + "\t"
                         + "\t".join(f"{diffs[x][k]:.4f}" for x in against) + "\n")

        p2 = f"{args.out_prefix}_tests.tsv"
        with open(p2, "w", newline="") as fh:
            fh.write("test\tdistance\tcomparison\tstatistic\tn\tp\tp_holm\n")
            adjw = holm([raw[l]["p_w"] for l in labs], labs)
            for l in labs:
                r = raw[l]
                fh.write(f"wilcoxon_signed_rank_exact\t-\t{l}\t-\t"
                         f"{len(r['d'])}\t{r['p_w']:.6g}\t{adjw[l][1]:.6g}\n")
                fh.write(f"sign_test_exact\t-\t{l}\t{r['k']}/{r['ns']}\t"
                         f"{r['ns']}\t{r['p_s']:.6g}\t\n")
                fh.write(f"paired_t\t-\t{l}\t{r['t']:.4f}\t{len(r['d'])}\t"
                         f"{r['p_t']:.6g}\t\n")
            for row in pgls_rows:
                fh.write(f"pgls_intercept\t{row['dist']}\t{row['comparison']}\t"
                         f"t={row['t']:.4f} beta={row['beta']:.4f}\t"
                         f"{row['df'] + 1}\t{row['p']:.6g}\t"
                         f"{row.get('p_adj', float('nan')):.6g}\n")
        outs = [p1, p2]
        if ident is not None:
            p3 = f"{args.out_prefix}_panel_identity_matrix.tsv"
            with open(p3, "w", newline="") as fh:
                fh.write("name\t" + "\t".join(names) + "\n")
                for i, nm in enumerate(names):
                    fh.write(nm + "\t" + "\t".join(f"{x:.2f}" for x in ident[i]) + "\n")
            outs.append(p3)
        print("\nwrote " + ", ".join(outs))
    print()


if __name__ == "__main__":
    main()
