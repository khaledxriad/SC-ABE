#!/usr/bin/env python3
"""Reproduces every measured number in the associated article. Writes results/results.json.

Usage:  python experiments/run_experiments.py [--quick]
Deterministic apart from wall-clock timings (seed 2026).
"""
import argparse, json, math, os, random, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from scabe import IBE, SCABE, anc, cover, cover_set, leaves, node_id, prime_below  # noqa: E402

SEED = 2026
random.seed(SEED)
rng = np.random.default_rng(SEED)


def exhaustive(b, d):
    L = b ** d
    worst = 0
    for a in range(L):
        for c in range(a, L):
            P = cover(a, c, b, d)
            seen = []
            for nd in P:
                lo, hi = leaves(nd, b, d)
                seen += list(range(lo, hi + 1))
            assert sorted(seen) == list(range(a, c + 1)), (a, c)       # exact partition
            for (j, i) in P:                                            # maximality
                if j > 0:
                    lo, hi = leaves((j - 1, i // b), b, d)
                    assert not (a <= lo and hi <= c)
            for x in range(L):                                          # path lemma
                hit = [nd for nd in P if nd in set(anc(x, b, d))]
                assert (len(hit) == 1) == (a <= x <= c) and len(hit) <= 1
            worst = max(worst, len(P))
    return L, worst, 2 * (b - 1) * d


def rand_runs(r, L):
    v = sorted(random.sample(range(L - r), 2 * r))
    return [(v[2 * i] + i, v[2 * i + 1] + i) for i in range(r)]


def cover_stats(trials):
    res = {}
    for (b, d) in [(2, 20), (4, 10), (16, 5)]:
        L = b ** d
        rows = []
        for r in [1, 2, 4, 8, 16, 32, 64, 128]:
            sizes = [sum(len(cover(a, c, b, d)) for a, c in rand_runs(r, L)) for _ in range(trials)]
            rows.append(dict(r=r, mean=float(np.mean(sizes)), mx=int(max(sizes)), bound=2 * (b - 1) * d * r))
        res[f"b{b}"] = dict(b=b, d=d, L=L, rows=rows)
    return res


def single_interval_vs_L(trials):
    rows = []
    for d in range(6, 31, 2):
        L = 2 ** d
        sz = [len(cover(*rand_runs(1, L)[0], 2, d)) for _ in range(trials)]
        rows.append(dict(d=d, mean=float(np.mean(sz)), mx=int(max(sz)), bound=2 * d))
    return rows


def arity_tradeoff(trials):
    L = 2 ** 20
    rows = []
    for b in [2, 4, 8, 16, 32, 64, 128, 1024]:
        d = round(math.log(L, b))
        if b ** d != L:
            continue
        sz = [len(cover(*rand_runs(1, L)[0], b, d)) for _ in range(trials)]
        rows.append(dict(b=b, d=d, keynodes=d + 1, ctmean=float(np.mean(sz)), ctmax=int(max(sz)), bound=2 * (b - 1) * d))
    return rows


def revocation(trials):
    L, d = 2 ** 20, 20
    rows = []
    for t in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        sz = []
        for _ in range(trials):
            rev = sorted(random.sample(range(L), t))
            S, prev = [], 0
            for x in rev:
                if x > prev:
                    S.append((prev, x - 1))
                prev = x + 1
            if prev <= L - 1:
                S.append((prev, L - 1))
            sz.append(sum(len(cover(a, c, 2, d)) for a, c in S))
        rows.append(dict(t=t, mean=float(np.mean(sz)), nnl=t * math.log2(L / t), runs=t + 1, bound=2 * d * (t + 1)))
    return rows



def fl_cohorts(trials):
    """Cohort-addressing scenarios for a fleet of 2^20 devices with hierarchical identifiers (b=2, d=20)."""
    L, d = 2 ** 20, 20
    rows = []
    def add(name, members, gen):
        sizes = [gen() for _ in range(trials)]
        rows.append(dict(name=name, members=members, mean=float(np.mean(sizes)), mx=int(max(sizes))))
    def interval(m):
        def g():
            a = random.randrange(L - m + 1)
            return len(cover(a, a + m - 1, 2, d))
        return g
    def blocks(k, m):          # k random aligned blocks of m devices (whole sub-fleets such as sites)
        def g():
            starts = random.sample(range(L // m), k)
            return sum(len(cover(s * m, s * m + m - 1, 2, d)) for s in starts)
        return g
    def all_but(t):
        def g():
            rev = sorted(random.sample(range(L), t)); tot, prev = 0, 0
            for x in rev:
                if x > prev: tot += len(cover(prev, x - 1, 2, d))
                prev = x + 1
            if prev <= L - 1: tot += len(cover(prev, L - 1, 2, d))
            return tot
        return g
    def sample(m):             # unstructured random sample of m devices
        def g():
            return len(cover_set(random.sample(range(L), m), 2, d))
        return g
    add("one site (aligned block of 2^10 devices)", 2 ** 10, blocks(1, 2 ** 10))
    add("one region (arbitrary interval of 2^16 devices)", 2 ** 16, interval(2 ** 16))
    add("8 sites (aligned blocks of 2^10 devices)", 8 * 2 ** 10, blocks(8, 2 ** 10))
    add("whole fleet except 16 compromised devices", L - 16, all_but(16))
    add("whole fleet except 256 compromised devices", L - 256, all_but(256))
    add("random sample of 100 devices", 100, sample(100))
    add("random sample of 1000 devices", 1000, sample(1000))
    return rows

def toy_correctness():
    ibe = IBE(16, prime_below(65537), 128, seed=SEED)
    ibe.setup()
    ok = bad = tests = false_pos = 0
    for (b, d) in [(2, 6), (3, 3), (4, 3)]:
        sc = SCABE(b, d, ibe, seed=SEED + b)
        L = b ** d
        for _ in range(12):
            S = set()
            for _ in range(random.randint(1, 3)):
                a = random.randrange(L)
                c = min(L - 1, a + random.randrange(0, L // 3 + 1))
                S |= set(range(a, c + 1))
            mu = rng.integers(0, 2, 128)
            ct = sc.encrypt(S, mu)
            for x in random.sample(range(L), 6):
                out = sc.decrypt(sc.keygen(x), ct, x)
                tests += 1
                if x in S:
                    good = out is not None and np.array_equal(out, mu)
                    ok += int(good)
                    bad += int(not good)
                else:
                    false_pos += int(out is not None)
    sc = SCABE(2, 6, ibe, seed=SEED)
    mu = rng.integers(0, 2, 128)
    hdr, body = sc.encrypt_nodes([(3, 2)], mu)
    K_wrong = ibe.dec(ibe.keygen(node_id((3, 3), 2)), hdr[(3, 2)])
    mism = int(np.sum(K_wrong != (body ^ mu)))
    return dict(trials=tests, ok=ok, failures=bad, false_accepts=false_pos, wrong_key_bit_mismatch=mism, kappa=128)


def bench(n, qbits, kappa=128, reps=3):
    q = prime_below(2 ** qbits)
    ibe = IBE(n, q, kappa, seed=SEED)
    t = time.perf_counter(); ibe.setup(); t_setup = time.perf_counter() - t
    t_kg = []
    for r in range(reps):
        t = time.perf_counter(); ibe.keygen(1000 + r); t_kg.append(time.perf_counter() - t)
    K = rng.integers(0, 2, kappa)
    t_en = []
    for r in range(reps + 2):
        t = time.perf_counter(); ct = ibe.enc(77, K); t_en.append(time.perf_counter() - t)
    E7 = ibe.keygen(77)
    t_de = []
    for r in range(reps + 2):
        t = time.perf_counter(); out = ibe.dec(E7, ct); t_de.append(time.perf_counter() - t)
    correct = bool(np.array_equal(out, K))
    sc = SCABE(2, 20, ibe, seed=SEED)
    sc.encrypt_nodes([(20, 1)], rng.integers(0, 2, kappa))   # warm-up
    series = []
    for N in [1, 2, 4, 8, 16, 32, 64]:
        P = [(20, 5 * i + 1) for i in range(N)]
        mu = rng.integers(0, 2, kappa)
        ts = []
        for _ in range(7):
            t = time.perf_counter(); sc.encrypt_nodes(P, mu); ts.append(time.perf_counter() - t)
        series.append(dict(N=N, ms=1000 * min(ts)))
    return dict(n=n, qbits=qbits, q=q, k=ibe.k, kappa=kappa, mbar=ibe.mb, w=ibe.w, mtot=ibe.mt, correct=correct,
                setup_s=t_setup, keygen_ms=1000 * min(t_kg), enc_ms=1000 * min(t_en), dec_ms=1000 * min(t_de),
                enc_series=series)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fewer samples (smoke test)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "results.json"))
    ap.add_argument("--only-fl", action="store_true", help="recompute only the federated-learning cohort table and merge it into --out")
    a = ap.parse_args()
    T = 20 if a.quick else 300
    if a.only_fl:
        res = json.load(open(a.out)); res["fl_cohorts"] = fl_cohorts(T)
        json.dump(res, open(a.out, "w"), indent=1); print("merged fl_cohorts into", a.out); return
    out = {}
    out["exhaustive"] = [dict(b=b, d=d, L=L, worst=w, bound=bd)
                         for (b, d) in [(2, 6), (2, 8), (3, 4), (4, 3), (5, 3)] for (L, w, bd) in [exhaustive(b, d)]]
    print("exhaustive ok", out["exhaustive"], flush=True)
    out["cover_stats"] = cover_stats(T)
    out["interval_vs_L"] = single_interval_vs_L(T)
    out["arity"] = arity_tradeoff(T)
    out["revocation"] = revocation(max(5, T // 5))
    cs = out["cover_stats"]   # reuse the r=1 rows so that the interval statistics share one sample
    for r in out["interval_vs_L"]:
        if r["d"] == 20:
            r["mean"] = cs["b2"]["rows"][0]["mean"]
    for r in out["arity"]:
        k = f"b{r['b']}"
        if k in cs:
            r["ctmean"], r["ctmax"] = cs[k]["rows"][0]["mean"], cs[k]["rows"][0]["mx"]
    out["fl_cohorts"] = fl_cohorts(T)
    out["toy"] = toy_correctness()
    print("toy", out["toy"], flush=True)
    out["bench"] = {"P1": bench(256, 16), "P2": bench(512, 20)}
    out["env"] = dict(python=sys.version.split()[0], numpy=np.__version__, cpus=os.cpu_count())
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
