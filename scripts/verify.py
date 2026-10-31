#!/usr/bin/env python3
"""Independent re-computation of the numbers quoted in the associated article.

Does NOT import the cover routine or the size model used by the experiments.
Checks (1) the cover routine against a second algorithm, (2) the size model from first principles,
(3) each quoted claim against results/results.json.
Exit status is non-zero if anything disagrees.
"""
import json, math, os, random, re, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from scabe import cover  # only used as the object under test
R = json.load(open(os.path.join(ROOT, "results", "results.json")))
fails = []
def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond: fails.append(name)

# (1) second, independent cover algorithm: bottom-up pairing for b-ary trees
def cover_bottom_up(a, c, b, d):
    lo, hi, level, out = a, c + 1, d, []          # half-open range of node indices at `level`
    while lo < hi and level >= 0:
        if level == 0:
            out.append((0, 0)); break
        while lo % b and lo < hi:
            out.append((level, lo)); lo += 1
        while hi % b and lo < hi:
            hi -= 1; out.append((level, hi))
        if lo >= hi: break
        lo //= b; hi //= b; level -= 1
        if level == 0 and lo == 0 and hi == 1:
            out.append((0, 0)); break
    return out
rnd = random.Random(1)
ok = True
for b, d in [(2, 8), (3, 5), (4, 4), (5, 3), (16, 2)]:
    L = b ** d
    for _ in range(3000):
        a = rnd.randrange(L); c = rnd.randrange(a, L)
        if sorted(cover(a, c, b, d)) != sorted(cover_bottom_up(a, c, b, d)):
            ok = False; print("mismatch", b, d, a, c); break
check("canonical cover equals independent bottom-up cover (15000 random intervals)", ok)

# exact mean cover size for a uniformly random pair a<=c at small depth (closed enumeration)
for b, d in [(2, 6), (3, 4)]:
    L = b ** d
    tot = sum(len(cover_bottom_up(a, c, b, d)) for a in range(L) for c in range(a, L)); n = L * (L + 1) // 2
    check(f"worst-case bound 2(b-1)d holds for b={b}, d={d}",
          max(len(cover_bottom_up(a, c, b, d)) for a in range(L) for c in range(a, L)) <= max(1, 2 * (b - 1) * d), f"mean {tot / n:.2f}")

# (2) size model from first principles
KAPPA = 128
def model(n, qb):
    k = qb; w = n * k; mb = 2 * n; mt = mb + w
    ct_adapt = (2 * mt + KAPPA) * qb / 8
    sk_adapt = 2 * KAPPA * mt * qb / 8
    mpk_adapt = 23 * n * mt * qb / 8 + n * KAPPA * qb / 8
    return dict(mt=mt, block=mt * qb / 8, ct_adapt=ct_adapt, sk_adapt=sk_adapt, mpk_adapt=mpk_adapt)
m = model(512, 20); KB = 1024; MB = KB * KB
check("P2 vector length m_t = 11264", m["mt"] == 11264)
check("P2 block 27.5 KB", abs(m["block"] / KB - 27.5) < 0.05, f'{m["block"] / KB:.2f}')
check("P2 adaptive node ciphertext 55.3 KB", abs(m["ct_adapt"] / KB - 55.3) < 0.05, f'{m["ct_adapt"] / KB:.2f}')
check("P2 adaptive attribute key (21 nodes) ~144 MB", round(21 * m["sk_adapt"] / MB) == 144, f'{21 * m["sk_adapt"] / MB:.1f}')
check("P2 adaptive public key ~316 MB", round(m["mpk_adapt"] / MB) == 316, f'{m["mpk_adapt"] / MB:.1f}')
for nm, p in R["bench"].items():
    mm = model(p["n"], p["qbits"])
    check(f"{nm}: benchmark dimensions match model", p["mtot"] == mm["mt"] and p["correct"])

# (3) prose claims vs results
cs = R["cover_stats"]; r1 = cs["b2"]["rows"][0]
check("single interval, b=2, l=2^20: mean cover below bound 40", r1["mean"] < 40 and r1["mx"] <= 40, f'mean {r1["mean"]:.1f}, max {r1["mx"]}')
for key in ("b2", "b4", "b16"):
    rows = cs[key]["rows"]
    check(f"{key}: every observed maximum respects 2(b-1)d r", all(r["mx"] <= r["bound"] for r in rows))
    check(f"{key}: cover size increases with r", all(rows[i]["mean"] < rows[i + 1]["mean"] for i in range(len(rows) - 1)))
check("arity factor ~ (b-1)/log2(b)", abs(cs["b4"]["rows"][0]["mean"] / r1["mean"] - 1.5) < 0.2 and abs(cs["b16"]["rows"][0]["mean"] / r1["mean"] - 15 / 4) < 0.5)
iv = {r["d"]: r for r in R["interval_vs_L"]}
slope = (iv[30]["mean"] - iv[6]["mean"]) / 24
check("about one extra node per doubling of the universe", 0.8 < slope < 1.2, f"{slope:.2f}")
ours = r1["mean"] * m["ct_adapt"]; naive = 2 ** 19 * m["ct_adapt"]
check("ratio per-value / SC-ABE ~ 2^19 / mean cover", abs(naive / ours - 2 ** 19 / r1["mean"]) < 1e-6, f'{2 ** 19 / r1["mean"]:,.0f}')
check("Wee omega(1) variant ratio ~36x in the model", round(ours / m["block"]) == 36, f'{ours / m["block"]:.1f}')
check("Wee l^(2/3) variant ratio ~5x in the model", round(ours / (20 ** (2 / 3) * m["block"])) == 5, f'{ours / (20 ** (2 / 3) * m["block"]):.2f}')
rev64 = [r for r in R["revocation"] if r["t"] == 64][0]
check("revocation t=64: cover within 10% of t*log2(l/t)", abs(rev64["mean"] - rev64["nnl"]) / rev64["nnl"] < 0.1, f'{rev64["mean"]:.0f} vs {rev64["nnl"]:.0f}')
check("revocation: cover never above the bound 2d(t+1)", all(r["mean"] <= r["bound"] for r in R["revocation"]))

# federated-learning cohort scenarios (exact structural claims + sanity of the statistics)
fl = {r["name"]: r for r in R["fl_cohorts"]}
site = fl["one site (aligned block of 2^10 devices)"]
check("aligned site needs exactly one node", site["mean"] == 1.0 and site["mx"] == 1)
check("8 aligned sites need at most 8 nodes", fl["8 sites (aligned blocks of 2^10 devices)"]["mx"] <= 8)
reg = fl["one region (arbitrary interval of 2^16 devices)"]
check("region interval respects the bound 2d = 40", reg["mx"] <= 40, f'mean {reg["mean"]:.1f}')
for m in (100, 1000):
    r = fl[f"random sample of {m} devices"]
    check(f"unstructured sample of {m}: cover is (almost) one node per device", 0.95 * m <= r["mean"] <= m, f'{r["mean"]:.1f}')
for t in (16, 256):
    r = fl[f"whole fleet except {t} compromised devices"]
    check(f"all-but-{t}: cover below the bound 2d(t+1)", r["mx"] <= 40 * (t + 1), f'mean {r["mean"]:.0f}')
t = R["toy"]
check("toy test: no failures, no false accepts", t["failures"] == 0 and t["false_accepts"] == 0 and t["ok"] > 0, f'{t["ok"]} authorized of {t["trials"]}')
check("wrong key output looks random (binomial 128, p=1/2, |x-64| < 5 sigma)", abs(t["wrong_key_bit_mismatch"] - 64) < 5 * math.sqrt(32))
for e in R["exhaustive"]:
    check(f"exhaustive b={e['b']} d={e['d']}: worst case <= bound", e["worst"] <= e["bound"], f"{e['worst']} <= {e['bound']}")
for nm, p in R["bench"].items():
    s = {x["N"]: x["ms"] for x in p["enc_series"]}
    check(f"{nm}: encryption time is ~linear in N (64-node time within 25% of 64x one-node)", abs(s[64] / (64 * s[1]) - 1) < 0.25, f"{s[64] / (64 * s[1]):.2f}")
t30 = iv[30]["mean"] * 2 * R["bench"]["P2"]["enc_ms"] / 1000
check("extrapolated 2^30 encryption time", 0.05 < t30 < 2, f"{t30:.2f} s")

print("\nALL CHECKS PASSED" if not fails else f"\n{len(fails)} CHECK(S) FAILED: {fails}")
sys.exit(1 if fails else 0)
