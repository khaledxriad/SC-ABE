# SC-ABE: reference implementation and experiments

Code repository associated with the article

> **SC-ABE: Compact lattice-based attribute-based encryption for cohort access control in federated cyber-physical systems**
> Khaled Riad, Computer Science Department, College of Computer Sciences & Information Technology,
> King Faisal University, Al-Ahsa, Saudi Arabia (kriad@kfu.edu.sa).

SC-ABE (subtree-cover attribute-based encryption) compiles any adaptively secure lattice identity-based
encryption (IBE) scheme into an ABE scheme for interval and subtree policies: the ciphertext carries one IBE
ciphertext per node of a canonical subtree cover, so an interval over `l = b^d` attribute values needs at
most `2(b-1)d` nodes. The article studies it for cohort access control in federated learning over cyber-physical systems; `experiments/run_experiments.py` includes the device-cohort scenarios (`fl_cohorts`) of that use case. No federated training is included.

> **Security warning.** `src/scabe/ibe.py` is a *tag-based dual-Regev stand-in* used for functional tests and
> per-node timing. It is **not** adaptively secure, uses a non-cryptographic random generator and a
> deterministic bit decomposition that leaks its trapdoor. Do not use it to protect real data.

## Layout

| Path | Content |
|---|---|
| `src/scabe/` | tree/cover routines (`tree.py`), prototype IBE (`ibe.py`), SC-ABE (`scheme.py`) |
| `tests/` | unit tests: cover exactness/minimality, path lemma, bound, end-to-end correctness |
| `experiments/run_experiments.py` | produces every measured number -> `results/results.json` |
| `results/results.json` | results reported in the article (seed 2026; timings are machine dependent) |
| `scripts/verify.py` | independent re-computation and consistency checks of the reported numbers |

## Reproduce

Requirements: Python >= 3.10, NumPy and pytest.

```bash
pip install -r requirements.txt   # or: pip install -e .
make test       # unit tests
make results    # about 1-2 minutes; `python experiments/run_experiments.py --quick` is a smoke test
make verify     # independent checks of the reported numbers
```

Counting results (cover sizes, test outcomes) are deterministic for a fixed seed; wall-clock timings depend
on the machine.

## Quick example

```python
import numpy as np
from scabe import IBE, SCABE, prime_below

ibe = IBE(n=16, q=prime_below(65537), kappa=128); ibe.setup()
abe = SCABE(b=2, d=5, ibe=ibe)
mu = np.random.default_rng(0).integers(0, 2, 128)
ct = abe.encrypt(set(range(5, 19)), mu)        # policy: attributes 5..18
assert np.array_equal(abe.decrypt(abe.keygen(7), ct, 7), mu)
assert abe.decrypt(abe.keygen(25), ct, 25) is None
```

## Citation and license

See `CITATION.cff` (add the Zenodo DOI once the archive exists). MIT license (see `LICENSE`).
