"""Tag-based dual-Regev style IBE stand-in (single trapdoor, gadget decomposition).

Used only to exercise the SC-ABE code path and to time realistic per-node costs.
Not adaptively secure; the deterministic bit-decomposition sampler leaks the trapdoor.
Identities must satisfy 0 <= id < q - 1 so that distinct identities get distinct tags.
"""
import numpy as np


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13):
        if n % p == 0:
            return n == p
    i = 17
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def prime_below(x: int) -> int:
    while not is_prime(x):
        x -= 1
    return x


class IBE:
    def __init__(self, n: int, q: int, kappa: int, sigma: float = 3.2, seed: int = 2026):
        if not is_prime(q):
            raise ValueError("q must be prime (scalar tags need invertible differences)")
        self.n, self.q, self.kappa, self.sigma = n, q, kappa, sigma
        self.k = q.bit_length()
        self.w = n * self.k
        self.mb = 2 * n
        self.mt = self.mb + self.w
        self.rng = np.random.default_rng(seed)
        self.pow2 = np.tile(2 ** np.arange(self.k, dtype=np.int64), n)

    def setup(self) -> None:
        n, q, rng = self.n, self.q, self.rng
        self.Abar = rng.integers(0, q, (n, self.mb), dtype=np.int64)
        self.R = rng.integers(-1, 2, (self.mb, self.w)).astype(np.float64)     # trapdoor
        B = np.mod(np.rint(self.Abar.astype(np.float64) @ self.R).astype(np.int64), q)
        self.AbarT = np.ascontiguousarray(self.Abar.T)
        self.BT = np.ascontiguousarray(B.T)
        self.U = rng.integers(0, q, (n, self.kappa), dtype=np.int64)
        self.UT = np.ascontiguousarray(self.U.T)

    def tag(self, ident: int) -> int:
        return (ident % (self.q - 1)) + 1

    def keygen(self, ident: int) -> np.ndarray:
        q, h = self.q, self.tag(ident)
        V = (pow(h, -1, q) * self.U) % q
        Z = ((V[:, None, :] >> np.arange(self.k)[None, :, None]) & 1).reshape(self.w, self.kappa)
        top = -np.rint(self.R @ Z.astype(np.float64)).astype(np.int64)
        return np.vstack([top, Z.astype(np.int64)])

    def enc(self, ident: int, K: np.ndarray):
        q, n, rng = self.q, self.n, self.rng
        h = self.tag(ident)
        s = rng.integers(0, q, n, dtype=np.int64)
        c1a = (self.AbarT @ s) % q
        c1b = (self.BT @ s + h * (np.repeat(s, self.k) * self.pow2)) % q
        x = np.rint(rng.normal(0, self.sigma, self.mt)).astype(np.int64)
        c1 = (np.concatenate([c1a, c1b]) + x) % q
        y = np.rint(rng.normal(0, self.sigma, self.kappa)).astype(np.int64)
        c0 = (self.UT @ s + y + np.asarray(K, dtype=np.int64) * (q // 2)) % q
        return c0, c1

    def dec(self, E: np.ndarray, ct) -> np.ndarray:
        c0, c1 = ct
        q = self.q
        v = (c0 - E.T @ c1) % q
        return ((v > q // 4) & (v < 3 * q // 4)).astype(np.int64)
