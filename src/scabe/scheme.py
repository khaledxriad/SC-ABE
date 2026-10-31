"""SC-ABE (Algorithms 2-5 of the paper) on top of an IBE object."""
import numpy as np
from .tree import anc, cover_set, node_id


class SCABE:
    def __init__(self, b: int, d: int, ibe, seed: int = 7):
        self.b, self.d, self.ibe = b, d, ibe
        self.rng = np.random.default_rng(seed)

    def keygen(self, x: int):
        return {nd: self.ibe.keygen(node_id(nd, self.b)) for nd in anc(x, self.b, self.d)}

    def encrypt_nodes(self, P, mu):
        K = self.rng.integers(0, 2, self.ibe.kappa)
        hdr = {nd: self.ibe.enc(node_id(nd, self.b), K) for nd in P}
        return hdr, (np.asarray(mu) ^ K)

    def encrypt(self, S, mu):
        return self.encrypt_nodes(cover_set(S, self.b, self.d), mu)

    def decrypt(self, sk, ct, x: int):
        hdr, body = ct
        for nd in anc(x, self.b, self.d):
            if nd in hdr:
                return body ^ self.ibe.dec(sk[nd], hdr[nd])
        return None
