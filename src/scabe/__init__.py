"""SC-ABE: subtree-cover attribute-based encryption (reference prototype).

The lattice IBE shipped here is a tag-based stand-in used for correctness checks and
timing only. It is NOT adaptively secure and NOT hardened; do not use it to protect data.
"""
from .tree import leaves, cover, anc, runs, cover_set, node_id
from .ibe import IBE, is_prime, prime_below
from .scheme import SCABE

__all__ = ["leaves", "cover", "anc", "runs", "cover_set", "node_id", "IBE", "is_prime", "prime_below", "SCABE"]
__version__ = "1.0.0"
