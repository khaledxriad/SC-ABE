import itertools, os, sys
import numpy as np
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from scabe import IBE, SCABE, anc, cover, cover_set, leaves, node_id, prime_below, runs


@pytest.mark.parametrize("b,d", [(2, 5), (3, 3), (4, 3)])
def test_cover_exact_minimal_and_path_lemma(b, d):
    L = b ** d
    for a, c in itertools.combinations_with_replacement(range(L), 2):
        P = cover(a, c, b, d)
        leaf_set = sorted(x for nd in P for x in range(leaves(nd, b, d)[0], leaves(nd, b, d)[1] + 1))
        assert leaf_set == list(range(a, c + 1))
        assert len(P) <= max(1, 2 * (b - 1) * d)
        for x in range(L):
            hits = set(P) & set(anc(x, b, d))
            assert len(hits) == (1 if a <= x <= c else 0)


def test_node_id_injective():
    for b, d in [(2, 8), (3, 5), (5, 3)]:
        ids = [node_id((j, i), b) for j in range(d + 1) for i in range(b ** j)]
        assert len(set(ids)) == len(ids)


def test_even_numbers_need_half_the_leaves():
    L = 64
    assert len(cover_set(range(0, L, 2), 2, 6)) == L // 2
    assert len(runs(range(0, L, 2))) == L // 2


def test_end_to_end_correctness_and_rejection():
    ibe = IBE(16, prime_below(65537), 128, seed=1)
    ibe.setup()
    sc = SCABE(2, 5, ibe, seed=1)
    mu = np.random.default_rng(3).integers(0, 2, 128)
    S = set(range(5, 19))
    ct = sc.encrypt(S, mu)
    for x in range(32):
        out = sc.decrypt(sc.keygen(x), ct, x)
        if x in S:
            assert np.array_equal(out, mu)
        else:
            assert out is None
