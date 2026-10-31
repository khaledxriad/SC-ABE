"""b-ary attribute tree, canonical covers and root-to-leaf paths (Algorithm 1 of the paper)."""
from typing import Iterable, List, Tuple

Node = Tuple[int, int]  # (level j, index i), 0 <= i < b**j


def leaves(node: Node, b: int, d: int) -> Tuple[int, int]:
    """Inclusive leaf interval below `node`."""
    j, i = node
    w = b ** (d - j)
    return i * w, (i + 1) * w - 1


def cover(a: int, c: int, b: int, d: int) -> List[Node]:
    """Canonical cover of the leaf interval [a, c]: maximal nodes whose leaves lie inside [a, c]."""
    if not (0 <= a <= c < b ** d):
        raise ValueError("interval outside the universe")
    out: List[Node] = []
    stack: List[Node] = [(0, 0)]
    while stack:
        nd = stack.pop()
        lo, hi = leaves(nd, b, d)
        if hi < a or lo > c:
            continue
        if a <= lo and hi <= c:
            out.append(nd)
            continue
        j, i = nd
        for t in range(b):
            stack.append((j + 1, i * b + t))
    return out


def anc(x: int, b: int, d: int) -> List[Node]:
    """The d+1 nodes on the root-to-leaf path of attribute x."""
    return [(j, x // b ** (d - j)) for j in range(d + 1)]


def runs(values: Iterable[int]) -> List[List[int]]:
    """Maximal runs [lo, hi] of consecutive integers in `values`."""
    out: List[List[int]] = []
    for x in sorted(set(values)):
        if out and x == out[-1][1] + 1:
            out[-1][1] = x
        else:
            out.append([x, x])
    return out


def cover_set(values: Iterable[int], b: int, d: int) -> List[Node]:
    """Canonical cover of an arbitrary attribute set (union of the covers of its runs)."""
    res: List[Node] = []
    for lo, hi in runs(values):
        res += cover(lo, hi, b, d)
    return res


def node_id(node: Node, b: int) -> int:
    """Injective heap-style encoding of a node as an IBE identity."""
    j, i = node
    return (b ** j - 1) // (b - 1) + i
