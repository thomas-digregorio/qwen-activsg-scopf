"""
Regression test: verify the contingency K matrices in dc_n1_opf against
independent post-contingency DC power-flow computations for every N-1
contingency (all line outages + all generator outages), at a random
dispatch.

History: a column-indexing bug in the theta->flow map M (M[bus, col]
indexed by position instead of bus) silently produced wrong post-
contingency flows for years-in-the-making N-1 constraint sets. This test
guards against regressions.
"""
import numpy as np

import dc_n1_opf as m
from dc_n1_opf import contingency_matrices, build_contingency_list

RNG = np.random.default_rng(42)


def independent_post_contingency(p_g, out_line=None, out_gen=None):
    """Direct DC PF with frozen pre-contingency dispatch (outage gen = 0).
    Returns dict {line_id: flow} for active lines, or None if disconnected."""
    active = [li for li in range(m.m) if li != out_line]
    adj = np.zeros((m.n, m.n), bool)
    for li in active:
        f, t, _ = m.br[li]
        adj[f, t] = adj[t, f] = True
    adj[m.slack, m.slack] = True
    from scipy.sparse.csgraph import connected_components
    from scipy import sparse
    _, labels = connected_components(sparse.csr_matrix(adj), directed=False)
    if not np.all(labels == labels[m.slack]):
        return None

    S = np.zeros((m.n, m.n))
    for li in active:
        f, t, bl = m.br[li]
        S[f, t] -= bl
        S[t, f] -= bl
        S[f, f] += bl
        S[t, t] += bl
    inj = m.E @ p_g - m.loads
    if out_gen is not None:
        inj[m.bidx[out_gen]] -= p_g[m.gcol[out_gen]]
    ns = m.ns
    th = np.zeros(m.n)
    th[ns] = np.linalg.solve(S[np.ix_(ns, ns)], inj[ns])
    return {li: m.br[li][2] * (th[m.br[li][0]] - th[m.br[li][1]]) for li in active}


def main():
    # random dispatch inside generator limits
    p_g = np.array([
        RNG.uniform(*m.C.GENS[b][:2]) * 0.8 + 0.1 * m.C.GENS[b][1]
        for b in m.gen_buses
    ])
    p_g = np.clip(p_g, 0, None)

    n_checked = 0
    n_skipped = 0
    for ct in build_contingency_list():
        out_line, out_gen = ct["out_line"], ct["out_gen"]
        ref = independent_post_contingency(p_g, out_line, out_gen)
        got = contingency_matrices(out_line, out_gen)

        if ref is None:
            assert got is None, f"{ct['name']}: expected disconnected, got K"
            n_skipped += 1
            continue
        assert got is not None, f"{ct['name']}: expected connected, got None"
        K, active, E_c = got
        F = K @ (E_c @ p_g - m.loads)
        assert F.shape[0] == len(active), f"{ct['name']}: row count mismatch"
        for row, li in enumerate(active):
            err = abs(F[row] - ref[li])
            assert err < 1e-7, (
                f"{ct['name']}: line {m.bus_ids[m.br[li][0]]}-"
                f"{m.bus_ids[m.br[li][1]]} K-flow {F[row]:.6f} != "
                f"direct-PF {ref[li]:.6f} (err {err:.2e})"
            )
        n_checked += 1

    print(f"PASS: {n_checked} contingencies verified vs independent DC PF, "
          f"{n_skipped} disconnected (correctly skipped)")


if __name__ == "__main__":
    main()
