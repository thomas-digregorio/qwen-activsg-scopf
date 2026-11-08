"""
Regression test: verify the ACTIVSg500 contingency K matrices against
independent post-contingency DC power-flow computations for a sample of
contingencies (line outages + generator outages), at a random dispatch.

Same guard as test_contingency_k.py for the 14-bus case.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

import activsg500_n1_opf as d

RNG = np.random.default_rng(42)


def independent_post_flow(M, p_g, out_line=None, out_gen=None):
    """Direct DC PF with frozen dispatch (outage gen set to 0), on the
    slack-side component. Returns {line_idx: flow_MW} for active
    slack-side lines, or None if the slack-side component is empty."""
    n, m = M["n"], M["m"]
    all_active = [li for li in range(m) if li != out_line]
    adj = np.zeros((n, n), bool)
    for li in all_active:
        f, t = M["br"][li][0], M["br"][li][1]
        adj[f, t] = adj[t, f] = True
    adj[M["slack"], M["slack"]] = True
    _, labels = connected_components(sparse.csr_matrix(adj), directed=False)
    sl = labels[M["slack"]]
    active = [li for li in all_active if labels[M["br"][li][0]] == sl]
    if not active:
        return None
    ns_c = [i for i in range(n) if i != M["slack"] and labels[i] == sl]
    pos = {i: k for k, i in enumerate(ns_c)}
    S = np.zeros((len(ns_c), len(ns_c)))
    for li in active:
        f, t, bl = M["br"][li][0], M["br"][li][1], M["br"][li][2]
        if f in pos:
            S[pos[f], pos[f]] += bl
        if t in pos:
            S[pos[t], pos[t]] += bl
        if f in pos and t in pos:
            S[pos[f], pos[t]] -= bl
            S[pos[t], pos[f]] -= bl
    inj = M["E"] @ p_g - M["loads"]
    if out_gen is not None:
        # outage unit's frozen output is lost at ITS BUS (not the gen column)
        inj[M["bidx"][M["gen_buses"][out_gen]]] -= p_g[out_gen]
    th = np.zeros(n)
    th[ns_c] = np.linalg.solve(S, inj[ns_c])
    # hybrid-unit flows: S in 1/pu, injections in MW -> flows in MW directly
    return {li: M["br"][li][2] * (th[M["br"][li][0]] - th[M["br"][li][1]])
            for li in active}


def main():
    M = d.build_model()
    Kme, Kload, active, valid, island, _, _ = d.precompute_contingencies(M)

    # random dispatch inside generator limits
    p_g = np.array([RNG.uniform(M["gen_pmin"][k], M["gen_pmax"][k])
                    for k in range(M["ng"])])
    p_g = np.clip(p_g, M["gen_pmin"], M["gen_pmax"])

    # sample: 25 random line contingencies + 10 random generator contingencies
    n_line = len(M["line_cts"])
    sample = list(RNG.choice(n_line, size=25, replace=False)) + \
        [n_line + i for i in RNG.choice(len(M["gen_cts"]), size=10, replace=False)]

    n_checked = n_skipped = 0
    for c in sample:
        out_line = M["line_cts"][c] if c < n_line else None
        out_gen = (M["gen_cts"][c - n_line] if c >= n_line else None)
        ref = independent_post_flow(M, p_g, out_line, out_gen)
        res = d.contingency_dispatch_map(M, out_line, out_gen)
        tag = f"line {M['bus_ids'][M['br'][out_line][0]]}-{M['bus_ids'][M['br'][out_line][1]]}" \
            if c < n_line else f"gen bus {M['gen_buses'][out_gen]}"
        if ref is None:
            assert res is None, f"{tag}: expected disconnected, got K"
            n_skipped += 1
            continue
        assert res is not None, f"{tag}: expected connected, got None"
        Kme_c, Kload_c, act, _, _, _, _ = res
        F = Kme_c @ p_g - Kload_c
        err = 0.0
        for r, li in enumerate(act):
            err = max(err, abs(F[r] - ref[li]))
        assert err < 1e-6 * max(1.0, max(abs(v) for v in ref.values())), \
            f"{tag}: K-flow vs independent PF error {err:.3e} MW"
        # precomputed table agrees with on-the-fly computation
        act3 = active[c]
        act3 = act3[act3 >= 0]
        assert list(act3) == list(act), f"{tag}: active line list mismatch"
        assert np.allclose(Kme[c, act3, :], Kme_c) and np.allclose(Kload[c, act3], Kload_c), \
            f"{tag}: precomputed K mismatch"
        n_checked += 1

    print(f"PASS: {n_checked} contingencies verified vs independent DC PF, "
          f"{n_skipped} disconnected (correctly skipped)")


if __name__ == "__main__":
    main()
