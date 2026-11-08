"""Validate LMP-from-duals sign convention on the 14-bus N-1 OPF via
highspy passLp (CSR). Compares row duals of the nodal-balance rows
against finite-difference LMPs from the scipy pipeline."""
import numpy as np
from scipy import sparse
import highspy
import dc_n1_opf as m

n, ng = m.n, m.ng
nseg = 12
nseg_total = ng * nseg
dim = ng + nseg_total + (n - 1)
s0, t0 = ng, ng + nseg_total
Fmax = 100.0

seg_cost, seg_width, seg_price = m.build_pl_costs(nseg)
c_obj = np.zeros(dim)
c_obj[s0:t0] = seg_price

Sall_ns = m.Sfull[:, m.ns]
a_eq, b_eq = [], []
for i in range(n):
    row = np.zeros(dim)
    row[0:ng] = m.E[i, :]
    row[t0:t0 + n - 1] = -Sall_ns[i, :]
    a_eq.append(row)
    b_eq.append(m.loads[i])
for k in range(ng):
    row = np.zeros(dim)
    row[k] = 1.0
    row[s0 + k * nseg:s0 + (k + 1) * nseg] = -1.0
    a_eq.append(row)
    b_eq.append(0.0)

a_ub, b_ub = [], []
for (f, t, bl) in m.br:
    row = np.zeros(dim)
    if f in m.ns_pos:
        row[t0 + m.ns_pos[f]] += bl
    if t in m.ns_pos:
        row[t0 + m.ns_pos[t]] -= bl
    a_ub += [row, -row]
    b_ub += [Fmax, Fmax]
for ct in m.build_contingency_list():
    res = m.contingency_matrices(ct["out_line"], ct["out_gen"])
    if res is None:
        continue
    K, active, E_c = res
    Kme = K @ E_c
    Kmloads = K @ m.loads
    for r in range(Kme.shape[0]):
        row = np.zeros(dim)
        row[0:ng] = Kme[r, :]
        a_ub += [row, -row]
        b_ub += [Fmax + Kmloads[r], Fmax - Kmloads[r]]

lb = np.zeros(dim)
ub = np.full(dim, np.inf)
for k, b in enumerate(m.gen_buses):
    pmin, pmax, _, _, _ = m.C.GENS[b]
    lb[k], ub[k] = pmin, pmax
for k in range(ng):
    for j in range(nseg):
        ub[s0 + k * nseg + j] = seg_width[k * nseg + j]
lb[t0:] = -np.inf

n_eq = len(a_eq)
A = sparse.vstack([sparse.csr_matrix(np.array(a_eq)),
                   sparse.csr_matrix(np.array(a_ub))], format="csr")
rlb = np.concatenate([np.array(b_eq), np.full(len(a_ub), -highspy.kHighsInf)])
rub = np.concatenate([np.array(b_eq), np.array(b_ub)])


def solve_highspy(c, lb, ub, A, rlb, rub):
    h = highspy.Highs()
    opts = highspy.HighsOptions()
    opts.output_flag = True
    h.passOptions(opts)
    lp = highspy.HighsLp()
    lp.num_row_ = A.shape[0]
    lp.num_col_ = A.shape[1]
    lp.sense_ = highspy.ObjSense.kMinimize
    lp.col_cost_ = c.tolist()
    lp.col_lower_ = lb.tolist()
    lp.col_upper_ = ub.tolist()
    lp.row_lower_ = rlb.tolist()
    lp.row_upper_ = rub.tolist()
    mat = highspy.HighsSparseMatrix()
    mat.format_ = highspy.MatrixFormat.kRowwise
    mat.num_row_ = A.shape[0]
    mat.num_col_ = A.shape[1]
    mat.start_ = A.indptr.tolist()
    mat.index_ = A.indices.tolist()
    mat.value_ = A.data.tolist()
    lp.a_matrix_ = mat
    model = highspy.HighsModel()
    model.lp_ = lp
    st = h.passModel(model)
    assert st in (highspy.HighsStatus.kOk, highspy.HighsStatus.kWarning), \
        f"passModel failed: {st}"
    h.run()
    ms = h.getModelStatus()
    sol = h.getSolution()
    x = np.array(sol.col_value, dtype=float)
    duals = np.array(sol.row_dual, dtype=float)
    ok = ms == highspy.HighsModelStatus.kOptimal
    return ok, x, duals, h


ok, x, duals, h = solve_highspy(c_obj, lb, ub, A, rlb, rub)
print("status optimal:", ok)
print("objective:", h.getObjectiveValue() if hasattr(h, "getObjectiveValue") else None)
print("gen:", np.round(x[:ng], 2), "sum:", round(float(x[:ng].sum()), 3))

y_eq = duals[:n_eq]          # first 14 rows = nodal balance
lmp_fd = m.compute_lmps(m.loads, Fmax, nseg)
print(f"max |dual - FD|  = {np.abs(y_eq[:14] - lmp_fd).max():.5f}")
print(f"max |-dual - FD| = {np.abs(-y_eq[:14] - lmp_fd).max():.5f}")
print("FD LMPs:", np.round(lmp_fd, 3))
print("duals  :", np.round(y_eq[:14], 3))
