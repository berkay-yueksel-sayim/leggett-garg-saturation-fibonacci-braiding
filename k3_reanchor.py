#!/usr/bin/env python3
"""
Paper 3 -- K3 argmax-word re-anchor (v1.2 correction; two independent routes)
============================================================================
The printed Z.299-301 sentence claims the L=11 exhaustive maximum
K3=1.499762 is "attained by the word"

    sigma_2 sigma_1 sigma_2^2 sigma_1^-1 sigma_2^2 sigma_1^-1 sigma_2^2 sigma_1^-1 sigma_2^2 sigma_1

which decodes to 14 letters (not 11) and, on inspection, is not the
word recorded for L=11 in the deposited lgi_results.json (produced by
lgi_fibonacci.py, this same deposit) -- that archived word is

    s1 s2 s2 S1 S1 S1 s2 s2 S1 s2 s1        (11 letters, K3=1.499761961900839)

This script independently re-derives the true L=11 argmax from the
paper's own engine code (lgi_fibonacci.py, same directory, NOT edited
in place -- it is a published Zenodo-live artifact) via two
algorithmically distinct routes that must agree bit-exactly:

  Route (i)  -- incremental engine-faithful enumeration: identical to
                best_K3_exhaustive() in lgi_fibonacci.py (build word
                L+1 from word L by right-multiplying every one of the
                4 generators). Cross-checked against the ALREADY
                archived lgi_results.json for L=1..10 (must match
                bit-exactly) before the L=11 output is trusted.
  Route (ii) -- meet-in-the-middle: enumerate all 4^5 length-5 and all
                4^6 length-6 words independently, combine every pair
                by matrix product (associativity guarantees this
                equals the direct L=11 product for the concatenated
                label sequence), and take the max over all 4^11
                combinations. This shares no incremental state with
                Route (i) and would not reproduce a Route-(i) bug
                (e.g. an off-by-one in the label bookkeeping).

Positive controls (must all pass before the L=11 result is trusted):
  - generators(0.0) unitary, Yang-Baxter s1 s2 s1 = s2 s1 s2, and
    (s1 s2)^3 scalar -- same 3 sanity checks as lgi_fibonacci.py.
  - C(U=I) = 1 exactly (L=0 boundary case: two-time correlator of a
    trivial evolution saturates the algebraic ceiling of the
    macrorealism bound, K3(I)=1).
  - Route (i) reproduces every archived best-K3(L) for L=1..10 in
    lgi_results.json bit-exactly (rtol/atol 1e-9).
  - Route (i) and Route (ii) agree on the L=11 maximum K3 bit-exactly
    (atol 1e-9), independent of whether they report the same argmax
    word (ties are expected -- the printed sentence itself says "and
    by several equivalent words").

Output: k3_reanchor.json (both routes' L=11 results, tie count,
archived-JSON cross-check, sigma-notation transcription).
"""
from __future__ import annotations
import json
import itertools
from pathlib import Path
import numpy as np

OUT = Path(__file__).parent
ARCHIVED = json.loads((OUT / "lgi_results.json").read_text())

# ------------------------------------------------------------------
# Generators -- identical construction to lgi_fibonacci.py
# ------------------------------------------------------------------
phi = (1 + np.sqrt(5)) / 2
F = np.array([[1/phi, 1/np.sqrt(phi)],
              [1/np.sqrt(phi), -1/phi]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def generators(delta: float = 0.0) -> np.ndarray:
    R1 = np.exp(-4j * np.pi / 5)
    Rt = np.exp(3j * np.pi / 5) * np.exp(1j * delta)
    s1 = np.diag([R1, Rt]).astype(complex)
    s2 = F @ s1 @ F
    return np.stack([s1, s1.conj().T, s2, s2.conj().T])  # s1, S1, s2, S2


def C_of(U: np.ndarray) -> np.ndarray:
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))
    return 0.5 * np.real(M[:, 0, 0] + M[:, 1, 1])


LABELS = ['s1', 'S1', 's2', 'S2']

# ------------------------------------------------------------------
# Positive controls
# ------------------------------------------------------------------
controls = {}
g = generators(0.0)
s1, S1, s2, S2 = g[0], g[1], g[2], g[3]

unit_err = max(float(np.max(np.abs(s1 @ s1.conj().T - np.eye(2)))),
               float(np.max(np.abs(s2 @ s2.conj().T - np.eye(2)))))
controls["unitarity_max_err"] = {"value": unit_err, "pass": unit_err < 1e-10}

yb_err = float(np.max(np.abs(s1 @ s2 @ s1 - s2 @ s1 @ s2)))
controls["yang_baxter_max_err"] = {"value": yb_err, "pass": yb_err < 1e-10}

delta2 = np.linalg.matrix_power(s1 @ s2, 3)
scalar_err = float(np.max(np.abs(delta2 - delta2[0, 0] * np.eye(2))))
controls["s1s2_cubed_scalar_max_err"] = {"value": scalar_err, "pass": scalar_err < 1e-10}

I2 = np.eye(2, dtype=complex)
CB_I = float(C_of(I2[None, :, :])[0])
K3_I = 2 * CB_I - CB_I
controls["C(U=I)==1_and_K3(I)==1"] = {
    "C_of_I": CB_I, "K3_of_I": K3_I,
    "pass": abs(CB_I - 1.0) < 1e-12 and abs(K3_I - 1.0) < 1e-12,
}

assert all(c["pass"] for c in controls.values()), f"POSITIVE CONTROL FAILED: {controls}"

# ------------------------------------------------------------------
# Route (i): incremental engine-faithful enumeration, L=1..11,
# cross-checked against archived lgi_results.json for L=1..10.
# ------------------------------------------------------------------
gens = generators(0.0)
words = gens.copy()
word_lbls = [[i] for i in range(4)]
route_i = {}
Lmax = 11
for L in range(1, Lmax + 1):
    B2 = np.matmul(words, words)
    CB = C_of(words)
    CB2 = C_of(B2)
    K3 = 2 * CB - CB2
    idx = int(np.argmax(K3))
    route_i[L] = dict(K3=float(K3[idx]), CB=float(CB[idx]), CB2=float(CB2[idx]),
                       word=' '.join(LABELS[i] for i in word_lbls[idx]),
                       n_words=len(words))
    if L == Lmax:
        K3_L11_full = K3          # keep full array for tie-counting
        lbls_L11_full = word_lbls
    if L < Lmax:
        words = np.matmul(words[:, None, :, :], gens[None, :, :, :]).reshape(-1, 2, 2)
        word_lbls = [w + [i] for w in word_lbls for i in range(4)]

cross_check = {}
for L in range(1, 11):
    arch = ARCHIVED["exhaustive_delta0"][str(L)]
    diff = abs(route_i[L]["K3"] - arch["K3"])
    cross_check[f"L={L}"] = {"route_i_K3": route_i[L]["K3"], "archived_K3": arch["K3"],
                              "diff": diff, "pass": diff < 1e-9}
assert all(c["pass"] for c in cross_check.values()), f"ROUTE-I / ARCHIVE MISMATCH: {cross_check}"

# Tie count at L=11 (within 1e-9 of the max)
k3_max_i = float(K3_L11_full.max())
tie_idx = np.where(np.abs(K3_L11_full - k3_max_i) < 1e-9)[0]
n_ties = int(tie_idx.size)
route_i_word = route_i[Lmax]["word"]

# ------------------------------------------------------------------
# Route (ii): meet-in-the-middle, L=11 = 5 + 6, independent of the
# incremental construction above (only associativity is shared).
# ------------------------------------------------------------------
def enumerate_words(L: int):
    mats = gens.copy()
    lbls = [[i] for i in range(4)]
    for _ in range(L - 1):
        mats = np.matmul(mats[:, None, :, :], gens[None, :, :, :]).reshape(-1, 2, 2)
        lbls = [w + [i] for w in lbls for i in range(4)]
    return mats, lbls

left_mats, left_lbls = enumerate_words(5)     # 4^5 = 1024
right_mats, right_lbls = enumerate_words(6)   # 4^6 = 4096

best_ii = {"K3": -np.inf, "i": None, "j": None}
for i in range(left_mats.shape[0]):
    combined = np.matmul(left_mats[i], right_mats)      # (4096,2,2)
    combined2 = np.matmul(combined, combined)
    CBc = C_of(combined)
    CB2c = C_of(combined2)
    K3c = 2 * CBc - CB2c
    j = int(np.argmax(K3c))
    if K3c[j] > best_ii["K3"]:
        best_ii = {"K3": float(K3c[j]), "i": i, "j": j}

route_ii_lbl = left_lbls[best_ii["i"]] + right_lbls[best_ii["j"]]
route_ii_word = ' '.join(LABELS[k] for k in route_ii_lbl)
k3_max_ii = best_ii["K3"]

route_mismatch = abs(k3_max_i - k3_max_ii)
assert route_mismatch < 1e-9, f"Route (i)/(ii) disagree: {k3_max_i} vs {k3_max_ii}"
assert len(route_ii_lbl) == 11 and len(lbls_L11_full[int(np.argmax(K3_L11_full))]) == 11

# ------------------------------------------------------------------
# Sigma-notation transcription (exponent-compressed) of the Route-(i)
# argmax word, for direct use in the paper.
# ------------------------------------------------------------------
def to_sigma_notation(lbl: list[int]) -> str:
    sym = {0: ('1', '+'), 1: ('1', '-'), 2: ('2', '+'), 3: ('2', '-')}
    runs = []
    for k in lbl:
        base, sign = sym[k]
        if runs and runs[-1][0] == base and runs[-1][1] == sign:
            runs[-1][2] += 1
        else:
            runs.append([base, sign, 1])
    out = []
    for base, sign, n in runs:
        if sign == '+':
            out.append(f"\\sigma_{{{base}}}" + (f"^{{{n}}}" if n > 1 else ""))
        else:
            out.append(f"\\sigma_{{{base}}}^{{-{n}}}" if n > 1 else f"\\sigma_{{{base}}}^{{-1}}")
    return ''.join(out)

sigma_str_i = to_sigma_notation(lbls_L11_full[int(np.argmax(K3_L11_full))])
sigma_str_ii = to_sigma_notation(route_ii_lbl)

out = {
    "meta": {
        "script": "k3_reanchor.py",
        "purpose": "v1.2 re-anchor: re-derive the true L=11 K3-argmax word from lgi_fibonacci.py's own generators",
        "route_i": "incremental engine-faithful enumeration (matches best_K3_exhaustive in lgi_fibonacci.py)",
        "route_ii": "meet-in-the-middle, L=11=5+6, algorithmically independent of route_i",
        "route_mismatch_K3": route_mismatch,
        "n_words_tied_at_max_L11": n_ties,
    },
    "positive_controls": controls,
    "archive_cross_check_L1_to_10": cross_check,
    "route_i_L11": {"K3": k3_max_i, "word": route_i_word, "n_letters": len(route_i_word.split()),
                     "sigma_notation": sigma_str_i},
    "route_ii_L11": {"K3": k3_max_ii, "word": route_ii_word, "n_letters": len(route_ii_word.split()),
                      "sigma_notation": sigma_str_ii},
    "archived_lgi_results_L11": ARCHIVED["exhaustive_delta0"]["11"],
}
(OUT / "k3_reanchor.json").write_text(json.dumps(out, indent=2))

print("=== Positive controls ===")
for k, v in controls.items():
    print(f"  {k}: {v}")
print("\n=== Archive cross-check L=1..10 (route i vs lgi_results.json) ===")
for k, v in cross_check.items():
    print(f"  {k}: diff={v['diff']:.2e}  pass={v['pass']}")
print(f"\n=== L=11 ===")
print(f"  Route (i)  K3={k3_max_i:.9f}  word='{route_i_word}'  ({len(route_i_word.split())} letters)")
print(f"             sigma: {sigma_str_i}")
print(f"  Route (ii) K3={k3_max_ii:.9f}  word='{route_ii_word}'  ({len(route_ii_word.split())} letters)")
print(f"             sigma: {sigma_str_ii}")
print(f"  Route mismatch: {route_mismatch:.2e}")
print(f"  # words tied at max (within 1e-9): {n_ties}")
print(f"\nWRITE {OUT / 'k3_reanchor.json'}")
