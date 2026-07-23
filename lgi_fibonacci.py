"""
Leggett-Garg test for Fibonacci-anyon braiding.

Setup: a single anyon pair (fusion-channel qubit, dim = 2) in the
standard B_3 Fibonacci representation on three tau anyons with total
charge tau. Topological charge Q is measured at three times t_1, t_2,
t_3 with Q = sigma_z, eigenvalues +/- 1. Between the measurements a
braid word B is applied (same evolution t_1 -> t_2 and t_2 -> t_3,
i.e. t_1 -> t_3 = B^2).

    K3 = C_12 + C_23 - C_13 = 2 C(B) - C(B^2)

    Macrorealism:  K3 <= 1
    Lueders:       K3 <= 3/2
    Algebraic:     K3 <= 3

Two-time correlator (projective Lueders measurement, maximally mixed
initial state):

    C(U) = (1/2) Re Tr[ Z U Z U^dag ],   Z = sigma_z

Standard Fibonacci F- and R-symbols (Nayak, Simon, Stern, Freedman,
Das Sarma, Rev. Mod. Phys. 80, 1083 (2008)):

    phi  = (1 + sqrt 5) / 2
    F    = [[1/phi, 1/sqrt(phi)], [1/sqrt(phi), -1/phi]]
    R_1  = exp(-4 pi i / 5)
    R_tau= exp(+3 pi i / 5)
    sigma_1 = diag(R_1, R_tau)
    sigma_2 = F sigma_1 F

The sector phase delta is an algebraic deformation
R_tau -> R_tau * exp(i delta); delta = 0 recovers the standard
Fibonacci theory.
"""
import numpy as np
import json

np.random.seed(20260524)

phi = (1 + np.sqrt(5)) / 2
F = np.array([[1/phi, 1/np.sqrt(phi)],
              [1/np.sqrt(phi), -1/phi]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def generators(delta=0.0):
    """sigma_1, sigma_2 and inverses in the pure B_3 Fibonacci representation (dim 2)."""
    R1 = np.exp(-4j*np.pi/5)
    Rt = np.exp(3j*np.pi/5) * np.exp(1j*delta)
    s1 = np.diag([R1, Rt]).astype(complex)
    s2 = F @ s1 @ F
    gens = np.stack([s1, s1.conj().T, s2, s2.conj().T])  # sigma1, sigma1^-1, sigma2, sigma2^-1
    return gens


def C_of(U):
    """Two-time correlator C(U) = 0.5 Re Tr[Z U Z U^dag], batched over (N,2,2)."""
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))
    return 0.5 * np.real(M[:, 0, 0] + M[:, 1, 1])


# ----- Sanity checks against standard Fibonacci MTC conventions -------------------------------
g = generators(0.0)
s1, s2 = g[0], g[2]
unit_err = max(np.max(np.abs(s1 @ s1.conj().T - np.eye(2))),
               np.max(np.abs(s2 @ s2.conj().T - np.eye(2))))
yb = s1 @ s2 @ s1 - s2 @ s1 @ s2                       # Yang-Baxter
delta2 = np.linalg.matrix_power(s1 @ s2, 3)             # (s1 s2)^3, must be scalar
scalar_phase = delta2[0, 0]
scalar_err = np.max(np.abs(delta2 - scalar_phase*np.eye(2)))
expected = np.exp(2j*np.pi/5)
print("=== Sanity ===")
print(f"  Unitarity        max err = {unit_err:.2e}")
print(f"  Yang-Baxter      max err = {np.max(np.abs(yb)):.2e}")
print(f"  (s1 s2)^3 scalar max err = {scalar_err:.2e}")
print(f"  Phase = {scalar_phase:.6f}  expected exp(2pi i/5) = {expected:.6f}  "
      f"diff = {abs(scalar_phase-expected):.2e}")

# ----- Exhaustive search: best K3 per braid-word length L ---------------------
def best_K3_exhaustive(delta, Lmax):
    """All 4^L braid words; returns per L: (best_K3, best_word, C(B), C(B^2))."""
    gens = generators(delta)
    labels = ['s1', 'S1', 's2', 'S2']  # S = inverse
    words = gens.copy()                # L=1
    word_lbls = [[i] for i in range(4)]
    out = {}
    for L in range(1, Lmax+1):
        B2 = words @ words
        CB = C_of(words)
        CB2 = C_of(B2)
        K3 = 2*CB - CB2
        idx = int(np.argmax(K3))
        out[L] = dict(K3=float(K3[idx]), CB=float(CB[idx]), CB2=float(CB2[idx]),
                      word=' '.join(labels[i] for i in word_lbls[idx]),
                      n_words=len(words))
        if L < Lmax:
            words = np.matmul(words[:, None, :, :], gens[None, :, :, :]
                              ).reshape(-1, 2, 2)
            word_lbls = [w + [i] for w in word_lbls for i in range(4)]
    return out

print("\n=== Exhaustive LGI search, sector delta=0 ===")
Lmax = 11
res0 = best_K3_exhaustive(0.0, Lmax)
print(f"{'L':>3} {'#words':>9} {'best K3':>9} {'C(B)':>8} {'C(B^2)':>8}  word")
for L in range(1, Lmax+1):
    r = res0[L]
    print(f"{L:>3} {r['n_words']:>9} {r['K3']:>9.5f} {r['CB']:>8.4f} "
          f"{r['CB2']:>8.4f}  {r['word']}")

K3_by_L = [res0[L]['K3'] for L in range(1, Lmax+1)]
first_violation = next((L for L in range(1, Lmax+1) if res0[L]['K3'] > 1.0), None)
best_overall = max(range(1, Lmax+1), key=lambda L: res0[L]['K3'])
print(f"\n  LGI violation (K3>1) first at L = {first_violation}")
print(f"  Best K3 = {res0[best_overall]['K3']:.6f} at L = {best_overall}")
print(f"  Lueders bound 1.5 -> gap = {1.5 - res0[best_overall]['K3']:.6f}")

# Monotonic?
mono = all(K3_by_L[i] <= K3_by_L[i+1] + 1e-12 for i in range(len(K3_by_L)-1))
dips = [L for L in range(2, Lmax+1) if res0[L]['K3'] < res0[L-1]['K3'] - 1e-9]
print(f"  Monotonic in L? {mono}.  Dips at L = {dips if dips else 'none'}")

# ----- Random search to longer lengths: convergence toward 1.5 ------------
print("\n=== Random search L=12..40 (convergence to Lueders bound) ===")
gens0 = generators(0.0)
rng = np.random.default_rng(20260524)
conv = {}
for L in range(12, 41, 2):
    nsamp = 400_000
    idx = rng.integers(0, 4, size=(nsamp, L))
    W = np.tile(np.eye(2, dtype=complex), (nsamp, 1, 1))
    for k in range(L):
        W = np.matmul(W, gens0[idx[:, k]])
    K3 = 2*C_of(W) - C_of(np.matmul(W, W))
    conv[L] = float(np.max(K3))
    print(f"  L={L:>3}  best K3 (random {nsamp:,}) = {conv[L]:.6f}   "
          f"gap to 1.5 = {1.5-conv[L]:.6f}")

# ----- Sector sweep: is the LGI violation sector-dependent? ------------------
print("\n=== Sector sweep: best K3(delta), exhaustive L<=9 ===")
deltas = np.linspace(0, 2*np.pi, 49)
sector = []
for d in deltas:
    r = best_K3_exhaustive(d, 9)
    sector.append(r[9]['K3'])
sector = np.array(sector)
print(f"  K3(delta) range: [{sector.min():.4f}, {sector.max():.4f}]")
print(f"  spread over sectors = {sector.max()-sector.min():.4f}")
imax = int(np.argmax(sector)); imin = int(np.argmin(sector))
print(f"  Maximum at delta/pi = {deltas[imax]/np.pi:.3f}  (K3={sector[imax]:.4f})")
print(f"  Minimum at delta/pi = {deltas[imin]/np.pi:.3f}  (K3={sector[imin]:.4f})")

# ----- Save --------------------------------------------------------------
results = dict(
    sanity=dict(unitarity_err=float(unit_err),
                yang_baxter_err=float(np.max(np.abs(yb))),
                delta2_scalar_err=float(scalar_err)),
    exhaustive_delta0={str(L): res0[L] for L in res0},
    first_violation_L=first_violation,
    best=dict(L=best_overall, K3=res0[best_overall]['K3']),
    monotone=mono, dips=dips,
    random_convergence={str(L): conv[L] for L in conv},
    sector_sweep=dict(delta_over_pi=list(deltas/np.pi),
                      K3=list(map(float, sector)),
                      spread=float(sector.max()-sector.min())),
)
with open('lgi_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\n-> lgi_results.json written")
