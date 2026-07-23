"""
Leggett-Garg test for Ising-anyon braiding.

Counterpart to lgi_fibonacci.py for the non-universal Ising-anyon model.
Setup is identical: one anyon-pair in the 2D fusion-space, observable
Q = sigma_z (eigenvalues +/-1), three projective Lueders measurements at
t_1, t_2, t_3 separated by braid word B (so t_1 -> t_3 = B^2).

  K3 = 2 C(B) - C(B^2),   C(U) = (1/2) Re Tr[Z U Z U^dag],   Z = sigma_z
  Macrorealism: K3 <= 1   Lueders: K3 <= 3/2   Algebraic: K3 <= 3

Ising-anyon generators in the 4-sigma fusion-space (Nayak/Simon/Stern/
Freedman/Das Sarma, Rev. Mod. Phys. 80, 1083 (2008), Section V.B):
  R^{sigma sigma}_1   = exp(-i pi/8)
  R^{sigma sigma}_psi = exp(+i 3 pi/8)
  F^{sigma sigma sigma sigma}_{sigma sigma} = (1/sqrt 2) [[1,1],[1,-1]]   (Hadamard)
  sigma_1 = diag(R_1, R_psi)
  sigma_2 = F sigma_1 F

Phase-deformation parameter delta (algebraic generalization, analogous to
the Fibonacci case): R_psi -> R_psi * exp(i delta). delta = 0 is the standard
Ising model.
"""
import numpy as np
import json

np.random.seed(20260525)

F = (1.0 / np.sqrt(2.0)) * np.array([[1, 1], [1, -1]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def generators(delta=0.0):
    R1 = np.exp(-1j * np.pi / 8)
    Rp = np.exp(3j * np.pi / 8) * np.exp(1j * delta)
    s1 = np.diag([R1, Rp]).astype(complex)
    s2 = F @ s1 @ F
    gens = np.stack([s1, s1.conj().T, s2, s2.conj().T])
    return gens


def C_of(U):
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))
    return 0.5 * np.real(M[:, 0, 0] + M[:, 1, 1])


# ----- Sanity checks against Ising-MTC consistency ---------------------------
g = generators(0.0)
s1, s2 = g[0], g[2]
unit_err = max(np.max(np.abs(s1 @ s1.conj().T - np.eye(2))),
               np.max(np.abs(s2 @ s2.conj().T - np.eye(2))))
yb = s1 @ s2 @ s1 - s2 @ s1 @ s2
# (sigma_1 sigma_2)^3 is a scalar; the Ising scalar phase is exp(-i pi/8) times
# a sign coming from the central charge c = 1/2: (s1 s2)^3 = exp(i pi/8) * I
# in the standard convention with R_1 = exp(-i pi/8). Verify it is scalar
# and report the actual phase to machine epsilon.
delta2 = np.linalg.matrix_power(s1 @ s2, 3)
scalar_phase = delta2[0, 0]
scalar_err = np.max(np.abs(delta2 - scalar_phase * np.eye(2)))
print("=== Sanity (Ising) ===")
print(f"  Unitarity        max err = {unit_err:.2e}")
print(f"  Yang-Baxter      max err = {np.max(np.abs(yb)):.2e}")
print(f"  (s1 s2)^3 scalar max err = {scalar_err:.2e}")
print(f"  Phase (s1 s2)^3 = {scalar_phase:.6f}  |phase| = {abs(scalar_phase):.6f}")

# ----- Exhaustive search for best K3 per braid length ------------------------
def best_K3_exhaustive(delta, Lmax):
    gens = generators(delta)
    labels = ['s1', 'S1', 's2', 'S2']
    words = gens.copy()
    word_lbls = [[i] for i in range(4)]
    out = {}
    for L in range(1, Lmax + 1):
        B2 = words @ words
        CB = C_of(words)
        CB2 = C_of(B2)
        K3 = 2 * CB - CB2
        idx = int(np.argmax(K3))
        out[L] = dict(K3=float(K3[idx]), CB=float(CB[idx]), CB2=float(CB2[idx]),
                      word=' '.join(labels[i] for i in word_lbls[idx]),
                      n_words=len(words))
        if L < Lmax:
            words = np.matmul(words[:, None, :, :],
                              gens[None, :, :, :]).reshape(-1, 2, 2)
            word_lbls = [w + [i] for w in word_lbls for i in range(4)]
    return out


print("\n=== Exhaustive LGI search (Ising), sector delta=0 ===")
Lmax = 11
res0 = best_K3_exhaustive(0.0, Lmax)
print(f"{'L':>3} {'#words':>9} {'best K3':>9} {'C(B)':>8} {'C(B^2)':>8}  word")
for L in range(1, Lmax + 1):
    r = res0[L]
    print(f"{L:>3} {r['n_words']:>9} {r['K3']:>9.5f} {r['CB']:>8.4f} "
          f"{r['CB2']:>8.4f}  {r['word']}")

K3_by_L = [res0[L]['K3'] for L in range(1, Lmax + 1)]
first_violation = next((L for L in range(1, Lmax + 1) if res0[L]['K3'] > 1.0 + 1e-9), None)
best_overall = max(range(1, Lmax + 1), key=lambda L: res0[L]['K3'])
print(f"\n  LGI violation (K3>1) first at L = {first_violation}")
print(f"  Best K3 = {res0[best_overall]['K3']:.6f} at L = {best_overall}")
print(f"  Lueders bound 1.5 -> distance = {1.5 - res0[best_overall]['K3']:.6f}")

# ----- Random search to longer lengths ---------------------------------------
print("\n=== Random search L=12..40 (convergence) ===")
gens0 = generators(0.0)
rng = np.random.default_rng(20260525)
conv = {}
for L in range(12, 41, 2):
    nsamp = 400_000
    idx = rng.integers(0, 4, size=(nsamp, L))
    W = np.tile(np.eye(2, dtype=complex), (nsamp, 1, 1))
    for k in range(L):
        W = np.matmul(W, gens0[idx[:, k]])
    K3 = 2 * C_of(W) - C_of(np.matmul(W, W))
    conv[L] = float(np.max(K3))
    print(f"  L={L:>3}  best K3 (random {nsamp:,}) = {conv[L]:.6f}   "
          f"distance to 1.5 = {1.5 - conv[L]:.6f}")

# ----- Sector sweep (delta) ---------------------------------------------------
print("\n=== Sector sweep (Ising): best K3(delta), exhaustive L<=9 ===")
deltas = np.linspace(0, 2 * np.pi, 49)
sector = []
for d in deltas:
    r = best_K3_exhaustive(d, 9)
    sector.append(r[9]['K3'])
sector = np.array(sector)
print(f"  K3(delta) range: [{sector.min():.4f}, {sector.max():.4f}]")
print(f"  spread over sectors = {sector.max()-sector.min():.4f}")
imax = int(np.argmax(sector)); imin = int(np.argmin(sector))
print(f"  Max at delta/pi = {deltas[imax]/np.pi:.3f}  (K3={sector[imax]:.4f})")
print(f"  Min at delta/pi = {deltas[imin]/np.pi:.3f}  (K3={sector[imin]:.4f})")

# ----- Save -------------------------------------------------------------------
results = dict(
    model="Ising",
    sanity=dict(unitarity_err=float(unit_err),
                yang_baxter_err=float(np.max(np.abs(yb))),
                delta2_scalar_err=float(scalar_err),
                delta2_phase_real=float(scalar_phase.real),
                delta2_phase_imag=float(scalar_phase.imag)),
    exhaustive_delta0={str(L): res0[L] for L in res0},
    first_violation_L=first_violation,
    best=dict(L=best_overall, K3=res0[best_overall]['K3']),
    random_convergence={str(L): conv[L] for L in conv},
    sector_sweep=dict(delta_over_pi=list(deltas / np.pi),
                      K3=list(map(float, sector)),
                      spread=float(sector.max() - sector.min())),
)
with open('lgi_ising_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\n-> lgi_ising_results.json written")
