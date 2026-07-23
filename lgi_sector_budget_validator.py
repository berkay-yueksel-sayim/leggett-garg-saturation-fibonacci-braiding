"""
Sector budget validator (Sec. 'Sector dependence and Fourier period').

Validates a single published statement at one fixed sector phase: that
the depression of the exhaustive L <= 9 envelope K3_max(delta) near the
algebraic singularity delta* = 3*pi/5 is a braid-budget effect, not a
sector exclusion.  The weakest NONsingular sectors of the deposited
240-point period scan (lgi_period_and_purestate_results.json) are the
immediate grid neighbors of the singularity, delta/pi = 0.5917 and
0.6083, where the L <= 9 envelope reaches only K3 ~ 1.051.  This script
performs seeded random sampling of longer braid words at the fixed
sector delta/pi = 73/120 = 0.60833... and reports the best K3 found.

This is a validation layer, not a discovery search: it evaluates one
pre-registered sector at longer word lengths, mirroring the random-
search construction already deposited in lgi_fibonacci.py (same
generator set {sigma1, sigma1^-1, sigma2, sigma2^-1}, same correlator
C(U) = 0.5 Re Tr[Z U Z U^dag], same K3 = 2 C(B) - C(B^2), same sample
count 400,000 per length).  As a positive control it reruns the same
code at delta = 0, L = 24, where the deposited random search
(lgi_results.json -> random_convergence) reached 1.49996.

Deterministic: fixed seed below.
"""
import numpy as np
import json

SEED = 20260723
NSAMP = 400_000
DELTA_OVER_PI = 73.0 / 120.0          # 0.608333..., grid neighbor of 3/5
LENGTHS = [12, 24]

phi = (1 + np.sqrt(5)) / 2
F = np.array([[1/phi, 1/np.sqrt(phi)],
              [1/np.sqrt(phi), -1/phi]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def generators(delta=0.0):
    R1 = np.exp(-4j*np.pi/5)
    Rt = np.exp(3j*np.pi/5) * np.exp(1j*delta)
    s1 = np.diag([R1, Rt]).astype(complex)
    s2 = F @ s1 @ F
    return np.stack([s1, s1.conj().T, s2, s2.conj().T])


def C_of(U):
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))
    return 0.5 * np.real(M[:, 0, 0] + M[:, 1, 1])


def best_K3_random(delta, L, nsamp, seed):
    gens = generators(delta)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, 4, size=(nsamp, L))
    W = np.tile(np.eye(2, dtype=complex), (nsamp, 1, 1))
    for k in range(L):
        W = np.matmul(W, gens[idx[:, k]])
    K3 = 2*C_of(W) - C_of(np.matmul(W, W))
    return float(np.max(K3))


def main():
    delta = DELTA_OVER_PI * np.pi
    print("=== Sector budget validator ===")
    print(f"  sector delta/pi = {DELTA_OVER_PI:.6f} (grid neighbor of the "
          f"singularity 3/5), seed = {SEED}, samples per length = {NSAMP:,}")

    results = {}
    for L in LENGTHS:
        k3 = best_K3_random(delta, L, NSAMP, SEED + L)
        results[str(L)] = k3
        print(f"  L = {L:>2}: best K3 (random) = {k3:.7f}   "
              f"gap to Lueders 1.5 = {1.5 - k3:.7f}")

    # deposited L <= 9 exhaustive envelope value at this grid point,
    # read from the deposited period scan for reference
    scan = json.load(open('lgi_period_and_purestate_results.json'))["period_scan"]
    i = int(np.argmin(np.abs(np.array(scan["delta_over_pi"]) - DELTA_OVER_PI)))
    env9 = scan["K3"][i]
    print(f"  deposited L<=9 envelope at this sector: K3 = {env9:.6f}")

    # positive control: same code at delta = 0, L = 24.  The deposited
    # random search (lgi_results.json -> random_convergence; different
    # seed, same construction) reached 1.4987 at L = 24 and 1.499964 as
    # its series maximum; a fresh seed must land in the same 1.498-1.500
    # class (maxima of 400,000 draws scatter at the 1e-3 level).
    ctrl = best_K3_random(0.0, 24, NSAMP, SEED)
    rc = json.load(open('lgi_results.json'))["random_convergence"]
    dep24, dep_max = rc["24"], max(rc.values())
    ctrl_ok = bool(abs(ctrl - dep_max) < 2e-3 or abs(ctrl - dep24) < 2e-3)
    print(f"  positive control (delta=0, L=24): {ctrl:.6f} vs deposited "
          f"L=24 {dep24:.6f} / series max {dep_max:.6f} -> "
          f"class match = {ctrl_ok}")

    out = dict(
        meta=dict(script="lgi_sector_budget_validator.py", seed=SEED,
                  nsamp=NSAMP, delta_over_pi=DELTA_OVER_PI,
                  note="Seeded random-search validation of the braid-budget "
                       "statement at the weakest nonsingular sector; not a "
                       "landscape sweep."),
        best_K3_random=results,
        deposited_L9_envelope_at_sector=env9,
        positive_control=dict(delta0_L24=ctrl, deposited_delta0_L24=dep24,
                              deposited_series_max=dep_max,
                              class_match=ctrl_ok),
    )
    with open('lgi_sector_budget_validator.json', 'w') as f:
        json.dump(out, f, indent=2)
    print("-> lgi_sector_budget_validator.json written")


if __name__ == "__main__":
    main()
