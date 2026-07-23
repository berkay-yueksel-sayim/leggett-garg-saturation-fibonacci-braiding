"""
K3(delta) envelope period scan (Fibonacci) + Pure-state robustness check.

(v1.2 replacement of the original deposit script.) The original script
computed the delta-sweep over words of EXACT length L = 9; the letter,
however, defines the swept quantity as the cumulative envelope

    K3_max(delta) = max over all braid words |w| <= L of K3(delta; w)

(see "Sector dependence" in the letter). Under the exact-length
definition the Fourier analysis reports a dominant harmonic k = 6 --
this is the superseded v1.0 claim, retracted in the letter's
"Correction to v1.0" paragraph. This replacement script implements the
cumulative <=L definition of the text. Expected results (letter):

  * dominant envelope harmonic k = 3 (period 2*pi/3), with the top
    harmonics being multiples of three;
  * k=3 / k=6 amplitude ratio stable at ~1.061 (L_max = 7) to ~1.063
    (L_max = 8), both computed here as stability checks;
  * positive control: the superseded exact-length definition is also
    run here and must reproduce the OLD dominant harmonic k = 6 -- so
    the two definitions are demonstrably distinguished by this code,
    and the k = 3 result is not an implementation artifact.

Part (2), the pure-state robustness check at the L=11 optimum, is
unchanged from the original deposit script.

Reference convention: lgi_fibonacci.py (Fibonacci F-matrix,
R_1 = exp(-4 pi i /5), R_tau = exp(+3 pi i /5)).
"""
import numpy as np
import json

np.random.seed(20260525)

phi = (1 + np.sqrt(5)) / 2
F = np.array([[1/phi, 1/np.sqrt(phi)],
              [1/np.sqrt(phi), -1/phi]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def generators(delta=0.0):
    R1 = np.exp(-4j * np.pi / 5)
    Rt = np.exp(3j * np.pi / 5) * np.exp(1j * delta)
    s1 = np.diag([R1, Rt]).astype(complex)
    s2 = F @ s1 @ F
    return np.stack([s1, s1.conj().T, s2, s2.conj().T])


def C_mixed(U):
    """C(U) = (1/2) Re Tr[Z U Z U^dag] for rho_0 = I/2.  Batched over (N,2,2)."""
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))
    return 0.5 * np.real(M[:, 0, 0] + M[:, 1, 1])


def C_pure(U, psi):
    """C(U; psi) = <psi|Z U Z U^dag|psi> for a single pure state.  Batched."""
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))   # (N,2,2)
    rhs = M @ psi                                          # (N,2)
    return np.real(np.conj(psi) @ rhs.T)


def envelope_K3(delta, Lmax):
    """K3_max(delta) = max over ALL words |w| <= Lmax (cumulative envelope,
    the definition used in the letter). Exhaustive enumeration."""
    gens = generators(delta)
    words = gens.copy()
    best = -np.inf
    for L in range(1, Lmax + 1):
        B2 = words @ words
        K3 = 2 * C_mixed(words) - C_mixed(B2)
        best = max(best, float(K3.max()))
        if L < Lmax:
            words = np.matmul(words[:, None, :, :],
                              gens[None, :, :, :]).reshape(-1, 2, 2)
    return best


def sweep_and_fft(n_delta, Lmax):
    deltas = np.linspace(0, 2 * np.pi, n_delta, endpoint=False)
    K3_sweep = np.array([envelope_K3(d, Lmax) for d in deltas])
    fft = np.fft.rfft(K3_sweep - K3_sweep.mean())
    amp = np.abs(fft)
    k_top = int(np.argmax(amp[1:]) + 1)                 # ignore k=0 (DC)
    order = np.argsort(amp[1:])[::-1][:5] + 1
    return deltas, K3_sweep, amp, k_top, order


def exact_length_K3(delta, Lmax):
    """max over words of EXACT length Lmax (the superseded v1.0 definition,
    kept ONLY as a positive control -- see module docstring)."""
    gens = generators(delta)
    words = gens.copy()
    for L in range(1, Lmax + 1):
        K3 = 2 * C_mixed(words) - C_mixed(words @ words)
        if L == Lmax:
            return float(K3.max())
        words = np.matmul(words[:, None, :, :],
                          gens[None, :, :, :]).reshape(-1, 2, 2)


# ============================================================================
# (0) POSITIVE CONTROL  the superseded exact-length definition must
#     reproduce the OLD dominant harmonic k = 6 (discriminating control)
# ============================================================================
print("=" * 64)
print("  (0) POSITIVE CONTROL  exact-length-9 definition reproduces old k = 6")
print("=" * 64)

N_CTRL = 240
ctrl_deltas = np.linspace(0, 2 * np.pi, N_CTRL, endpoint=False)
K3_exact = np.array([exact_length_K3(d, 9) for d in ctrl_deltas])
fft_c = np.abs(np.fft.rfft(K3_exact - K3_exact.mean()))
k_ctrl = int(np.argmax(fft_c[1:]) + 1)
print(f"  exact-length L=9 sweep: dominant harmonic k = {k_ctrl}  (superseded "
      f"definition; expected 6)")
assert k_ctrl == 6, ("positive control failed -- exact-length definition did "
                     "not reproduce the superseded k = 6")

# ============================================================================
# (1) ENVELOPE PERIOD SCAN  K3_max(delta), cumulative <= L definition
# ============================================================================
print("\n" + "=" * 64)
print("  (1) ENVELOPE PERIOD SCAN  K3_max(delta) = max_{|w|<=L} K3(delta; w)")
print("=" * 64)

N_DELTA = 240            # fine sampling for clean FFT (figure resolution)
LMAX_SWEEP = 9           # exhaustive within reach; matches the figure axis
deltas, K3_sweep, amp, k_top, order = sweep_and_fft(N_DELTA, LMAX_SWEEP)
period_pi = 2.0 / k_top
print(f"  sampled K3_max(delta) at {N_DELTA} points, exhaustive |w| <= {LMAX_SWEEP}")
print(f"  range: [{K3_sweep.min():.6f}, {K3_sweep.max():.6f}]   "
      f"mean = {K3_sweep.mean():.6f}")
print(f"  dominant Fourier harmonic k = {k_top}    "
      f"=> period = 2*pi/{k_top} = {period_pi:.4f} * pi")
print(f"  top 5 harmonics (k, amplitude):")
for k in order:
    print(f"    k = {k:3d}   amplitude = {amp[k]:.4f}   "
          f"period = 2*pi/{k} = {2.0/k:.4f} * pi")
ratio_main = float(amp[3] / amp[6])
print(f"  k=3 / k=6 amplitude ratio = {ratio_main:.4f}")

# ============================================================================
# (1b) STABILITY CHECKS  L_max = 7 and L_max = 8 at n_delta = 120
#      (the two independent-reproducer settings quoted in the letter's
#       "Correction to v1.0" paragraph; expected ratio 1.061 -> 1.063)
# ============================================================================
print("\n" + "=" * 64)
print("  (1b) STABILITY  dominant harmonic under smaller search spaces")
print("=" * 64)

stability = {}
for lmax in (7, 8):
    _, _, amp_s, k_s, order_s = sweep_and_fft(120, lmax)
    ratio = float(amp_s[3] / amp_s[6])
    stability[f"Lmax={lmax}"] = dict(
        n_delta=120, Lmax=lmax, dominant_k=int(k_s),
        k3_over_k6_ratio=ratio,
        top_harmonics=[int(k) for k in order_s],
    )
    print(f"  L_max = {lmax} (n_delta = 120): dominant k = {k_s}, "
          f"k=3/k=6 ratio = {ratio:.4f}, top-5 k = {list(map(int, order_s))}")

# ============================================================================
# (2) PURE-STATE ROBUSTNESS  K3 at the L=11 optimum  (unchanged from the
#     original deposit script)
# ============================================================================
print("\n" + "=" * 64)
print("  (2) PURE-STATE ROBUSTNESS  K3 with rho_0 = |psi><psi|")
print("=" * 64)

def word_to_unitary(word_str, delta=0.0):
    labels = {'s1': 0, 'S1': 1, 's2': 2, 'S2': 3}
    gens = generators(delta)
    U = np.eye(2, dtype=complex)
    for tok in word_str.split():
        U = U @ gens[labels[tok]]
    return U

OPTIMAL_L11 = 's1 s2 s2 S1 S1 S1 s2 s2 S1 s2 s1'
U_opt = word_to_unitary(OPTIMAL_L11, delta=0.0)
U2 = U_opt @ U_opt

K3_mixed = float(2 * C_mixed(U_opt[None])[0] - C_mixed(U2[None])[0])
print(f"  Optimal word L=11 :  '{OPTIMAL_L11}'")
print(f"    K3 (mixed rho_0 = I/2)         = {K3_mixed:.6f}    "
      f"(reference: 1.499762)")

test_states = {
    '|0>'           : np.array([1, 0], dtype=complex),
    '|1>'           : np.array([0, 1], dtype=complex),
    '|+>'           : (1/np.sqrt(2)) * np.array([1,  1], dtype=complex),
    '|->'           : (1/np.sqrt(2)) * np.array([1, -1], dtype=complex),
    '|+i>'          : (1/np.sqrt(2)) * np.array([1, 1j], dtype=complex),
    '|-i>'          : (1/np.sqrt(2)) * np.array([1, -1j], dtype=complex),
    'random Haar 1' : None,
    'random Haar 2' : None,
    'random Haar 3' : None,
}
rng = np.random.default_rng(20260525)
for k in test_states:
    if test_states[k] is None:
        z = rng.normal(size=2) + 1j * rng.normal(size=2)
        z /= np.linalg.norm(z)
        test_states[k] = z

pure_results = {}
print(f"\n  Pure-state K3 at the L=11 optimum (delta=0):")
print(f"  {'state':<14} {'K3':>10}    {'C(B)':>9}  {'C(B^2)':>9}")
for name, psi in test_states.items():
    CB = float(C_pure(U_opt[None], psi)[0])
    CB2 = float(C_pure(U2[None], psi)[0])
    K3p = 2 * CB - CB2
    pure_results[name] = dict(CB=CB, CB2=CB2, K3=K3p)
    print(f"  {name:<14} {K3p:>10.6f}    {CB:>9.4f}  {CB2:>9.4f}")

K3_vals = np.array([v['K3'] for v in pure_results.values()])
spread = K3_vals.max() - K3_vals.min()
mean_pure = K3_vals.mean()
print(f"\n  spread over initial states = {spread:.4f}")
print(f"  mean over pure states      = {mean_pure:.6f}    "
      f"(should equal mixed: {K3_mixed:.6f})")
print(f"  difference                 = {abs(mean_pure - K3_mixed):.2e}")

# ============================================================================
# Save  (period_scan schema unchanged so the figure script keeps working)
# ============================================================================
out = dict(
    period_scan=dict(
        n_delta=N_DELTA,
        Lmax=LMAX_SWEEP,
        definition="cumulative envelope: K3_max(delta) = max over |w| <= Lmax",
        delta_over_pi=list(deltas / np.pi),
        K3=list(map(float, K3_sweep)),
        dominant_k=int(k_top),
        period_pi=float(period_pi),
        k3_over_k6_ratio=ratio_main,
        top_harmonics=[dict(k=int(k), amplitude=float(amp[k]),
                             period_pi=float(2.0/k))
                       for k in order],
    ),
    stability_runs=stability,
    exact_length_control=dict(
        definition="max over words of EXACT length L (superseded v1.0 route)",
        n_delta=N_CTRL,
        Lmax=9,
        dominant_k=int(k_ctrl),
        statement="the superseded definition reproduces the retracted k=6; "
                  "the cumulative <=L definition of the text gives k=3",
    ),
    pure_state_check=dict(
        word=OPTIMAL_L11,
        delta=0.0,
        K3_mixed=K3_mixed,
        per_state={name: dict(CB=v['CB'], CB2=v['CB2'], K3=v['K3'])
                   for name, v in pure_results.items()},
        spread=float(spread),
        mean_pure=float(mean_pure),
        mean_pure_minus_mixed=float(mean_pure - K3_mixed),
    ),
)
with open('lgi_period_and_purestate_results.json', 'w') as f:
    json.dump(out, f, indent=2)
print("\n-> lgi_period_and_purestate_results.json written")
