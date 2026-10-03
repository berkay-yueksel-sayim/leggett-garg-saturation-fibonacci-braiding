r"""
Braid relation and mirror symmetry of the deformed Fibonacci generators.

The letter deforms the Fibonacci braid generators by a phase delta,

    sigma_1(delta) = diag(R_1, R_tau exp(i delta)),   sigma_2 = F sigma_1 F,

with R_1 = exp(-4 pi i/5), R_tau = exp(+3 pi i/5) (lgi_fibonacci.py
convention). This script reproduces the two statements of the letter that
rest on this family (Secs. II B and III D):

  (1) Braid relation.  sigma_1 sigma_2 sigma_1 = sigma_2 sigma_1 sigma_2
      holds in [0, 2 pi) only at delta = 0, 3 pi/5 and 6 pi/5 -- exactly,
      not merely up to a global phase -- and at no other delta, not even
      up to a global phase. At 3 pi/5 both generators are the scalar
      exp(-4 pi i/5) times the identity; at 6 pi/5 they are the complex
      conjugates of the delta = 0 generators times a global phase.
  (2) Mirror symmetry.  K3(delta; w) = K3(6 pi/5 - delta; w) for every
      word w, hence also for the envelope K3_max(delta); the mirror axes
      are delta = 3 pi/5 and 8 pi/5.

Norm: Frobenius. "Up to a global phase" means min over theta of
|| A - exp(i theta) B ||, evaluated with theta = arg Tr(B^dagger A)
(no square-root cancellation, so no precision floor near zero).

Checks:
  * positive control at delta = 0;
  * 240-point grid delta_k = 2 pi k / 240 (the grid of the letter);
  * fine scan with 200,000 points on [0, 2 pi): local minima of the
    distance (a zero between grid points would show up as a minimum);
  * the scalar point and the conjugate point, with a negative control
    one grid step away;
  * envelope mirror symmetry on the deposited envelope
    (lgi_period_and_purestate_results.json, period_scan.K3), both axes,
    with a negative control about a wrong axis;
  * word level: all 340 words of length 1..4 at all 240 grid points.

Deterministic: no random numbers are drawn. Only numpy is used. The result
file is written next to this script.
"""
import os
import json
import time
import numpy as np

_HIER = os.path.dirname(os.path.abspath(__file__))

N_GRID = 240
N_FINE = 200000
TOL_EXACT = 1e-12
WORD_LMAX = 4

phi = (1 + np.sqrt(5)) / 2
F = np.array([[1/phi, 1/np.sqrt(phi)],
              [1/np.sqrt(phi), -1/phi]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
R1 = np.exp(-4j * np.pi / 5)
RT = np.exp(3j * np.pi / 5)


def gens(deltas):
    """sigma_1, sigma_2 for an array of deltas, shape (n, 2, 2) each."""
    d = np.atleast_1d(np.asarray(deltas, dtype=float))
    s1 = np.zeros((d.size, 2, 2), dtype=complex)
    s1[:, 0, 0] = R1
    s1[:, 1, 1] = RT * np.exp(1j * d)
    s2 = F @ s1 @ F
    return s1, s2


def fro(A):
    return np.sqrt(np.sum(np.abs(A) ** 2, axis=(-2, -1)))


def braid_distances(deltas):
    """exact and up-to-global-phase distance of s1 s2 s1 and s2 s1 s2."""
    s1, s2 = gens(deltas)
    A = s1 @ s2 @ s1
    B = s2 @ s1 @ s2
    exact = fro(A - B)
    ov = np.einsum('nij,nij->n', B.conj(), A)          # Tr(B^dagger A)
    ph = np.where(np.abs(ov) > 0, ov / np.abs(ov), 1.0)
    up_to_phase = fro(A - ph[:, None, None] * B)
    return exact, up_to_phase


def local_minima(y):
    """indices i with y[i] < y[i-1] and y[i] < y[i+1] (circular)."""
    return np.where((y < np.roll(y, 1)) & (y < np.roll(y, -1)))[0]


def C_mixed(U):
    return 0.5 * np.real(np.einsum('ij,njk,kl,nli->n', Z, U, Z, U.conj().transpose(0, 2, 1)))


t0 = time.time()
out = dict(description="Braid relation and mirror symmetry of the deformed "
                       "Fibonacci generators (Secs. II B, III D of the letter).",
           norm="Frobenius", deterministic="no random numbers")

# (1a) positive control
e0, p0 = braid_distances([0.0])
out['positive_control_delta0'] = dict(exact=float(e0[0]), up_to_phase=float(p0[0]))
print("positive control delta = 0: exact %.2e, up to phase %.2e" % (e0[0], p0[0]))

# (1b) the 240-point grid
kg = np.arange(N_GRID)
dg = 2 * np.pi * kg / N_GRID
eg, pg = braid_distances(dg)
zero = np.where(eg <= TOL_EXACT)[0]
rest = np.setdiff1d(kg, zero)
out['grid'] = dict(n=N_GRID, tolerance=TOL_EXACT,
                   exact_zero_k=[int(k) for k in zero],
                   exact_zero_delta_over_pi=[float(dg[k] / np.pi) for k in zero],
                   exact_at_zeros=[float(eg[k]) for k in zero],
                   n_rest=int(rest.size),
                   rest_up_to_phase_min=float(pg[rest].min()),
                   rest_up_to_phase_max=float(pg[rest].max()),
                   rest_exact_min=float(eg[rest].min()))
print("grid %d: exact zeros at k = %s (delta/pi %s)" % (N_GRID, out['grid']['exact_zero_k'],
                                                          [round(x, 4) for x in out['grid']['exact_zero_delta_over_pi']]))
print("   other %d points, up to phase: min %.4f, max %.4f" % (rest.size, pg[rest].min(), pg[rest].max()))

# (1c) fine scan
df = 2 * np.pi * np.arange(N_FINE) / N_FINE
ef, pf = braid_distances(df)
me, mp = local_minima(ef), local_minima(pf)
out['fine_scan'] = dict(n=N_FINE,
                        local_minima_exact=[dict(delta_over_pi=float(df[i] / np.pi), value=float(ef[i])) for i in me],
                        local_minima_up_to_phase=[dict(delta_over_pi=float(df[i] / np.pi), value=float(pf[i])) for i in mp])
print("fine scan %d: local minima (exact) %s" % (N_FINE, ["%.6f: %.2e" % (df[i] / np.pi, ef[i]) for i in me]))
print("   local minima (up to phase) %s" % ["%.6f: %.2e" % (df[i] / np.pi, pf[i]) for i in mp])

# (1d) the scalar point
s1, s2 = gens([3 * np.pi / 5])
c = R1
out['scalar_point'] = dict(delta_over_pi=0.6, c=[float(c.real), float(c.imag)],
                           dev_sigma1=float(fro(s1 - c * np.eye(2))[0]),
                           dev_sigma2=float(fro(s2 - c * np.eye(2))[0]))
print("delta = 3pi/5: ||s1 - c I|| %.2e, ||s2 - c I|| %.2e, c = exp(-4 pi i/5)"
      % (out['scalar_point']['dev_sigma1'], out['scalar_point']['dev_sigma2']))


# (1e) the conjugate point, and one grid step away
def conj_fit(delta):
    a1, a2 = gens([delta])
    b1, b2 = gens([0.0])
    b1, b2 = b1.conj(), b2.conj()
    ov = np.einsum('nij,nij->n', b1.conj(), a1)[0]
    lam = ov / abs(ov)
    return lam, float(fro(a1 - lam * b1)[0]), float(fro(a2 - lam * b2)[0])


lam, r1, r2 = conj_fit(6 * np.pi / 5)
lam_n, r1n, r2n = conj_fit(2 * np.pi * 145 / N_GRID)
out['conjugate_point'] = dict(delta_over_pi=1.2, arg_lambda_over_pi=float(np.angle(lam) / np.pi),
                              residual_sigma1=r1, residual_sigma2=r2,
                              negative_control=dict(delta_over_pi=float(2 * 145 / N_GRID),
                                                    residual_sigma1=r1n, residual_sigma2=r2n))
print("delta = 6pi/5: sigma_i = lambda conj(sigma_i(0)), arg lambda/pi = %.4f, residuals %.2e / %.2e"
      % (np.angle(lam) / np.pi, r1, r2))
print("   negative control k = 145: residuals %.4f / %.4f" % (r1n, r2n))

# (2a) envelope mirror symmetry on the deposited envelope
with open(os.path.join(_HIER, 'lgi_period_and_purestate_results.json')) as f:
    K = np.array(json.load(f)['period_scan']['K3'])
ax1 = float(np.max(np.abs(K - K[(144 - kg) % N_GRID])))
ax2 = float(np.max(np.abs(K[(192 + kg) % N_GRID] - K[(192 - kg) % N_GRID])))
neg = float(np.max(np.abs(K - K[(120 - kg) % N_GRID])))
out['envelope_mirror'] = dict(source="lgi_period_and_purestate_results.json: period_scan.K3",
                              axis_3pi_over_5=ax1, axis_8pi_over_5=ax2,
                              negative_control_axis_pi_over_2=neg,
                              K3_at_0=float(K[0]), K3_at_6pi_over_5=float(K[144]),
                              K3_at_0p5917pi=float(K[71]), K3_at_0p6083pi=float(K[73]))
print("envelope mirror: axis 3pi/5 %.2e, axis 8pi/5 %.2e, wrong axis %.4f" % (ax1, ax2, neg))

# (2b) word level: all words of length 1..WORD_LMAX at all grid points
s1g, s2g = gens(dg)
alph = np.stack([s1g, s1g.conj().transpose(0, 2, 1), s2g, s2g.conj().transpose(0, 2, 1)])
words = [()]
worst, n_words = 0.0, 0
mirror_idx = (144 - kg) % N_GRID
for L in range(1, WORD_LMAX + 1):
    words = [w + (a,) for w in words for a in range(4)]
    for w in words:
        B = np.broadcast_to(np.eye(2, dtype=complex), (N_GRID, 2, 2)).copy()
        for a in w:
            B = B @ alph[a]
        K3w = 2 * C_mixed(B) - C_mixed(B @ B)
        worst = max(worst, float(np.max(np.abs(K3w - K3w[mirror_idx]))))
        n_words += 1
out['word_level_mirror'] = dict(n_words=n_words, lengths="1..%d" % WORD_LMAX,
                                n_delta=N_GRID, max_abs_dev=worst)
print("word level: %d words x %d points, max |K3(d;w) - K3(6pi/5-d;w)| = %.2e" % (n_words, N_GRID, worst))

# the runtime is printed, not stored: the result file stays byte-reproducible
print("runtime %.1f s" % (time.time() - t0))
with open(os.path.join(_HIER, 'lgi_braid_relation_check_results.json'), 'w') as f:
    json.dump(out, f, indent=2)
print("\n-> lgi_braid_relation_check_results.json written")
