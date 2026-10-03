"""
Envelope word audit for the K3(delta) sector-phase sweep (Fibonacci).

Question: the letter states that fixed-word traces K3(delta; w) are generically
only 2*pi-periodic, that the envelope's Fourier structure comes from
optimal-word reshuffling across the sweep, and that the k = 3 harmonic is a
property of the maximization, not of any single word. This script measures
the four quantities behind those statements, on the main run of the letter:

  engine      : identical to lgi_period_and_purestate.py (Fibonacci F-matrix,
                R_1 = exp(-4 pi i/5), R_tau = exp(+3 pi i/5) exp(i delta),
                rho_0 = I/2, K3 = 2 C(B) - C(B^2))
  words       : ALL braid words over {s1, s1^-1, s2, s2^-1} with |w| <= 9
                (349,524 words), same enumeration order as the engine
  grid        : n_delta = 240, delta in [0, 2 pi), endpoint excluded
  tolerance   : a word "reaches the envelope" at a point if
                K3(delta; w) >= K3_max(delta) - 1e-12

  (a) Does ONE word reach the envelope at all 240 points?
  (b) W* = all words that reach the envelope at >= 1 point. For each: the
      Fourier spectrum of K3(delta; w) on the same 240 points (k = 0
      excluded) -- which component is the largest?
  (c) For W*: smallest period (invariance under delta -> delta + 2 pi/m,
      m = 2..6, i.e. a shift by 240/m grid points, max deviation <= 1e-12).
  (d) A uniform random sample of 1000 of all 349,524 words (fixed seed):
      fraction whose smallest period is 2 pi. Constant curves are listed
      separately.

Degenerate grid points: at delta = 3 pi/5 (grid point 72) R_tau exp(i delta)
equals R_1, both generators are multiples of the identity and EVERY word has
K3 = 1 -- so every word "reaches" the envelope there and W* as defined above
contains all words. The script detects such points (all words equal within the
tolerance) and reports (b), (c) twice: for W* as defined (all 240 points) and
for W*' (words reaching the envelope at >= 1 NONdegenerate point). Per-word
rows are written for W*' only (up to ROW_LIMIT words); for W* the counts.

Positive control: the envelope assembled here must reproduce the deposited
envelope in lgi_period_and_purestate_results.json (period_scan.K3) at all 240
points. Cross-check: a sample of words is decoded to its letter string and
re-evaluated by an independent sequential product.

The result file is written next to this script.
"""
import os
import json
import time
import numpy as np

_HIER = os.path.dirname(os.path.abspath(__file__))

SEED = 20260928
N_DELTA = 240
LMAX = 9
TOL = 1e-12
N_SAMPLE = 1000
M_SHIFTS = (2, 3, 4, 5, 6)
ROW_LIMIT = 5000
CHUNK = 20000

phi = (1 + np.sqrt(5)) / 2
F = np.array([[1/phi, 1/np.sqrt(phi)],
              [1/np.sqrt(phi), -1/phi]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
LABELS = ('s1', 'S1', 's2', 'S2')      # generator order of the engine


def generators(delta=0.0):
    R1 = np.exp(-4j * np.pi / 5)
    Rt = np.exp(3j * np.pi / 5) * np.exp(1j * delta)
    s1 = np.diag([R1, Rt]).astype(complex)
    s2 = F @ s1 @ F
    return np.stack([s1, s1.conj().T, s2, s2.conj().T])


def C_mixed(U):
    M = Z @ U @ Z @ np.conj(np.transpose(U, (0, 2, 1)))
    return 0.5 * np.real(M[:, 0, 0] + M[:, 1, 1])


def all_K3(delta):
    """K3(delta; w) for all words |w| <= LMAX, levels concatenated in the
    engine's enumeration order (length 1 first)."""
    gens = generators(delta)
    words = gens.copy()
    parts = []
    for L in range(1, LMAX + 1):
        parts.append(2 * C_mixed(words) - C_mixed(words @ words))
        if L < LMAX:
            words = np.matmul(words[:, None, :, :],
                              gens[None, :, :, :]).reshape(-1, 2, 2)
    return np.concatenate(parts)


OFFSETS = np.cumsum([0] + [4 ** L for L in range(1, LMAX + 1)])
N_WORDS = int(OFFSETS[-1])


def decode(idx):
    """global index -> word string (first letter = most significant digit)."""
    L = int(np.searchsorted(OFFSETS, idx, side='right'))
    i = int(idx - OFFSETS[L - 1])
    digits = []
    for _ in range(L):
        digits.append(i % 4)
        i //= 4
    return ' '.join(LABELS[d] for d in reversed(digits))


def K3_of_word(word, delta):
    """independent sequential evaluation of one word (cross-check)."""
    gens = generators(delta)
    U = np.eye(2, dtype=complex)
    for tok in word.split():
        U = U @ gens[LABELS.index(tok)]
    return float(2 * C_mixed(U[None])[0] - C_mixed((U @ U)[None])[0])


def analyse(curves):
    """per row: constant?, largest Fourier component k (k=0 excluded),
    smallest period as m (period 2 pi/m; m = 1 means 2 pi), top-3 k."""
    const = np.ptp(curves, axis=1) <= TOL
    amp = np.abs(np.fft.rfft(curves - curves.mean(axis=1, keepdims=True), axis=1))
    kmax = np.argmax(amp[:, 1:], axis=1) + 1
    top3 = np.argsort(amp[:, 1:], axis=1)[:, ::-1][:, :3] + 1
    m_best = np.ones(curves.shape[0], dtype=int)
    devs = {}
    for m in M_SHIFTS:
        dev = np.max(np.abs(curves - np.roll(curves, N_DELTA // m, axis=1)), axis=1)
        devs[m] = dev
        m_best = np.where(dev <= TOL, np.maximum(m_best, m), m_best)
    return const, kmax, top3, amp, m_best, devs


def zaehle(werte, maske=None):
    w = werte if maske is None else werte[maske]
    u, c = np.unique(w, return_counts=True)
    return {str(int(a)): int(b) for a, b in zip(u, c)}


t0 = time.time()
deltas = np.linspace(0, 2 * np.pi, N_DELTA, endpoint=False)
print("=" * 72)
print("  ENVELOPE WORD AUDIT  |w| <= %d (%d words), n_delta = %d, tol = %.0e"
      % (LMAX, N_WORDS, N_DELTA, TOL))
print("=" * 72)

# ---------------------------------------------------------------- pass 1: envelope, hits, degenerate points
envelope = np.empty(N_DELTA)
degenerate = np.zeros(N_DELTA, dtype=bool)
hits_all = np.zeros(N_WORDS, dtype=np.int32)
hits_nd = np.zeros(N_WORDS, dtype=np.int32)
n_tied = np.zeros(N_DELTA, dtype=np.int64)
for j, d in enumerate(deltas):
    k3 = all_K3(d)
    envelope[j] = k3.max()
    degenerate[j] = np.ptp(k3) <= TOL
    top = k3 >= envelope[j] - TOL
    n_tied[j] = int(top.sum())
    hits_all[top] += 1
    if not degenerate[j]:
        hits_nd[top] += 1
deg_idx = [int(j) for j in np.flatnonzero(degenerate)]
print("  pass 1 done (%.0f s); degenerate grid points (all words equal): %s  delta/pi = %s"
      % (time.time() - t0, deg_idx, [round(float(deltas[j] / np.pi), 12) for j in deg_idx]))

# positive control against the deposited envelope
with open(os.path.join(_HIER, 'lgi_period_and_purestate_results.json')) as fh:
    dep = json.load(fh)['period_scan']
assert dep['n_delta'] == N_DELTA and dep['Lmax'] == LMAX
pc_dev = float(np.max(np.abs(envelope - np.array(dep['K3']))))
pc_ok = pc_dev <= TOL
print("  POSITIVE CONTROL  max |envelope - deposited envelope| = %.2e  -> %s"
      % (pc_dev, "OK" if pc_ok else "FAILED"))
assert pc_ok, "engine does not reproduce the deposited envelope"

n_all = int(np.sum(hits_all == N_DELTA))
n_all_nd = int(np.sum(hits_nd == N_DELTA - len(deg_idx)))
print("  (a) words on the envelope at ALL %d points: %d (best coverage %d points);"
      " at all %d nondegenerate points: %d (best %d)"
      % (N_DELTA, n_all, int(hits_all.max()), N_DELTA - len(deg_idx), n_all_nd, int(hits_nd.max())))
print("      words tied on the envelope per point: min %d, median %d, max %d (nondegenerate points)"
      % (int(n_tied[~degenerate].min()), int(np.median(n_tied[~degenerate])), int(n_tied[~degenerate].max())))

rng = np.random.default_rng(SEED)
sample = np.sort(rng.choice(N_WORDS, size=N_SAMPLE, replace=False))

# ---------------------------------------------------------------- pass 2: all curves (needed for W* = all words)
curves = np.empty((N_WORDS, N_DELTA))
for j, d in enumerate(deltas):
    curves[:, j] = all_K3(d)
print("  pass 2 done (%.0f s)" % (time.time() - t0))

const = np.empty(N_WORDS, dtype=bool)
kmax = np.empty(N_WORDS, dtype=int)
m_best = np.empty(N_WORDS, dtype=int)
top3 = np.empty((N_WORDS, 3), dtype=int)
amp3 = np.empty(N_WORDS)
for a in range(0, N_WORDS, CHUNK):
    b = min(a + CHUNK, N_WORDS)
    c_, k_, t_, amp_, m_, _ = analyse(curves[a:b])
    const[a:b], kmax[a:b], top3[a:b], m_best[a:b] = c_, k_, t_, m_
    amp3[a:b] = amp_[:, 3]
kmax_eff = np.where(const, 0, kmax)          # 0 = constant curve (no Fourier component)
print("  analysis done (%.0f s)" % (time.time() - t0))

# cross-check: decode and re-evaluate independently
xc_idx = list(np.flatnonzero(hits_nd > 0)[:300]) + list(sample[:50])
xc_dev = 0.0
for idx in xc_idx:
    w = decode(idx)
    for j in (0, 57, 131, 239):
        xc_dev = max(xc_dev, abs(K3_of_word(w, deltas[j]) - curves[idx, j]))
print("  CROSS-CHECK %d decoded words re-evaluated independently: max deviation %.2e" % (len(xc_idx), xc_dev))

# ---------------------------------------------------------------- (b), (c) -- W* as defined, and W*'
res_bc = {}
for name, mask in (("W_star_as_defined", hits_all > 0), ("W_star_nondegenerate", hits_nd > 0)):
    kc = zaehle(kmax_eff, mask)
    pc = zaehle(m_best, mask & ~const)
    res_bc[name] = dict(size=int(mask.sum()),
                        n_constant=int((mask & const).sum()),
                        largest_component_counts_k=kc,
                        n_with_k3_largest=int((mask & ~const & (kmax == 3)).sum()),
                        smallest_period_counts_m_nonconstant=pc,
                        n_shorter_than_2pi=int((mask & ~const & (m_best > 1)).sum()))
    print("  (b)(c) %-22s |W|=%d const=%d  largest-k counts %s  k3 largest %d | period m counts %s  shorter than 2pi %d"
          % (name, res_bc[name]['size'], res_bc[name]['n_constant'], kc, res_bc[name]['n_with_k3_largest'],
             pc, res_bc[name]['n_shorter_than_2pi']))

Wnd = np.flatnonzero(hits_nd > 0)
rows = []
for idx in Wnd[:ROW_LIMIT]:
    rows.append(dict(word=decode(idx), points_on_envelope_nondegenerate=int(hits_nd[idx]),
                     constant=bool(const[idx]),
                     largest_component_k=None if const[idx] else int(kmax[idx]),
                     top3_k=[int(k) for k in top3[idx]],
                     smallest_period_over_pi=2.0 / int(m_best[idx])))

# ---------------------------------------------------------------- (d)
sc = const[sample]
sm = m_best[sample]
n_const = int(sc.sum())
n_2pi = int((~sc & (sm == 1)).sum())
frac_all = n_2pi / N_SAMPLE
frac_nc = n_2pi / (N_SAMPLE - n_const)
s_per = zaehle(sm, ~sc)
exc = [dict(word=decode(i), constant=bool(const[i]), smallest_period_over_pi=2.0 / int(m_best[i]))
       for i in sample if const[i] or m_best[i] > 1]
print("  (d) sample %d (seed %d): period 2pi %d, constant %d, other periods %s"
      % (N_SAMPLE, SEED, n_2pi, n_const, {k: v for k, v in s_per.items() if k != '1'}))
print("      fraction period 2pi: %.4f of all 1000, %.4f of the %d nonconstant"
      % (frac_all, frac_nc, N_SAMPLE - n_const))
# same over ALL words (not asked; context for (d))
n_c_all = int(const.sum())
print("      (context, all %d words: constant %d, period 2pi %d = %.4f of nonconstant)"
      % (N_WORDS, n_c_all, int((~const & (m_best == 1)).sum()),
         int((~const & (m_best == 1)).sum()) / (N_WORDS - n_c_all)))

# ---------------------------------------------------------------- rule (fixed in advance by the order)
bc = res_bc["W_star_as_defined"]
cond = dict(a_no_single_word_everywhere=(n_all == 0),
            b_no_word_in_Wstar_with_k3_largest=(bc['n_with_k3_largest'] == 0),
            c_no_word_in_Wstar_shorter_than_2pi=(bc['n_shorter_than_2pi'] == 0),
            d_at_least_99pct_of_sample_period_2pi=(frac_all >= 0.99),
            d_same_with_constants_excluded=(frac_nc >= 0.99))
bcn = res_bc["W_star_nondegenerate"]
cond_nd = dict(b_nondegenerate=(bcn['n_with_k3_largest'] == 0),
               c_nondegenerate=(bcn['n_shorter_than_2pi'] == 0))
print("  RULE CONDITIONS (as defined): %s" % cond)
print("  (b)(c) with W*' instead of W*: %s" % cond_nd)
print("  runtime %.0f s" % (time.time() - t0))

out = dict(
    description="Envelope word audit: which words carry the K3(delta) envelope, "
                "their Fourier spectra and periods.",
    engine="identical to lgi_period_and_purestate.py",
    n_delta=N_DELTA, Lmax=LMAX, n_words=N_WORDS, tolerance=TOL, seed=SEED,
    positive_control=dict(max_abs_dev_to_deposited_envelope=pc_dev, ok=pc_ok),
    cross_check=dict(n_words=len(xc_idx), max_abs_dev=xc_dev),
    degenerate_grid_points=dict(index=deg_idx, delta_over_pi=[float(deltas[j] / np.pi) for j in deg_idx],
                                note="all words have the same K3 within the tolerance"),
    words_tied_on_envelope_nondegenerate=dict(min=int(n_tied[~degenerate].min()),
                                              median=int(np.median(n_tied[~degenerate])),
                                              max=int(n_tied[~degenerate].max())),
    a_single_word_everywhere=dict(n_all_points=n_all, best_coverage_all=int(hits_all.max()),
                                  n_all_nondegenerate_points=n_all_nd,
                                  best_coverage_nondegenerate=int(hits_nd.max())),
    b_c=res_bc,
    W_star_nondegenerate_rows=dict(n_rows=len(rows), truncated=bool(Wnd.size > ROW_LIMIT), rows=rows),
    d_sample=dict(n=N_SAMPLE, n_constant=n_const, n_period_2pi=n_2pi,
                  fraction_period_2pi_of_all=frac_all,
                  fraction_period_2pi_of_nonconstant=frac_nc,
                  smallest_period_counts_m_nonconstant=s_per,
                  exceptions=exc),
    rule_conditions_as_defined=cond,
    rule_b_c_with_nondegenerate_Wstar=cond_nd,
)
with open(os.path.join(_HIER, 'lgi_envelope_word_audit_results.json'), 'w') as f:
    json.dump(out, f, indent=2)
print("\n-> lgi_envelope_word_audit_results.json written")
