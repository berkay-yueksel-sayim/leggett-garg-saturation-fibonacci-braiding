# LGI Letter v1.2 — Zenodo Build

**Title:** Leggett–Garg K₃ Values Above 1 in Fibonacci-Anyon Braiding: 99.998 % of the Lüders Bound, and Exactly 1 for Ising Braiding
**Author:** Berkay Yüksel Sayim
**ORCID:** [0009-0004-4993-7352](https://orcid.org/0009-0004-4993-7352)
**Affiliation:** Independent Research, Germany
**Version:** 1.2
**Date:** 2026-07-23
**Resource type:** Preprint
**License:** paper, figures, and data — [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), see `LICENSE`; source code (`*.py`) — [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), see `LICENSE-CODE`

## v1.2 changes (this release)

- Two orphan bibliography entries resolved: `EmaryLambertNori2014` is now cited
  (moving-bound formula $K_n^{\max}=n\cos(\pi/n)$, Open Question O1); `Fine1982`
  is removed (no body citation existed for it in v1.1).
- Three new citations added: Fritz 2010 (closed-form temporal-CHSH correlator),
  Emary 2013 (decoherence/noise-threshold framework), Kofler--Brukner 2008
  (conditions for quantum violation of macroscopic realism).
- The bare "saturates ... already at $L=9$" table caption now carries the
  $99.998\%$/never-exact qualifier used elsewhere in the text.
- A companion-work citation to the $\mathrm{SU}(2)_k$ Leggett--Garg study
  (Concept-DOI 10.5281/zenodo.20531124) is added at Open Question O3, plus a
  one-sentence limitation noting the result is for projective L\"uders
  measurement (weak/non-projective protocols untested).
- Title hyphen corrected to an en dash ("Leggett--Garg") for series
  consistency with the companion papers; PDF metadata title is left as ASCII
  hyphen.
- The v1.1 in-document correction-note block is removed (its content is
  preserved here and in the Zenodo version history); `%% v1.1` scaffolding
  markers are stripped.
- No numerical result, figure, or existing body claim is changed.

## Additional corrections (builder-fidelity review, wave R-B, same v1.2 release)

A pre-release builder-fidelity review wave (each simulation strand checking
the faithful rendition of its own contribution) produced five further
corrections; none changes a numerical result, table, or figure.

- **Attribution of the Ising no-violation result (three sites).** The
  impossibility of a spatial-CHSH violation by Ising braiding alone is now
  attributed to its primary source, Howard and Vala (Phys. Rev. A 85,
  022304, 2012), consistent with the companion papers; Clarke, Sau, and
  Das Sarma (Phys. Rev. X 6, 021005, 2016) are repositioned as supplying
  the enabling non-Clifford phase gate for Majorana wires (a constructive
  result, not the prohibition).
- **Fibonacci CHSH saturation attribution.** The saturation statement is
  now carried by braid-representation density (Nayak et al.), with Brennen
  et al. cited for their explicit --- sub-Tsirelson --- CHSH-violating
  settings.
- **Open Question O2 restated in two stages:** the proximate cause of the
  $L=5,7,10$ dips is the reachability geometry (established); whether the
  divisibility structure of the central element $(\sigma_1\sigma_2)^3$ is
  the deeper algebraic origin remains open.
- **Wording precision:** two qubit involutions "have an anticommutator
  proportional to the identity" (previously "anticommute up to a multiple
  of the identity"); the formula $\{Q_1,Q_2\}\propto I$ was always correct.
- **Cross-reference added:** the scalar-collapse threshold
  $\delta_\star=3\pi/5$ is identified with the intrinsic $\sigma_1$
  rotation angle $\delta^{\ast}(k{=}3)=\pi k/(k+2)$ of the companion
  $\mathrm{SU}(2)_k$ analysis; both are the relative phase
  $\arg(R_1/R_\tau)$ of the two braid eigenvalues.
- Release date set (2026-07-17).

## Additional corrections (internal review pass, same v1.2 release)

A subsequent fresh-context review pass found two further items inherited
from v1.0–v1.2 and not caught by the changes above; both are corrected here.

- **Argmax-word transcription (Sec.~III.A, "L\"uders saturation").** The
  printed word said to attain $K_3^{\max}=1.499762$ at $L=11$ was a
  mistranscribed 14-letter string (it decodes to the $L=10$ value
  $1.497834$, not the $L=11$ value). The true $L=11$ argmax word is
  re-derived directly from the deposited search code
  (`lgi_fibonacci.py`) via two algorithmically independent enumeration
  routes (incremental engine-faithful enumeration, cross-checked
  bit-exactly against the archived `lgi_results.json` for $L=1\ldots10$;
  and an independent meet-in-the-middle $L=5+6$ decomposition), which
  agree to $2\times10^{-16}$ and confirm $352$ words tie the maximum.
  The corrected 11-letter word
  ($\sigma_1\sigma_2^2\sigma_1^{-3}\sigma_2^2\sigma_1^{-1}\sigma_2\sigma_1$)
  replaces the mistranscribed one; **$K_3^{\max}=1.499762$ itself is
  unchanged** and was already correct in Table~I. Script and archived
  output: `k3_reanchor.py` / `k3_reanchor.json`.
- **Retracted length-divisibility claim (Sec.~III.A).** A sentence
  claiming braid lengths divisible by $6$ are "structurally favored"
  and that primes such as $11$ are "structurally further away" directly
  contradicted the paper's own Table~I, where $L=11$ (not divisible by
  $6$) ties the exhaustive maximum with $L=9$ while $L=6$ reaches only
  $K_3=1.47744$. The unsupported claim is retracted and replaced with
  the data-grounded statement.
- **State-independence mechanism (Sec.~III.C, "Initial-state
  independence").** The sentence attributing the pure-state /
  mixed-state agreement to $ZUZU^{\dagger}$ and
  $ZU^2Z(U^2)^{\dagger}$ being "proportional to the identity for the
  optimal $U$" was mathematically false (this holds, at most, for the
  Hermitian part, and for *every* unitary $U$, not a property singled
  out by the optimal word). It is replaced by the correct, general
  argument: the two-time L\"uders correlator for dichotomic observables
  is $C(Q_1,Q_2)=\mathrm{Tr}[\rho\,\tfrac12\{Q_1,Q_2\}]$ in any
  dimension, and on the 2D fusion space used here $Z$ and
  $UZU^{\dagger}$ are traceless $\pm1$-eigenvalue qubit operators for
  *every* $U$, so $\{Q_1,Q_2\}\propto I$ — a generic $d=2$ identity, not
  a finding. The subsection and the corresponding abstract sentence
  (finding (v)) are reframed as a consistency check rather than an
  independent result; no numerical value changes.
- Minor consistency/hedging fixes: "artefact"→"artifact" (AmE, 3
  occurrences); the Ising-sampling range is stated consistently as
  "exhaustive $L\leq11$, random even $L\in\{12,14,\ldots,40\}$" in the
  abstract, figure caption, and body (previously "every $L\leq40$" in
  the abstract overstated the sampled range); the Ising $L=11$
  $(C(B),C(B^2))=(0,-1)$ statement is scoped to "the $L=11$ argmax
  word" rather than implied as unique; two unhedged priority claims
  gain "to our knowledge".

## Additional corrections (pre-publication verification pass, same v1.2 release, 2026-07-23)

A final pre-publication verification pass re-anchored every claim of the
sector-dependence subsection against the deposited data; three items are
corrected.

- **Sector-dependence paragraph (Sec. III.D) re-anchored to the deposited
  scan.** The previous text stated that $K_3 = 1$ is reached also at
  nonsingular $\delta$ values and that random sampling at
  $\delta/\pi = 1.556$ with $L=24$ recovers $K_3 = 1.4999996$; neither
  statement reproduces from the deposited 240-point period scan (the only
  $K_3 = 1.0$ point is the singularity $\delta^{\ast} = 3\pi/5$ itself,
  and at $\delta/\pi = 1.5583$ the $L\leq9$ envelope already reaches
  $1.49995$). The paragraph now states the deposited facts: the weakest
  nonsingular sectors are the two grid neighbors of the singularity
  ($K_3 \approx 1.051$), and a new seeded validator
  (`lgi_sector_budget_validator.py`, seed 20260723) shows the monotone
  budget recovery there ($K_3 = 1.075$ at $L=12$, $1.148$ at $L=24$),
  well short of the L\"uders bound within the deposited budgets.
- **Abstract precision:** the initial-state-independence agreement is
  stated as $\sim 10^{-15}$ (machine precision), matching the deposited
  spread $1.1\times10^{-15}$ (previously "$10^{-16}$").
- **This README:** the residual claims that all sanity checks pass "in
  every reproduction script" (Build section) and that fixed-word traces
  carry a $2\pi/6$ period (three places) are corrected to match the
  deposited scripts and the paper body ($2\pi$-generic fixed-word
  periodicity); the Budroni-Emary reference gains its arXiv ID.

## v1.1 Changes (corrections to v1.0)

This release corrects four items in v1.0 (Zenodo Concept-DOI
[10.5281/zenodo.20372744](https://doi.org/10.5281/zenodo.20372744)).
**No numerical results, no figures, and no other body claims are affected.**

1. **Fourier-harmonic correction (envelope vs fixed-word).**
   v1.0 reported the dominant Fourier harmonic of $K_3(\delta)$ as $k=6$
   (period $2\pi/6$), based on a sample at $L_{\max} = 9$. Independent
   reproducer runs at $L_{\max} = 7$ and a pre-v1.1 sanity re-run at
   $L_{\max} = 8$ both find $k=3$ as the dominant harmonic of the
   *envelope* $K_{3,\max}(\delta) = \max_{|w| \leq L} K_3(\delta;\,w)$,
   with stable $k=3/k=6$ amplitude ratio 1.061 → 1.063. The algebraic
   $(\sigma_1\sigma_2)^3$ identity does not impose a $2\pi/6$ period on
   individual fixed words: fixed-word traces $K_3(\delta;\,w)$ are
   generically only $2\pi$-periodic, while the envelope shows
   optimal-word reshuffling on a finer scale. The envelope period
   ($2\pi/3$) and the generic fixed-word periodicity ($2\pi$) are now
   stated separately in the abstract, overview, figure caption,
   and §III.D body.

2. **Hou et al. → Clarke-Sau-Das Sarma (PRX 6, 021005).** The
   bibliography entry previously labeled "Hou et al." at PRX 6, 021005
   (2016) is corrected to its true attribution: D. J. Clarke, J. D. Sau,
   and S. Das Sarma, *A practical phase gate for producing Bell
   violations in Majorana wires*. Body text "Hou–Shtengel split" is
   updated to "Clarke–Sau–Das-Sarma split" throughout. Body content
   (Clifford-only Ising braiding, missing non-Clifford direction) is
   unchanged — only the attribution is fixed.

3. **Sorella 2023 bibliography entry corrected.** Title corrected from
   the v1.0 placeholder *"Representation dependence of the Tsirelson
   bound"* to the verified title *"On the representations of Bell's
   operators in Quantum Mechanics"*; publication venue updated to
   Foundations of Physics 53, 59 (2023), DOI
   [10.1007/s10701-023-00699-6](https://doi.org/10.1007/s10701-023-00699-6).

4. **Minev 2025 bibliography entry corrected.** Author list expanded to
   the full 8-author form (with K. Najafi as second author, not "S.
   Najafi" as in v1.0); title updated to *"Realizing string-net
   condensation: Fibonacci anyon braiding for universal gates and
   sampling chromatic polynomials"*; article number 6225 added.

## Abstract

We numerically test the three-time Leggett–Garg inequality
$K_3 \leq 1$ for the standard B$_3$ Fibonacci-anyon braiding
representation on the two-dimensional fusion space of three $\tau$
anyons. Exhaustive enumeration over all $4^L$ braid words up to
$L = 11$ and random sampling to $L = 40$ show that $K_3$ saturates
the Lüders bound $3/2$ to $99.998\%$, with the first violation
already at $L = 3$. Two structural signatures accompany the saturation.
First, replacing Fibonacci by Ising generators on the
same 2D fusion space gives $K_3 = 1$ exactly for every $L \leq 11$ in the
exhaustive search and every random word tested at even $L \in \{12, 14,
\ldots, 40\}$, a sharp split mirroring the Howard–Vala no-Bell-violation
result for Ising braiding in the spatial CHSH setting. Second, an
algebraic phase-deformation parameter $\delta$ tunes a singular point
$\delta = 3\pi/5$ at
which the generator $\sigma_1$ collapses to a scalar to machine
precision, so that every braid word reduces to a global phase. The
Fourier-period analysis of the $\delta$ dependence reported in earlier
versions of this letter is withdrawn. $K_3$ at
the optimal $L=11$ word is initial-state independent.

---

## What's in this archive

| File | Role |
|---|---|
| `main_v1.3.tex` | LaTeX source (RevTeX 4-2, PRX style) |
| `main_v1.3.pdf` | Compiled preprint |
| `lgi_letter_figure.png` | Two-panel preprint figure |
| `lgi_letter_figure.py` | Preprint-figure builder |
| `lgi_fibonacci.py` | Main Fibonacci LGI computation |
| `lgi_results.json` | Fibonacci raw results |
| `lgi_ising.py` | Ising LGI computation (universality split) |
| `lgi_ising_results.json` | Ising raw results |
| `lgi_period_and_purestate.py` | $K_3(\delta)$ period scan + pure-state robustness |
| `lgi_period_and_purestate_results.json` | period and pure-state raw results |
| `k3_reanchor.py` | Independent $L=11$ argmax-word re-derivation (2 routes) |
| `k3_reanchor.json` | `k3_reanchor.py` output (both routes, tie count, archive cross-check) |
| `scalar_collapse_reanchor.py` | Re-anchor of the Sec. III.C scalar-collapse value at $\delta_\star = 3\pi/5$ |
| `scalar_collapse_reanchor.json` | `scalar_collapse_reanchor.py` output (negative control, algebraic identity, collapse deviation) |
| `lgi_sector_budget_validator.py` | Seeded random-search validation of the braid-budget statement at the weakest nonsingular grid point ($\delta/\pi = 0.6083$) |
| `lgi_sector_budget_validator.json` | `lgi_sector_budget_validator.py` output (best $K_3$ at $L=12/24$, positive control at $\delta=0$) |
| `lgi_envelope_word_audit.py` | Exhaustive audit of the words that carry the envelope $K_{3,\max}(\delta)$: coverage, Fourier spectra and periods of their own curves |
| `lgi_envelope_word_audit_results.json` | `lgi_envelope_word_audit.py` output (positive control against the deposited envelope, coverage, harmonics, periods, random sample) |
| `lgi_braid_relation_check.py` | Braid relation and mirror symmetry of the deformed generators (Secs. II B, III D) |
| `lgi_braid_relation_check_results.json` | `lgi_braid_relation_check.py` output (grid, fine scan, scalar and conjugate points, envelope and word-level mirror symmetry) |
| `LICENSE` | CC BY 4.0 — paper, figures, data |
| `LICENSE-CODE` | Apache License 2.0 — source code (`*.py`) |
| `README.md` | This file |

## Reproduction

Deterministic. NumPy and Matplotlib only.

```bash
python lgi_fibonacci.py              # Fibonacci main, seed 20260524 (~2 min)
python lgi_ising.py                  # Ising, seed 20260525 (~2 min)
python lgi_period_and_purestate.py   # period + pure-state, seed 20260525 (~5.5 min)
python lgi_letter_figure.py          # preprint figure (~4 s)
python lgi_sector_budget_validator.py  # sector budget validator, seed 20260723 (~9 s)
python k3_reanchor.py                # L=11 argmax-word re-derivation (~23 s)
python scalar_collapse_reanchor.py   # scalar-collapse re-anchor (~1 s)
python lgi_envelope_word_audit.py    # envelope word audit, seed 20260928 (~5.7 min)
python lgi_braid_relation_check.py   # braid relation + mirror symmetry, no random numbers (~1 s)
```

`lgi_fibonacci.py` and `lgi_ising.py` open with sanity checks (unitarity,
Yang–Baxter, $(\sigma_1\sigma_2)^3$ scalar) that must pass at machine
precision; if they do not, the generator conventions have been altered.
`lgi_period_and_purestate.py` carries a separate structural positive
control for its Fourier fit (`k_ctrl == 6`) but not the three checks
above; `lgi_letter_figure.py` only loads the three `*_results.json`
files and plots them, with no numerical checks of its own.
The times in the comments are wall-clock times from one run each on a single
desktop machine (Python 3.12.10, NumPy 2.4.3); each run reproduced its
deposited output file byte for byte. The two stability runs in
`lgi_period_and_purestate.py` ($L_{\max} = 7$ and $8$, $n_\delta = 120$) take
the maximum over all words with $|w| \leq L_{\max}$, i.e. 21,844 and 87,380
words.

## Key numerical claims (independently verifiable from the JSON files)

| Claim | Value | File / field |
|---|---|---|
| Fibonacci $K_3^{\max}$ (exhaustive, $L=11$) | 1.499762 | `lgi_results.json` → `best` |
| Fibonacci $K_3^{\max}$ (random, $L \geq 22$) | 1.499964 | `lgi_results.json` → `random_convergence` |
| First Fibonacci LGI violation | $L=3$, $K_3=1.292$ | `lgi_results.json` → `exhaustive_delta0.3` |
| Ising $K_3^{\max}$ at $\delta=0$ over $L=1\ldots 40$ | 1.000000 (exact) | `lgi_ising_results.json` |
| Ising $K_3^{\max}$ over deformed $\delta$ | 1.500 at $\delta/\pi \approx 1.167$ | `lgi_ising_results.json` → `sector_sweep` |
| Envelope $K_{3,\max}(\delta)$ not $2\pi/3$-periodic | $K_{3,\max} = 1$ only at $\delta/\pi = 0.6$; $1.2047$ at $\delta/\pi = 0.6 + 2/3$ | `lgi_period_and_purestate_results.json` → `period_scan.K3[72]`, `period_scan.K3[152]` |
| No single word attains the envelope at every grid point | best word: 10 of the 239 nonsingular points | `lgi_envelope_word_audit_results.json` → `a_single_word_everywhere.best_coverage_nondegenerate` |
| Envelope words with $k=3$ as the strongest harmonic of their own curve | 4,420 of 19,362 | `lgi_envelope_word_audit_results.json` → `b_c.W_star_nondegenerate` |
| Fixed-word curves $K_3(\delta;\,w)$ with a period shorter than $2\pi$ (random sample of 1000 words) | 42 of 924 nonconstant (76 constant, 882 $2\pi$-periodic) | `lgi_envelope_word_audit_results.json` → `d_sample` |
| Braid relation $\sigma_1\sigma_2\sigma_1 = \sigma_2\sigma_1\sigma_2$ | exact only at $\delta/\pi = 0, 0.6, 1.2$; elsewhere $\geq 0.0291$ (Frobenius), also up to a global phase | `lgi_braid_relation_check_results.json` → `grid`, `fine_scan` |
| Mirror symmetry $K_{3,\max}(\delta) = K_{3,\max}(6\pi/5 - \delta)$ | $6.7 \times 10^{-15}$ about both axes ($3\pi/5$, $8\pi/5$); wrong axis $0.4996$ | `lgi_braid_relation_check_results.json` → `envelope_mirror` |
| Pure-state $K_3$ at $L=11$ optimum (9 states) | 1.499762 ± 0 | `lgi_period_and_purestate_results.json` → `pure_state_check` |
| Scalar-collapse deviation at $\delta_\star = 3\pi/5$ | $3.1 \times 10^{-16}$ | `scalar_collapse_reanchor.json` → `collapse_deviation` |
| Sanity checks (unitarity, Yang–Baxter, $(\sigma_1\sigma_2)^3$ scalar) | $\leq 6 \times 10^{-16}$ | `lgi_results.json` and `lgi_ising_results.json` → `sanity`; `lgi_period_and_purestate_results.json` carries no `sanity` field (see Reproduction section) |


---

## Build

- LaTeX engine: MiKTeX pdfTeX-1.40.29 (3 × pdflatex passes, converged).
- Python: 3.12, NumPy and Matplotlib only.
- Sanity checks (unitarity, Yang–Baxter, $(\sigma_1\sigma_2)^3$ scalar)
  pass at machine precision in `lgi_fibonacci.py` and `lgi_ising.py`;
  the other reproduction scripts carry no such checks
  (see Reproduction section).

*Supplementary files tidied 2026-06-22; no result, value, table, or figure changed.*
