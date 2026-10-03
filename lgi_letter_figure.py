"""
Build the two-panel Letter figure from the three result JSONs:
  lgi_results.json                          (Fibonacci main)
  lgi_ising_results.json                    (Ising)
  lgi_period_and_purestate_results.json     (Fibonacci period scan
                                             + pure-state robustness)

Output: lgi_letter_figure.png

  Panel A  best K3 versus braid length L for Fibonacci and Ising.
           Horizontal markers at K3 = 1 (macrorealism) and 3/2 (Lueders).
  Panel B  K3(delta) deformation scan for Fibonacci, with the
           scalar-singularity point delta = 3 pi / 5 marked. The Fourier-
           period annotation of earlier versions is withdrawn.
"""
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))


def _find(name):
    """Look for a result JSON next to this script first, then one level up."""
    for cand in (os.path.join(HERE, name),
                 os.path.join(os.path.abspath(os.path.join(HERE, os.pardir)),
                              name)):
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(name)


with open(_find('lgi_results.json')) as f:
    FIB = json.load(f)
with open(_find('lgi_ising_results.json')) as f:
    ISI = json.load(f)
with open(_find('lgi_period_and_purestate_results.json')) as f:
    EXT = json.load(f)


# ---------- gather data --------------------------------------------------------
fib_L = sorted(int(k) for k in FIB['exhaustive_delta0'])
fib_K3 = [FIB['exhaustive_delta0'][str(L)]['K3'] for L in fib_L]
fib_random_L = sorted(int(k) for k in FIB['random_convergence'])
fib_random_K3 = [FIB['random_convergence'][str(L)] for L in fib_random_L]

isi_L = sorted(int(k) for k in ISI['exhaustive_delta0'])
isi_K3 = [ISI['exhaustive_delta0'][str(L)]['K3'] for L in isi_L]
isi_random_L = sorted(int(k) for k in ISI['random_convergence'])
isi_random_K3 = [ISI['random_convergence'][str(L)] for L in isi_random_L]

delta_over_pi = np.array(EXT['period_scan']['delta_over_pi'])
K3_sweep = np.array(EXT['period_scan']['K3'])


# ---------- plot ---------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 9,
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
})

fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), constrained_layout=True)


# --- Panel A:  K3(L) for Fibonacci vs. Ising ---
axA = axes[0]
axA.axhline(1.0, color='gray', linestyle=':', linewidth=0.7, zorder=1)
axA.axhline(1.5, color='gray', linestyle=':', linewidth=0.7, zorder=1)
axA.text(40.0, 1.01, 'macrorealism  $K_3 \\leq 1$', ha='right', va='bottom',
         fontsize=8, color='gray')
axA.text(40.0, 1.51, u'Lüders  $K_3 \\leq 3/2$', ha='right', va='bottom',
         fontsize=8, color='gray')

axA.plot(fib_L, fib_K3, 'o-', markersize=4, linewidth=1.0,
         color='black', label='Fibonacci (exhaustive)', zorder=4)

# Mark the non-monotonicity dips at L = 5, 7, 10
DIP_LS = [5, 7, 10]
for L_dip in DIP_LS:
    K_dip = fib_K3[L_dip - 1]
    axA.plot(L_dip, K_dip, 'o', markersize=8, markerfacecolor='none',
             markeredgecolor='#3366cc', markeredgewidth=0.9, zorder=5)
axA.annotate('dips at $L=5,7,10$',
             xy=(5, fib_K3[4]), xytext=(13, 1.13),
             fontsize=7.5, color='#3366cc', ha='left', va='center',
             arrowprops=dict(arrowstyle='-', color='#3366cc',
                             linewidth=0.5, shrinkA=0, shrinkB=2))
axA.plot(fib_random_L, fib_random_K3, 's', markersize=3.5,
         markerfacecolor='white', markeredgecolor='black', markeredgewidth=0.7,
         label='Fibonacci (random, $4 \\cdot 10^5$)', zorder=3)

axA.plot(isi_L, isi_K3, '^-', markersize=4, linewidth=1.0,
         color='#cc4040', label='Ising (exhaustive)', zorder=4)
axA.plot(isi_random_L, isi_random_K3, 'v', markersize=3.5,
         markerfacecolor='white', markeredgecolor='#cc4040', markeredgewidth=0.7,
         label='Ising (random, $4 \\cdot 10^5$)', zorder=3)

axA.set_xlabel('braid word length $L$')
axA.set_ylabel(r'best $K_3$')
axA.set_xlim(0, 41)
axA.set_ylim(0.95, 1.56)
axA.xaxis.set_major_locator(MultipleLocator(5))
axA.xaxis.set_minor_locator(MultipleLocator(1))
axA.yaxis.set_major_locator(MultipleLocator(0.1))
axA.yaxis.set_minor_locator(MultipleLocator(0.05))
axA.legend(loc='center right', fontsize=7.5, framealpha=0.95,
           edgecolor='black', fancybox=False).get_frame().set_linewidth(0.5)
axA.set_title('(a) Universality split', fontsize=9, loc='left')


# --- Panel B:  K3(delta) deformation scan, Fibonacci ---
axB = axes[1]
axB.axhline(1.0, color='gray', linestyle=':', linewidth=0.7, zorder=1)
axB.axhline(1.5, color='gray', linestyle=':', linewidth=0.7, zorder=1)

# the scalar-singularity point delta = 3 pi / 5  =  0.6 pi
axB.axvline(0.6, color='#3366cc', linestyle='--', linewidth=0.7, zorder=2)
# label sits high on the dashed line, clear of the curve
axB.text(0.62, 1.53, r'scalar gen.  $\delta = 3\pi/5$',
         fontsize=7.5, color='#3366cc', ha='left', va='top')

axB.plot(delta_over_pi, K3_sweep, '-', linewidth=1.0, color='black', zorder=4)

axB.set_xlabel(r'deformation parameter  $\delta / \pi$')
axB.set_ylabel(r'best $K_3$  (exhaustive $L \leq 9$)')
axB.set_xlim(0, 2)
axB.set_ylim(0.95, 1.56)
axB.xaxis.set_major_locator(MultipleLocator(0.5))
axB.xaxis.set_minor_locator(MultipleLocator(0.1))
axB.yaxis.set_major_locator(MultipleLocator(0.1))
axB.yaxis.set_minor_locator(MultipleLocator(0.05))

axB.set_title('(b) Deformation dependence (Fibonacci)', fontsize=9, loc='left')


outpath = os.path.join(HERE, 'lgi_letter_figure.png')
fig.savefig(outpath, dpi=300, bbox_inches='tight')
print(f"-> {outpath}")
print(f"   panel A : Fibonacci vs Ising K3(L)")
print("   panel B : Fibonacci K3(delta), scalar singularity marked")
