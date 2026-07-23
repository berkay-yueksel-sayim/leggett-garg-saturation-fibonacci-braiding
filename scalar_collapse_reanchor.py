#!/usr/bin/env python3
"""
Paper 3 -- scalar-collapse re-anchor (Sec. IV.B, delta = 3 pi/5)
=================================================================
The printed sentence (Sec. IV.B) claims that at the sector-phase
singularity delta_star = 3 pi/5, sigma_1 becomes a scalar times the
identity to machine precision:

    ||sigma_1(delta_star) - e^{i theta} I|| = 3.1e-16

This script independently re-derives that number from the paper's own
generator construction (identical to lgi_fibonacci.py, this same
deposit, NOT edited in place) and checks it against a negative control
(a generic delta far from the singularity, where no such collapse is
expected).

Positive controls (must all pass before the re-anchored value is
trusted):
  - at delta=0 (the standard, undeformed Fibonacci theory), the two
    sigma_1 eigenvalues R_1 and R_tau are NOT equal (no accidental
    collapse) -- this is the negative canary for the singularity check.
  - the algebraic identity R_tau * e^{i delta_star} = R_1 holds to
    machine precision, matching the paper's own derivation
    (Sec. IV.B, Eq. after "so sigma_1 becomes a scalar").

Output: scalar_collapse_reanchor.json.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).parent

# ------------------------------------------------------------------
# Generators -- identical construction to lgi_fibonacci.py / k3_reanchor.py
# ------------------------------------------------------------------
R1 = np.exp(-4j * np.pi / 5)
Rt = np.exp(3j * np.pi / 5)


def sigma1(delta: float) -> np.ndarray:
    return np.diag([R1, Rt * np.exp(1j * delta)]).astype(complex)


# ------------------------------------------------------------------
# Negative control: at delta=0, no collapse (eigenvalues stay distinct)
# ------------------------------------------------------------------
s1_zero = sigma1(0.0)
distinct_err = float(abs(s1_zero[0, 0] - s1_zero[1, 1]))
negative_control = {
    "delta": 0.0,
    "abs_eigenvalue_difference": distinct_err,
    "pass_no_accidental_collapse": distinct_err > 0.1,
}

# ------------------------------------------------------------------
# Algebraic identity check: R_tau * e^{i delta_star} == R_1
# ------------------------------------------------------------------
delta_star = 3 * np.pi / 5
algebraic_err = float(abs(Rt * np.exp(1j * delta_star) - R1))
def c2list(z):
    z = complex(z)
    return [z.real, z.imag]


algebraic_identity = {
    "delta_star": delta_star,
    "R_tau_times_exp_i_delta_star_re_im": c2list(Rt * np.exp(1j * delta_star)),
    "R_1_re_im": c2list(R1),
    "abs_diff": algebraic_err,
    "pass": algebraic_err < 1e-10,
}

# ------------------------------------------------------------------
# Main re-anchor: ||sigma_1(delta_star) - e^{i theta} I||, theta = arg(R_1)
# ------------------------------------------------------------------
s1_star = sigma1(delta_star)
theta = float(np.angle(R1))
scalar = np.exp(1j * theta) * np.eye(2)
collapse_err = float(np.max(np.abs(s1_star - scalar)))

assert negative_control["pass_no_accidental_collapse"], f"NEGATIVE CONTROL FAILED: {negative_control}"
assert algebraic_identity["pass"], f"ALGEBRAIC IDENTITY FAILED: {algebraic_identity}"

paper_value = 3.1e-16
matches_paper = abs(round(collapse_err, 17) - 0) >= 0 and abs(collapse_err / paper_value - 1) < 0.5

out = {
    "meta": {
        "script": "scalar_collapse_reanchor.py",
        "purpose": "v1.2 re-anchor of the Sec. IV.B scalar-collapse claim at delta_star = 3 pi/5",
        "paper_reported_value": "3.1e-16",
    },
    "negative_control_delta_zero": negative_control,
    "algebraic_identity_R_tau_exp_i_delta_star_eq_R_1": algebraic_identity,
    "theta_arg_R1": theta,
    "collapse_deviation": {
        "delta_star": delta_star,
        "value": collapse_err,
        "rounded_2sig": float(f"{collapse_err:.1e}"),
        "matches_paper_3p1e-16": round(collapse_err, 18) > 0 and abs(collapse_err - paper_value) / paper_value < 0.1,
    },
}
(OUT / "scalar_collapse_reanchor.json").write_text(json.dumps(out, indent=2))

print("=== Negative control (delta=0, no collapse expected) ===")
print(f"  |R_1 - R_tau| = {distinct_err:.4f}  pass={negative_control['pass_no_accidental_collapse']}")
print("\n=== Algebraic identity R_tau * e^{i delta_star} == R_1 ===")
print(f"  abs diff = {algebraic_err:.2e}  pass={algebraic_identity['pass']}")
print("\n=== Scalar-collapse re-anchor at delta_star = 3 pi/5 ===")
print(f"  ||sigma_1(delta_star) - e^(i theta) I|| = {collapse_err:.10e}")
print(f"  rounded (2 sig figs) = {collapse_err:.1e}")
print(f"  paper reports 3.1e-16 -- match: {out['collapse_deviation']['matches_paper_3p1e-16']}")
print(f"\nWRITE {OUT / 'scalar_collapse_reanchor.json'}")
