# Why the Sequential EKF Needs σ = 5° Despite 1″ Observations

## The Problem

The USNO transit circle observations (1963–1971) have ~1″ precision, yet the
sequential Extended Kalman Filter diverges with σ = 1″ and only converges
with σ = 5° (18,000× larger than the actual measurement noise). This is not
a bug — it's a fundamental limitation of sequential filtering for this class
of problem.

## What Happens with σ = 1″

The EKF processes observations **one at a time**. For each observation it:

1. Predicts RA/Dec from the current state (15 orbital parameters)
2. Computes the innovation (observed − predicted)
3. Computes the Kalman gain K = P·Hᵀ·(H·P·Hᵀ + R)⁻¹
4. Updates the state: x ← x + K·innovation

With σ = 1″, the measurement noise R = (1/3600°)² is tiny. The innovation
covariance S = H·P·Hᵀ + R is dominated by the prior term H·P·Hᵀ, making
the Kalman gain K ≈ 1. This means the filter **trusts each single
observation completely**.

Even starting 0.002° from truth, the first observation has a ~40″ innovation.
The filter applies a 40″ correction distributed across 15 parameters using a
linearized Jacobian. But the mapping from orbital elements to RA/Dec is
highly nonlinear — the Jacobian is only accurate for infinitesimal
perturbations, not for 40″ corrections in a 15-dimensional space. The
correction overshoots, the state gets worse, the next observation has an even
larger innovation, and the filter spirals into divergence.

**Tested:** Even with perturbations as small as 0.002° in angles and 0.01%
in scales, σ = 1″ diverges by observation #200 (eccentricity error reaches
1500%). The problem is not the size of the initial perturbation — it's that
any finite perturbation produces innovations much larger than 1″, and the
per-observation update mechanism cannot handle this.

## Why σ = 5° Works

With σ = 5°, R = (5°)² dominates S, so the Kalman gain K ≈ 0.01. Each
observation nudges the state by only ~1% of the innovation. After processing
1260 observations, the cumulative effect converges to a good estimate.

The tradeoff: the posterior covariance reflects 5° noise, not 1″. The
reported uncertainties (e.g., σ_aE = 0.17%) are ~18,000× larger than what
the data actually supports. The point estimates are correct; the error bars
are not.

## Current Results (σ = 5°, 5 iterations)

From 745 Sun + 515 Moon transit observations spanning 8 years:

| Quantity       | Estimated          | True               | Error   |
|----------------|--------------------|--------------------|---------|
| GM_sun         | 1.3287 × 10¹¹     | 1.3271 × 10¹¹     | 0.12%   |
| Earth period   | 364.43 days        | 364.81 days        | 0.11%   |
| Moon period    | 27.521 days        | 27.520 days        | 0.004%  |
| a_Earth        | 149,432,000 km     | 149,476,000 km     | 0.03%   |
| a_Moon         | 386,537 km         | 386,610 km         | 0.02%   |

Known degeneracies (expected with geocentric angular data):
- ω_E–M₀_E correlation: −0.90 (mean longitude constrained, not individual angles)
- μ_SE–a_E correlation: 0.997 (mean motion constrained, not individual μ and a)
- GM_moon poorly determined (81% error) — see below

## Geocentric Data Limitations

The USNO published **geocentric** RA/Dec — corrected from the observatory's
topocentric position using a pre-computed lunar ephemeris. This correction
removed the diurnal parallax signal (~57′ for the Moon) that would directly
measure distances. We cannot recover absolute distances from this data;
the μ–a degeneracy is only weakly broken by the **parallactic inequality**
(the Sun's perturbation of the Moon's orbit, amplitude ~125″, depends on
the ratio a_Moon/a_Earth).

GM_moon is essentially unobservable from geocentric data because it enters
only through the Earth-Moon barycenter offset (~4700 km), which produces
a ~0.02° signal in the Sun's apparent position — well below the σ = 5° noise
floor.

## The Fix: Batch Least Squares

The correct algorithm for orbit determination with precise observations is
**batch least squares** (Gauss-Newton or Levenberg-Marquardt):

1. Compute residuals for ALL observations at the current state
2. Compute the Jacobian H_i for each observation (same state for all)
3. Form normal equations: (Σ Hᵢᵀ R⁻¹ Hᵢ) Δx = Σ Hᵢᵀ R⁻¹ rᵢ
4. Solve for Δx, update x ← x + Δx
5. Repeat until convergence

This works with σ = 1″ because:
- The correction Δx averages over all 1260 observations simultaneously
- A single large residual cannot destroy the state (it's one of 1260 terms)
- The linearization only needs to be accurate for the single correction Δx,
  not for 1260 sequential updates
- Levenberg-Marquardt damping provides additional robustness

The computational cost per iteration is the same (1260 × 16 Jacobian
evaluations), but batch LS converges in 5–10 iterations with honest
uncertainties reflecting the actual 1″ measurement precision.

## Summary

| Approach            | σ_obs | Converges? | Honest σ? | Handles large perturbations? |
|---------------------|-------|------------|-----------|------------------------------|
| Sequential EKF      | 1″    | No         | —         | No                           |
| Sequential EKF      | 5°    | Yes        | No        | Yes (up to ~1° angles)       |
| Batch LS (planned)  | 1″    | Yes        | Yes       | Yes (with LM damping)        |
