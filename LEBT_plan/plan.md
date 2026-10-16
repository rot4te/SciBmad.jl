# Plan: Space-charge (low-energy) beam support in SciBmad

Goal: make the SciBmad stack (SciBmad.jl, BeamTracking.jl, Beamlines.jl, SpaceCharge.jl)
able to track low-energy, space-charge-dominated beams.

Status: survey only, no code changed. Items marked **(verify)** are inferences that were not
confirmed by reading the full code.

Paths below are relative to `~/.julia/dev/`.

---

## Current state

- SpaceCharge.jl is a standalone solver: CIC deposition, an FFT free-space Poisson solve with
  integrated Green functions, and field interpolation. Nothing in SciBmad, BeamTracking,
  Beamlines or SciBmadSimUtils references it.
- `SciBmad.jl/README.md` lists "Space charge" as a planned feature.
- The only collective effect, IBS, is applied **once per element**. It uses beam moments
  computed on the host before the kernel launches
  (`BeamTracking.jl/ext/BeamTrackingBeamlinesExt/general_bl.jl`, `ibs_kick`).

---

## Step 1: Add a collective-effect synchronization point (BeamTracking)

**Problem.** `launch!` (`BeamTracking.jl/src/kernel.jl`) runs one kernel per element. Each
particle goes through all integrator steps of that element on its own
(`order_two_integrator!` etc. in `src/kernels/integrators.jl`). Callbacks
(`src/callbacks.jl`) are per-particle closures that run inside the kernel, so they can't see
the other particles. Space charge needs deposit → solve → kick across the whole bunch
several times per element.

**Tasks.**
- [ ] Choose an approach:
  - **A. Split launches:** split each element's `n_steps` into separate kernel launches, with
    a host-side collective step between them (split-operator: half map, SC kick, half map).
    Keep the current fused single-launch path when space charge is off.
  - **B. Element slicing in Beamlines:** add slicing to Beamlines (`src/superimpose.jl` is
    empty; there is no slicing today). Fringes, bend geometry, alignment and apertures must
    be correct on thin slices.
- [ ] Give drifts multiple steps too (today they are one `exact_drift!`).
- [ ] Define a general "collective kick" hook so IBS can move to the same per-step mechanism.

## Step 2: Coordinate and frame conversion (BeamTracking ↔ SpaceCharge)

**Problem.** BeamTracking uses s as the independent variable: `z` is position relative to the
reference particle at fixed s. SpaceCharge.jl expects positions at a fixed time. At low
energy, velocity spread makes these noticeably different.

**Tasks.**
- [ ] Before depositing, drift particles to a common time. `compute_time` in
  `src/utils/z_to_time.jl` gives the relation. After kicking, map back to fixed-s
  coordinates.
- [ ] Deposit only `STATE_ALIVE` particles. `Mesh3D` takes its bounds from the `extrema` of
  all particle coordinates, so lost particles and halo outliers would stretch the grid.
- [ ] Decide how to handle outliers: clipping, a separate treatment, or a fixed-size mesh.
- [ ] Handle curvature: the solver geometry is straight and ignores `g` inside bends.

## Step 3: Energy spread and the single-γ assumption (SpaceCharge)

**Problem.** `Mesh3D` has one `gamma` for the whole bunch. Its Green function stretches z by
γ and scales the transverse fields by γ (a rest-frame electrostatic solve). This fails with
large relative energy spread, as near a cathode or in chirped bunches.

**Tasks.**
- [ ] Add energy binning (one solve per energy slice, then sum the fields), or document the
  limits where a single γ is acceptable.
- [ ] Check the Green-function caching (keyed on `delta`, `gamma`, `offset`) still pays off
  when the mesh is rebuilt every step.

## Step 4: Fields → momentum kicks (BeamTracking)

**Tasks.**
- [ ] Write a space-charge kick kernel. SpaceCharge.jl returns **E only**, so the kernel must
  also apply the magnetic part (B = β×E/c), giving the 1/γ² cancellation, or use the
  rest-frame/boost form directly.
- [ ] Convert to s-based kicks using each particle's dt/ds, following the pattern in
  `ibs_damping_and_diffusion!` (`src/kernels/ibs_kick.jl`).
- [ ] Define the charge normalization. `Bunch` has no bunch-charge or
  charge-per-macroparticle field, and IBS treats `weight` as a particle count. Settle one
  convention for `weight` and add bunch charge if needed.

## Step 5: Boundary conditions (SpaceCharge)

**Problem.** Only free space and an image-charge cathode at z = 0
(`solve!(…; at_cathode=true)`) are supported.

**Tasks.**
- [ ] Add conducting-pipe boundaries (round/rectangular image fields) for low-energy rings
  and transport lines.

## Step 6: Guns and cathodes (only if in scope)

**Problem.** At β ≈ 0, s-based tracking breaks down. Particles are marked lost when the
longitudinal momentum squared drops to zero or below (`exact_drift!`), or when they move
backward (the dt/ds check in IBS).

**Tasks.**
- [ ] Add a time-based tracking mode (t as the independent variable) for the gun region.
- [ ] Emit particles from the cathode using `STATE_PREBORN` → `STATE_ALIVE`. **(verify:** is
  `STATE_PREBORN` used anywhere today?**)**
- [ ] Consider `RungeKutta` with an external field function `(x, y, s, t, params)` as a host
  for space-charge fields. It still steps each particle on its own, so Step 1 still applies.

## Step 7: RF at low β (BeamTracking) **(verify)**

- [ ] `cavity!` takes fixed `beta_0` and `tilde_m` per element. Check whether
  `cavity_kick!` and `SaganCavity` handle large velocity change and phase slip within a
  cavity. Add tests at low β.

## Step 8: SciBmad analysis layer

**Problem.** Twiss, closed orbit, normal forms and dynamic aperture all run single-particle
GTPSA maps. Space charge makes the map depend on the beam distribution. The deposit/FFT path
also can't be pushed through GTPSA as it stands (mesh bounds come from `extrema`, and the
solver is restricted to plain floats).

**Tasks.**
- [ ] Add a frozen analytic (e.g. Gaussian) space-charge kick that GTPSA can differentiate
  through, for optics with space charge, or a self-consistent envelope model.
- [ ] Extend `TrackingConfig` / `_track` (`SciBmad.jl/src/track.jl`) with space-charge
  options (grid size, steps per element, boundary type, solver choice) and pass them down
  to `track!`.

## Step 9: Bug fix found during the survey (BeamTracking)

- [ ] Fix `R_to_v` in `BeamTracking.jl/src/utils/energy.jl:9`.
  - Current: `|q|·c / sqrt(1 + (m/(R·c))²)`
  - Correct: `c / sqrt(1 + (m/(R·q·c))²)`
  - Checked in the REPL: for C⁶⁺ at R = 1 T·m, `R_to_v` gives 4.824e7 m/s, while
    `beta_gamma_to_v` gives 4.764e7 m/s. At high energy `R_to_v` tends to |q|·c, which is
    faster than light.
  - Impact: the harmonic-number RF frequency (`rf_omega_calc` in
    `ext/BeamTrackingBeamlinesExt/utils_bl.jl`), which matters for low-energy ion rings.
  - Add a test for |q| ≠ 1.

---

## Suggested order

1. Step 9: a small, independent bug fix.
2. Steps 1 and 4: the synchronization point plus a kick kernel. Validate with a frozen
   analytic Gaussian field first.
3. Step 2, then hook up SpaceCharge.jl free-space solves. Benchmark against analytic results
   (e.g. KV or Gaussian beam expansion in a drift).
4. Steps 3 and 5: energy binning and boundaries.
5. Steps 7 and 8: RF at low β and the analysis layer.
6. Step 6: only if the scope includes guns.
