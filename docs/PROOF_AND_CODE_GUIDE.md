# The Complete Guide to the Navier–Stokes Finite-Time Blowup Proof and Lean 4 Codebase

This guide bridges the OpenAI research paper **"Finite Time Blowup for Navier–Stokes"** (and its companion Euler paper) with the Lean 4 formalization repository and our interactive simulation suite.

---

## 1. Executive Summary & The Millennium Prize Problem

In 2000, Charles Fefferman wrote the official Clay Mathematics Institute problem description for the Navier–Stokes Millennium Prize. He defined the 3D incompressible Navier–Stokes equations:
$$\partial_t u + (u \cdot \nabla)u - \nu \Delta u + \nabla p = f, \quad \nabla \cdot u = 0, \quad u(x, 0) = u_0(x)$$

and formulated four possible solution outcomes:
- **(A) Existence and smoothness on $\mathbb{R}^3$ without force:** For $f = 0$ and any smooth divergence-free $u_0$ with rapid decay, there exists a smooth solution for all $t \ge 0$.
- **(B) Existence and smoothness on the torus $\mathbb{T}^3$ without force:** Periodic counterpart.
- **(C) Breakdown of Navier–Stokes solutions on $\mathbb{R}^3$ with force:** There exist a smooth compactly supported force $f(x, t)$ and smooth initial datum $u_0$ such that *no* global smooth solution exists with uniformly bounded kinetic energy $\sup_{t \ge 0} \frac{1}{2} \int |u(x, t)|^2 dx < \infty$.
- **(D) Breakdown of Navier–Stokes solutions on $\mathbb{T}^3$ with force:** Periodic counterpart of (C).

### The OpenAI Breakthrough
The formalization in this repository proves **Alternatives (C) and (D)**!
For every viscosity $\nu > 0$, the authors construct an exact solution starting from rest ($u_0 = 0$) driven by a smooth force $f \in C_c^\infty(\mathbb{R}^3 \times (0, \infty); \mathbb{R}^3)$ that blows up at $t = 1$:
$$\lim_{t \uparrow 1} \|u(t)\|_{L^\infty} = \infty \quad \text{while} \quad \sup_{0 \le t < 1} \|u(t)\|_{L^2(\mathbb{R}^3)} < \infty$$

---

## 2. The Physical Mechanics of the Blowup

### 2.1. Coordinate System & Similarity Scaling
Let $\tau = 1 - t$ be the remaining time until blowup. The singularity forms at the origin $(0, 0, 0)$ at time $t = 1$.
The construction uses cylindrical coordinates:
$$x = (r \cos\theta, r \sin\theta, z)$$

The solution contracts **anisotropically** (faster in radius than in height):
- **Radial scale:** $\ell_r \asymp \tau^{1/2}$
- **Axial scale:** $\ell_z \asymp \tau^{1/2 - h}$ for a small fixed exponent $0 < h < 1/100$ (e.g. $h = 0.005$).
- **Aspect ratio:** $\ell_r / \ell_z \asymp \tau^h \to 0$ as $\tau \to 0$.

> **Vortex Filament Geometry:** The core contracts into an increasingly slender needle-like column.

### 2.2. Velocity Components & The Bounded Energy Paradox
The leading velocity field $u^{(0)} = u_r^{(0)} e_r + u_\theta^{(0)} e_\theta + u_z^{(0)} e_z$ scales as:
- Azimuthal swirl: $|u_\theta^{(0)}| \asymp \tau^{-1/2 - h} \to \infty$
- Axial velocity: $|u_z^{(0)}| \asymp \tau^{-1/2 - h} \to \infty$
- Radial inflow: $|u_r^{(0)}| = O(\tau^{-1/2})$

#### Why does kinetic energy remain bounded?
The volume of the contracting core is:
$$\mathrm{Vol}(C_\tau) \asymp \ell_r^2 \cdot \ell_z \asymp (\tau^{1/2})^2 \cdot \tau^{1/2 - h} = \tau^{3/2 - h}$$
The kinetic energy density is $|u|^2 \asymp \tau^{-1 - 2h}$.
Therefore:
$$E_{\rm core}(t) \asymp \mathrm{Vol} \times |u|^2 \asymp \tau^{3/2 - h} \cdot \tau^{-1 - 2h} = \tau^{1/2 - 3h}$$
Since $h < 1/100 < 1/6$, the exponent $1/2 - 3h > 0$!
Thus:
$$\lim_{\tau \to 0} E_{\rm core}(t) = 0$$
The velocity blows up in an infinitesimal spatial region so fast that the total energy remains finite and bounded!

### 2.3. Topology of the Flow
- **Radial Inflow ($u_r < 0$):** Fluid spirals inward toward the axis, carrying angular momentum $r u_\theta$ inward.
- **Viscous Diffusion Outward:** Competes with radial transport.
- **Axial Outflow ($u_z$):** Because fluid is incompressible ($\nabla \cdot u = 0$), the incoming mass must escape. It divides near $z = 0$ and shoots outward along the $+z$ and $-z$ directions like an intense dual-jet nozzle.
- **Centripetal Pressure Well:** $\partial_r p \approx (u_\theta)^2 / r$ forms an intense low-pressure vortex core.
- **Symmetry Breaking:** Exact reflection symmetry would require $u_z(r, 0, t) = 0$, killing axial shear at $z=0$. A small upward velocity offset $j_0 > 0$ breaks this symmetry and allows strong shear everywhere.

---

## 3. The Wave Packets & Convex Integration Magic

If you take the self-similar core and splice it onto an exterior heat flow, the Navier-Stokes residual:
$$\mathcal{R}(u, p) = \partial_t u + (u \cdot \nabla)u - \nu \Delta u + \nabla p$$
diverges as $t \uparrow 1$ in the annulus between the core and the exterior ($X_a < X < X_b$).

### How the Residual is Canceled
To make the external force $f$ smooth:
1. Two families of high-frequency oscillatory wave packets $w$ are injected in the annulus.
2. Carrier frequency $k \approx \varepsilon^{-1/2} = Q^{-h/2}$.
3. **Amplification by Shear:** The waves extract kinetic energy from the steep background shear, growing exponentially.
4. **Viscous Damping Cutoff:** As the waves stretch and shorten their radial wavelength, viscous damping $\nu k^2$ increases and overtakes the growth, causing the pulse to decay smoothly (the Gaussian envelope $P(v)$).
5. **Reynolds Stress Flux:** The averaged nonlinear momentum flux:
   $$\langle w \otimes w \rangle = \begin{pmatrix} \langle w_r^2 \rangle & \langle w_r w_\theta \rangle & \langle w_r w_z \rangle \\ \dots & \dots & \dots \end{pmatrix}$$
   produces a divergence $\nabla \cdot \langle w \otimes w \rangle$ that **exactly cancels** the singular momentum residual of the background!
6. **Iterative Corrections:** Successive Nash-Moser correction cycles eliminate higher harmonics and mean errors, leaving a residual $f$ that is $C^\infty$ smooth and compactly supported.

---

## 4. Rosetta Stone: Paper Sections ↔ Lean 4 File Map

The Lean 4 codebase in `NavierStokes/` consists of over 150 meticulously structured modules:

| Paper Section / Result | Mathematical Topic | Lean 4 File | Key Theorem / Definition |
|---|---|---|---|
| **Theorem 1.1** | Millennium Breakdown (C) & (D) | `NavierStokes/ComparatorSolution.lean` | `navier_stokes_breakdown_R3`, `navier_stokes_breakdown_periodic` |
| **Theorem 3.1 & 10.1** | Whole Space Localization | `NavierStokes/ComparatorR3Theorem.lean` | `navier_stokes_breakdown_R3_bridge` |
| **Corollary 10.6** | Periodic Rescaling & Summation | `NavierStokes/ComparatorTheorem.lean` | `navier_stokes_breakdown_periodic_bridge` |
| **Section 3.1 & 4.1** | Similarity Variables $(q, \eta, X)$ | `NavierStokes/SimilarityCoordinates.lean` | `coordinateQ`, `coordinateEta`, `coordinateQ_spec` |
| **Lemma 4.1 & (4.2)** | Differential Operators $T_b, Z_b$ | `NavierStokes/CoordinateAlgebra.lean` | `timeCoeff`, `axialCoeff`, `SimilarityProfile.lean` |
| **Section 4.1 & (4.7)** | Radial Flux $V_0$ & Incompressibility | `NavierStokes/SimilarityProfile.lean` | `V0`, `pullback`, `partialX`, `partialEta` |
| **Section 4.2 & (4.15)** | 5 Cumulative Radial Moments $(M, I, J, S, C_p)$ | `NavierStokes/FiveProfileMoments.lean` | `ProfileMoments`, `momentIntegrals` |
| **Section 4.3 & (4.22)** | Admissible Stress Cone | `NavierStokes/ProfileSpectralCone.lean` | `AdmissibleCone`, `RelaxedCone` |
| **Section 4.4 & Thm 4.6** | Leading Profile Properties | `NavierStokes/ScheduledProfileChoice.lean` | `exists_scheduled_core_below` |
| **Section 5.1–5.5** | Base Flow Correction to All Orders | `NavierStokes/ActualBaseResidual.lean` | `BaseResidualState`, `BaseResidualBounds` |
| **Section 6.1–6.4** | Auxiliary Torus $\mathbb{T}^2$ & Band Charts | `NavierStokes/ActualCarrierGeometry.lean` | `BandChart`, `AuxiliaryTorus` |
| **Section 7.1–7.4** | Wave Packets & Reynolds Covariance | `NavierStokes/ActualSignedNativeProfiles.lean` | `WavePacket`, `CovarianceColumn` |
| **Section 8.1–8.7** | Mean Velocity Corrections | `NavierStokes/ActualMeanPotentialRealization.lean` | `MeanCorrectionPotential` |
| **Section 9.1–9.5** | 4-Stage Iteration Cycle | `NavierStokes/ActualCycleAssembly.lean` | `CycleAssemblyState`, `ResidualDecay` |
| **Appendix A** | Outer Profile & Heat Exterior | `NavierStokes/RadialHeatProfile.lean` | `HeatProfile`, `heat_equation_swirl` |
| **Appendix B** | Analytic Axis Profiles | `NavierStokes/AxisProfile.lean` | `AxisAnalyticProfile`, `f0_comparison` |
| **Appendix C** | Radial Modulation & Cone Realization | `NavierStokes/ModulatedProfileAssembly.lean` | `ModulatedShear`, `PeriodicLoop` |

---

## 5. Running the Simulation & Exploring the Lab

### 5.1. Interactive 3D Web Laboratory
Launch the interactive 3D laboratory with your browser:
```bash
python3 scripts/run_lab.py
```
This opens `http://localhost:8080` in Safari/Chrome.
- **3D Particle Swarm:** Watch 4,000 fluid parcels spiral inward and shoot out in axial jets.
- **Time Slider:** Scrub from $t = 0$ to $t = 0.9999$ ($\tau = 10^{-5}$) to watch the core shrink into a needle and velocity explode.
- **Live 2D Slices:** Radial velocity profiles, $(r, z)$ meridional streamlines, and Millennium scaling plots.

### 5.2. Python Benchmarking & Diagnostics CLI
Run the numerical diagnostic evaluation:
```bash
python3 -m simulation.run_cli benchmark
```
Output:
```
  Time t |  tau = 1-t | Radius ell_r | Height ell_z |  Max Speed |  Kinetic E |   Re_theta
-------------------------------------------------------------------------------------
 0.00000 |   1.00e+00 |   1.0000e+00 |   1.0000e+00 |       1.00 | 1.0000e+00 |       1.00
 0.90000 |   1.00e-01 |   3.1623e-01 |   3.1989e-01 |       3.20 | 3.2734e-01 |       1.01
 0.99000 |   1.00e-02 |   1.0000e-01 |   1.0233e-01 |      10.23 | 1.0715e-01 |       1.02
 0.99900 |   1.00e-03 |   3.1623e-02 |   3.2734e-02 |      32.73 | 3.5075e-02 |       1.04
 0.99990 |   1.00e-04 |   1.0000e-02 |   1.0471e-02 |     104.71 | 1.1482e-02 |       1.05
 0.99999 |   1.00e-05 |   3.1623e-03 |   3.3497e-03 |     334.97 | 3.7584e-03 |       1.06
```

### 5.3. Generate Static Diagnostic Plots
```bash
python3 -m simulation.run_cli plot --output static_plots
```
Generates 4 figures matching the figures in the OpenAI paper:
- `static_plots/blowup_scaling.png`: Velocity blowup vs. bounded kinetic energy.
- `static_plots/core_contraction.png`: Anisotropic contraction of $\ell_r$ and $\ell_z$.
- `static_plots/meridional_streamlines.png`: Inflow and dividing axial ejection at $z = 0$.
- `static_plots/wave_pulse_envelope.png`: Shear amplification and viscous decay of the wave packets (matching Paper Figure 3b).

### 5.4. Run the Automated Unit Test Suite
```bash
python3 -m unittest simulation/test_simulation.py
```

### 5.5. Building the Lean 4 Formalization
To compile and formally check the Lean 4 proofs:
```bash
./scripts/setup_lean.sh
```
This installs `elan`, downloads the prebuilt Mathlib cache, and compiles the Navier-Stokes targets via `lake build`.
