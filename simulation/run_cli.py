"""
Command-Line Interface for Navier-Stokes Blowup Analysis & Plot Generation
Usage:
  python3 -m simulation.run_cli benchmark
  python3 -m simulation.run_cli plot --output static_plots
"""

import os
import argparse
import numpy as np

# Set writable cache directory for Matplotlib to avoid macOS permission warnings
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulation.similarity import SimilarityCoordinates
from simulation.fields import NavierStokesBlowupField
from simulation.pulses import AnnulusWavePulses
from simulation.tracer import ParticleTracer
from simulation.diagnostics import BlowupDiagnostics


def run_benchmark():
    print("=" * 70)
    print("Navier-Stokes Finite-Time Blowup: Physical & Mathematical Scaling")
    print("Reference: 'Finite Time Blowup for Navier-Stokes' (OpenAI)")
    print("=" * 70)

    field = NavierStokesBlowupField(nu=1.0, h=0.005)
    diag = BlowupDiagnostics(field)

    t_vals = [0.0, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999]
    print(f"{'Time t':>8} | {'tau = 1-t':>10} | {'Radius ell_r':>12} | {'Height ell_z':>12} | {'Max Speed':>10} | {'Kinetic E':>10} | {'Re_theta':>10}")
    print("-" * 85)

    for t in t_vals:
        res = diag.compute_scaling_at_time(t)
        print(f"{res['t']:8.5f} | {res['tau']:10.2e} | {res['ell_r']:12.4e} | {res['ell_z']:12.4e} | "
              f"{res['u_scale_tangential']:10.2f} | {res['kinetic_energy_core']:10.4e} | {res['re_theta']:10.2f}")

    print("=" * 70)
    print("Key Highlights:")
    print("1. As tau -> 0, Velocity scale ~ tau^(-1/2-h) DIVERGES to infinity.")
    print("2. Core Kinetic Energy ~ tau^(1/2-3h) TENDS TO ZERO (bounded!).")
    print("3. Slenderness ratio ell_r / ell_z ~ tau^h -> 0 (forms razor-thin needle).")
    print("4. Angular Reynolds number Re_theta -> infinity (super-fast rotation).")
    print("=" * 70)


def generate_plots(output_dir: str = "static_plots"):
    os.makedirs(output_dir, exist_ok=True)
    print(f"Generating diagnostic figures in '{output_dir}'...")

    field = NavierStokesBlowupField(nu=1.0, h=0.005)
    diag = BlowupDiagnostics(field)

    # 1. Scaling Plot: Velocity Divergence vs Bounded Kinetic Energy
    tau_vals = np.logspace(-5, 0, 100)
    t_vals = 1.0 - tau_vals
    res = diag.time_series(t_vals)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: Velocity blowup
    ax1.loglog(tau_vals, res["u_scale_tangential"], 'r-', lw=2.5, label=r'$\|u_\theta\|_{L^\infty} \sim \tau^{-1/2-h}$')
    ax1.set_xlabel(r'Time remaining $\tau = 1 - t$', fontsize=12)
    ax1.set_ylabel('Velocity Scale', fontsize=12)
    ax1.set_title('Finite-Time Velocity Blowup (Divergence)', fontsize=13, fontweight='bold')
    ax1.grid(True, which="both", ls="--", alpha=0.5)
    ax1.legend(fontsize=11)

    # Right: Kinetic Energy boundedness
    ax2.loglog(tau_vals, res["kinetic_energy_core"], 'b-', lw=2.5, label=r'$E_{\rm core}(t) \sim \tau^{1/2-3h} \to 0$')
    ax2.set_xlabel(r'Time remaining $\tau = 1 - t$', fontsize=12)
    ax2.set_ylabel('Core Kinetic Energy', fontsize=12)
    ax2.set_title('Kinetic Energy Boundedness (Millennium Prop C)', fontsize=13, fontweight='bold')
    ax2.grid(True, which="both", ls="--", alpha=0.5)
    ax2.legend(fontsize=11)

    plt.tight_layout()
    plot1_path = os.path.join(output_dir, "blowup_scaling.png")
    fig.savefig(plot1_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {plot1_path}")

    # 2. Geometry & Anisotropic Core Contraction (ell_r vs ell_z)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(tau_vals, res["ell_r"], 'g-', lw=2.5, label=r'Radial scale $\ell_r \sim \tau^{1/2}$')
    ax.loglog(tau_vals, res["ell_z"], 'm--', lw=2.5, label=r'Axial scale $\ell_z \sim \tau^{1/2-h}$')
    ax.loglog(tau_vals, res["aspect_ratio"], 'k:', lw=2, label=r'Aspect ratio $\ell_r / \ell_z \sim \tau^h \to 0$')
    ax.set_xlabel(r'Time remaining $\tau = 1 - t$', fontsize=12)
    ax.set_ylabel('Length scale', fontsize=12)
    ax.set_title('Anisotropic Contraction of the Vortex Core', fontsize=13, fontweight='bold')
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plot2_path = os.path.join(output_dir, "core_contraction.png")
    fig.savefig(plot2_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {plot2_path}")

    # 3. Meridional Flow Field (r-z slice streamlines at t = 0.95)
    t_snap = 0.95
    tau_snap = 1.0 - t_snap
    ell_r_snap = np.sqrt(tau_snap)
    ell_z_snap = np.power(tau_snap, 0.495)

    r_grid = np.linspace(0.01 * ell_r_snap, 3.0 * ell_r_snap, 50)
    z_grid = np.linspace(-1.5 * ell_z_snap, 1.5 * ell_z_snap, 50)
    R_grid, Z_grid = np.meshgrid(r_grid, z_grid)

    eval_meridional = field.evaluate_cylindrical(R_grid.ravel(), Z_grid.ravel(), t_snap)
    Ur = eval_meridional["u_r"].reshape(R_grid.shape)
    Uz = eval_meridional["u_z"].reshape(R_grid.shape)
    Utheta = eval_meridional["u_theta"].reshape(R_grid.shape)

    fig, ax = plt.subplots(figsize=(8, 6))
    strm = ax.streamplot(r_grid / ell_r_snap, z_grid / ell_z_snap, Ur, Uz,
                         color=Utheta, cmap='plasma', density=1.5, linewidth=1.2)
    cbar = fig.colorbar(strm.lines, ax=ax)
    cbar.set_label(r'Swirl Velocity $u_\theta$', fontsize=11)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.6, label='Dividing layer (z=0)')
    ax.set_xlabel(r'Normalized radius $r / \ell_r$', fontsize=12)
    ax.set_ylabel(r'Normalized height $z / \ell_z$', fontsize=12)
    ax.set_title(f'Meridional Inflow & Axial Outflow (t = {t_snap})', fontsize=13, fontweight='bold')
    ax.legend()
    plt.tight_layout()
    plot3_path = os.path.join(output_dir, "meridional_streamlines.png")
    fig.savefig(plot3_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {plot3_path}")

    # 4. Wave Packet Envelope P(v) (matching Figure 3b from the paper)
    pulses = AnnulusWavePulses(Xa=1.0, Xb=5.0, h=0.005)
    v_norm = np.linspace(0, 1, 300)
    envelope = pulses.envelope_P(v_norm, Ls=1.0)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(v_norm, envelope, 'k-', lw=2.5, label=r'Reference Amplitude $P(v)$')
    ax.fill_between(v_norm, 0, envelope, where=(v_norm < 0.2), color='orange', alpha=0.3, label='Shear Amplification')
    ax.fill_between(v_norm, 0, envelope, where=(v_norm >= 0.2) & (v_norm <= 0.8), color='green', alpha=0.2, label='Peak Wave Activity')
    ax.fill_between(v_norm, 0, envelope, where=(v_norm > 0.8), color='red', alpha=0.3, label='Viscous Decay')
    ax.set_xlabel('Slot Time coordinate v / Ls', fontsize=12)
    ax.set_ylabel('Normalized Amplitude P(v)', fontsize=12)
    ax.set_title('Annulus Wave Pulse Evolution (Paper Figure 3b)', fontsize=13, fontweight='bold')
    ax.grid(True, ls="--", alpha=0.5)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plot4_path = os.path.join(output_dir, "wave_pulse_envelope.png")
    fig.savefig(plot4_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {plot4_path}")

    print("All diagnostic plots successfully generated!")


def main():
    parser = argparse.ArgumentParser(description="Navier-Stokes Blowup Analysis CLI")
    parser.add_argument("command", choices=["benchmark", "plot"], help="Command to run")
    parser.add_argument("--output", default="static_plots", help="Directory for plots")
    args = parser.parse_args()

    if args.command == "benchmark":
        run_benchmark()
    elif args.command == "plot":
        generate_plots(args.output)


if __name__ == "__main__":
    main()
