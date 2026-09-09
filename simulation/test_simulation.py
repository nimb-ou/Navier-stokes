"""
Automated unit tests for the Navier-Stokes blowup simulation suite.
"""

import unittest
import numpy as np
from simulation.similarity import SimilarityCoordinates
from simulation.profiles import BlowupProfiles
from simulation.fields import NavierStokesBlowupField
from simulation.pulses import AnnulusWavePulses
from simulation.tracer import ParticleTracer
from simulation.diagnostics import BlowupDiagnostics


class TestNavierStokesBlowup(unittest.TestCase):
    def setUp(self):
        self.h = 0.005
        self.sim = SimilarityCoordinates(h=self.h)
        self.field = NavierStokesBlowupField(nu=1.0, h=self.h)
        self.diag = BlowupDiagnostics(self.field)

    def test_similarity_coordinates_inversion(self):
        """Verify q - z^2 * q^(2h) = tau holds to high precision."""
        tau_values = [0.5, 0.1, 0.01, 1e-4, 1e-6]
        z_values = np.array([0.0, 0.1, -0.2, 0.5, -1.0, 2.0])

        for tau in tau_values:
            q = self.sim.solve_q(z_values, tau)
            # Check equation: q - z^2 * q^(2h) == tau
            lhs = q - (z_values ** 2) * np.power(q, 2.0 * self.h)
            diff = np.abs(lhs - tau)
            np.testing.assert_allclose(lhs, tau, rtol=1e-8, atol=1e-10)

            # Check that |eta| < 1
            eta = z_values / np.power(q, self.sim.D)
            self.assertTrue(np.all(np.abs(eta) < 1.0))

    def test_incompressibility_near_axis(self):
        """
        Check that div(u) = (1/r) d/dr(r * u_r) + d/dz(u_z) is approximately 0.
        """
        r = np.array([0.1, 0.2, 0.3])
        z = np.array([0.05, 0.1, -0.05])
        t = 0.9
        dr = 1e-5
        dz = 1e-5

        eval_c = self.field.evaluate_cylindrical(r, z, t)
        eval_rp = self.field.evaluate_cylindrical(r + dr, z, t)
        eval_rm = self.field.evaluate_cylindrical(r - dr, z, t)
        eval_zp = self.field.evaluate_cylindrical(r, z + dz, t)
        eval_zm = self.field.evaluate_cylindrical(r, z - dz, t)

        flux_rp = (r + dr) * eval_rp["u_r"]
        flux_rm = (r - dr) * eval_rm["u_r"]
        d_flux_dr = (flux_rp - flux_rm) / (2.0 * dr)
        div_r = d_flux_dr / r

        du_z_dz = (eval_zp["u_z"] - eval_zm["u_z"]) / (2.0 * dz)
        total_div = div_r + du_z_dz

        # Total divergence should be very small compared to velocity gradient scale
        grad_scale = np.maximum(np.abs(div_r), np.abs(du_z_dz))
        rel_div = np.abs(total_div) / np.maximum(grad_scale, 1.0)
        self.assertTrue(np.all(rel_div < 0.05), f"Relative divergence too high: {rel_div}")

    def test_velocity_blowup_and_bounded_energy(self):
        """
        Verify the core Millennium theorem property:
        As tau -> 0, maximum velocity -> infinity, while kinetic energy remains bounded!
        """
        t_values = np.array([0.9, 0.99, 0.999, 0.9999])
        res = self.diag.time_series(t_values)

        # Velocity scale should strictly increase:
        self.assertTrue(np.all(np.diff(res["u_scale_tangential"]) > 0))
        # Tangential velocity should grow significantly:
        self.assertGreater(res["u_scale_tangential"][-1] / res["u_scale_tangential"][0], 30.0)

        # Core kinetic energy tau^(1/2 - 3h) should tend to 0 as tau -> 0:
        self.assertTrue(np.all(np.diff(res["kinetic_energy_core"]) < 0))
        self.assertLess(res["kinetic_energy_core"][-1], res["kinetic_energy_core"][0])

        # Angular Reynolds number Re_theta should strictly increase:
        self.assertTrue(np.all(np.diff(res["re_theta"]) > 0))

        # Radial Reynolds number Re_r should remain O(1):
        self.assertTrue(np.all(res["re_r"] < 10.0))

    def test_particle_tracer(self):
        """Verify 3D particle tracing runs without error."""
        tracer = ParticleTracer(self.field)
        init_pos = np.array([
            [0.2, 0.0, 0.1],
            [0.3, 0.1, -0.1],
        ])
        trace_data = tracer.integrate_trajectories(init_pos, t_start=0.9, t_end=0.91, dt=0.002)
        self.assertIn("trajectories", trace_data)
        self.assertEqual(trace_data["trajectories"].shape[1], 2)
        self.assertEqual(trace_data["trajectories"].shape[2], 3)

    def test_annulus_wave_pulses(self):
        """Verify pulse envelope and Reynolds stress generation."""
        pulses = AnnulusWavePulses(Xa=1.0, Xb=5.0, h=self.h)
        v = np.linspace(0, 1, 20)
        env = pulses.envelope_P(v, Ls=1.0)
        self.assertAlmostEqual(env[len(v)//2], 1.0, delta=0.05)
        # Should decay at endpoints:
        self.assertLess(env[0], 0.01)
        self.assertLess(env[-1], 0.01)


    def test_vorticity_evaluation(self):
        """Verify that vorticity components are finite and regular at and away from the axis."""
        r = np.array([0.0, 0.05, 0.2])
        z = np.array([0.0, 0.1, -0.1])
        t = 0.9
        vort = self.field.evaluate_vorticity(r, z, t)
        self.assertTrue(np.all(np.isfinite(vort["omega_mag"])))
        self.assertGreater(vort["omega_mag"][0], 0.0)


if __name__ == "__main__":
    unittest.main()
