"""
Edge-Case, Boundary Regularity, and Scale-Invariance Tests for Navier-Stokes Simulation Suite
Validates:
1. Axisymmetric smoothness at r = 0 (no division-by-zero or coordinate singularities).
2. Viscosity scaling invariance u_nu(r, z, t) = sqrt(nu) * u_1(r/sqrt(nu), z/sqrt(nu), t).
3. Extreme coordinate regimes (|z| >> 1 and deep blowup tau -> 0).
4. Adaptive particle tracer robustness for on-axis and far-field particles.
5. Incompressible velocity and pressure positivity and finiteness.
"""

import unittest
import numpy as np
from simulation.similarity import SimilarityCoordinates
from simulation.profiles import BlowupProfiles
from simulation.fields import NavierStokesBlowupField
from simulation.tracer import ParticleTracer
from simulation.diagnostics import BlowupDiagnostics


class TestNavierStokesEdgeCases(unittest.TestCase):
    def setUp(self):
        self.h = 0.005
        self.field = NavierStokesBlowupField(nu=1.0, h=self.h)
        self.sim = SimilarityCoordinates(h=self.h)

    def test_regularity_at_cylindrical_axis(self):
        """Verify that at r = 0, u_r = 0, u_theta = 0, and Cartesian fields are regular."""
        r_zero = np.array([0.0, 0.0, 0.0])
        z = np.array([-0.5, 0.0, 0.5])
        t = 0.95

        cyl = self.field.evaluate_cylindrical(r_zero, z, t)
        np.testing.assert_allclose(cyl["u_r"], 0.0, atol=1e-12)
        # u_theta at r=0 has E(0, eta) = 0 because profile is O(X) near X=0
        np.testing.assert_allclose(cyl["u_theta"], 0.0, atol=1e-6)
        self.assertTrue(np.all(np.isfinite(cyl["u_z"])))
        self.assertTrue(np.all(np.isfinite(cyl["p"])))

        # In Cartesian coordinates at (0, 0, z):
        cart = self.field.evaluate_cartesian(np.zeros(3), np.zeros(3), z, t)
        np.testing.assert_allclose(cart["u_x"], 0.0, atol=1e-6)
        np.testing.assert_allclose(cart["u_y"], 0.0, atol=1e-6)
        self.assertTrue(np.all(np.isfinite(cart["u_z"])))

    def test_viscosity_scaling_identity(self):
        """
        Verify the exact Navier-Stokes Navier-Leray invariance:
        u_nu(r, z, t) = sqrt(nu) * u_1(r / sqrt(nu), z / sqrt(nu), t)
        for multiple viscosities nu in {0.1, 2.5, 10.0}.
        """
        r = np.array([0.1, 0.2, 0.5])
        z = np.array([0.0, 0.1, -0.2])
        t = 0.9

        eval_base = self.field.evaluate_cylindrical(r, z, t)

        for nu in [0.1, 0.5, 2.5, 10.0]:
            field_nu = NavierStokesBlowupField(nu=nu, h=self.h)
            sqrt_nu = np.sqrt(nu)

            # Coordinates scaled by sqrt(nu)
            r_scaled = r * sqrt_nu
            z_scaled = z * sqrt_nu

            eval_nu = field_nu.evaluate_cylindrical(r_scaled, z_scaled, t)

            # u_nu(r*sqrt(nu), z*sqrt(nu), t) should equal sqrt(nu) * u_1(r, z, t)
            np.testing.assert_allclose(
                eval_nu["u_theta"],
                sqrt_nu * eval_base["u_theta"],
                rtol=1e-7,
                atol=1e-9
            )
            np.testing.assert_allclose(
                eval_nu["u_z"],
                sqrt_nu * eval_base["u_z"],
                rtol=1e-7,
                atol=1e-9
            )
            np.testing.assert_allclose(
                eval_nu["u_r"],
                sqrt_nu * eval_base["u_r"],
                rtol=1e-7,
                atol=1e-9
            )

    def test_extreme_similarity_coordinate_inversion(self):
        """
        Verify that similarity coordinates solve accurately for extreme z values
        (|z| up to 50.0) and deep blowup tau (1e-8).
        """
        tau_extreme = 1e-8
        z_extreme = np.array([-50.0, -10.0, 0.0, 10.0, 50.0])

        q = self.sim.solve_q(z_extreme, tau_extreme)
        # Verify q > 0 strictly
        self.assertTrue(np.all(q > 0.0))

        # Check equation: q - z^2 * q^(2h) == tau
        lhs = q - (z_extreme ** 2) * np.power(q, 2.0 * self.h)
        np.testing.assert_allclose(lhs, tau_extreme, rtol=1e-7, atol=1e-9)

        # Coordinate eta = z / q^D must satisfy |eta| < 1
        eta = z_extreme / np.power(q, self.sim.D)
        self.assertTrue(np.all(np.abs(eta) < 1.0))

    def test_particle_tracer_on_axis_and_extreme_points(self):
        """Verify particle tracer on axis and with adaptive steps without crashing or NaNs."""
        tracer = ParticleTracer(self.field)
        positions = np.array([
            [0.0, 0.0, 0.05],    # Exactly on axis
            [0.1, -0.1, 0.2],    # Off axis
            [1.0, 1.0, -1.0],    # Far field
        ])
        trace_data = tracer.integrate_trajectories(
            positions,
            t_start=0.90,
            t_end=0.91,
            dt=0.005
        )

        self.assertIn("trajectories", trace_data)
        traj = trace_data["trajectories"]
        self.assertFalse(np.any(np.isnan(traj)))
        self.assertFalse(np.any(np.isinf(traj)))
        self.assertEqual(traj.shape[1], 3)  # 3 particles

    def test_divergence_free_property_multiple_regimes(self):
        """Verify divergence div(u) remains bounded and close to zero across multiple time slices."""
        r = np.array([0.15, 0.25])
        z = np.array([0.05, -0.1])
        dr = 1e-5
        dz = 1e-5

        for t in [0.5, 0.8, 0.95]:
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

            grad_scale = np.maximum(np.abs(div_r), np.abs(du_z_dz))
            rel_div = np.abs(total_div) / np.maximum(grad_scale, 1.0)
            self.assertTrue(np.all(rel_div < 0.05))


if __name__ == "__main__":
    unittest.main()
