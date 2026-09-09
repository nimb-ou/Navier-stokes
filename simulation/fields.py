"""
3D Physical Velocity, Pressure, and Vorticity Fields
Transforms similarity variables into full 3D Cartesian and cylindrical fields
for any viscosity nu > 0 and time t in [0, 1).
"""

import numpy as np
from .similarity import SimilarityCoordinates
from .profiles import BlowupProfiles


class NavierStokesBlowupField:
    def __init__(self, nu: float = 1.0, h: float = 0.005, C: float = 2.0, P_star: float = 3.0,
                 Xa: float = 1.0, Xb: float = 5.0, c_inf: float = 1.5, j0: float = 0.03):
        """
        nu: kinematic viscosity > 0
        h: core scaling exponent (0 < h < 0.01)
        """
        assert nu > 0, "Viscosity nu must be strictly positive"
        self.nu = nu
        self.sqrt_nu = np.sqrt(nu)
        self.h = h
        self.sim = SimilarityCoordinates(h=h)
        self.profiles = BlowupProfiles(h=h, C=C, P_star=P_star, Xa=Xa, Xb=Xb, c_inf=c_inf, j0=j0)

    def evaluate_cylindrical(self, r: np.ndarray, z: np.ndarray, t: float):
        """
        Evaluates cylindrical velocity components (u_r, u_theta, u_z) and pressure p.
        Includes viscosity-scaling: u_nu(r, z, t) = sqrt(nu) * u_1(r / sqrt(nu), z / sqrt(nu), t).
        """
        r_arr = np.asarray(r, dtype=np.float64)
        z_arr = np.asarray(z, dtype=np.float64)

        # Scale to viscosity 1 reference coordinates:
        r_ref = r_arr / self.sqrt_nu
        z_ref = z_arr / self.sqrt_nu

        sim_vars = self.sim.compute_all(r_ref, z_ref, t)
        q = sim_vars["q"]
        eta = sim_vars["eta"]
        X = sim_vars["X"]
        d = sim_vars["d"]
        L = sim_vars["L"]
        tau = sim_vars["tau"]

        q_A = np.power(q, self.profiles.A)
        q_2A = np.power(q, 2.0 * self.profiles.A)

        # Profile functions:
        E = self.profiles.E_profile(X, eta)
        U = self.profiles.U_profile(X, eta)
        V0 = self.profiles.V0_flux(X, eta, d, L)
        Pi = self.profiles.Pi_pressure(X, eta)

        # Cylindrical velocities at nu = 1:
        u_theta_ref = E / q_A
        u_z_ref = U / q_A

        # u_r = V0 / r_ref:
        # Near r = 0, V0 / r = (r / 2qL) * [...] which is regular and O(r)
        safe_r = np.where(r_ref > 1e-12, r_ref, 1e-12)
        u_r_ref = np.where(r_ref > 1e-12, V0 / safe_r, 0.0)

        p_ref = Pi / q_2A

        # Rescale by viscosity:
        u_r = self.sqrt_nu * u_r_ref
        u_theta = self.sqrt_nu * u_theta_ref
        u_z = self.sqrt_nu * u_z_ref
        p = self.nu * p_ref

        # Characteristic core scales:
        ell_r = self.sqrt_nu * sim_vars["ell_r"]
        ell_z = self.sqrt_nu * sim_vars["ell_z"]

        return {
            "u_r": u_r,
            "u_theta": u_theta,
            "u_z": u_z,
            "p": p,
            "q": q,
            "eta": eta,
            "X": X,
            "tau": tau,
            "ell_r": ell_r,
            "ell_z": ell_z,
            "speed": np.sqrt(u_r ** 2 + u_theta ** 2 + u_z ** 2),
        }

    def evaluate_cartesian(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, t: float):
        """
        Evaluates 3D Cartesian velocity (u_x, u_y, u_z) and pressure p at (x, y, z, t).
        """
        x_arr = np.asarray(x, dtype=np.float64)
        y_arr = np.asarray(y, dtype=np.float64)
        z_arr = np.asarray(z, dtype=np.float64)

        r = np.sqrt(x_arr ** 2 + y_arr ** 2)
        cyl = self.evaluate_cylindrical(r, z_arr, t)

        safe_r = np.where(r > 1e-14, r, 1.0)
        cos_theta = np.where(r > 1e-14, x_arr / safe_r, 1.0)
        sin_theta = np.where(r > 1e-14, y_arr / safe_r, 0.0)

        u_r = cyl["u_r"]
        u_theta = cyl["u_theta"]

        u_x = u_r * cos_theta - u_theta * sin_theta
        u_y = u_r * sin_theta + u_theta * cos_theta
        u_z = cyl["u_z"]

        return {
            "u_x": u_x,
            "u_y": u_y,
            "u_z": u_z,
            "speed": np.sqrt(u_x ** 2 + u_y ** 2 + u_z ** 2),
            "p": cyl["p"],
            "u_r": u_r,
            "u_theta": u_theta,
            "r": r,
            "tau": cyl["tau"],
            "ell_r": cyl["ell_r"],
            "ell_z": cyl["ell_z"],
        }

    def evaluate_vorticity(self, r: np.ndarray, z: np.ndarray, t: float, dr: float = 1e-4, dz: float = 1e-4):
        """
        Computes 3D vorticity omega = curl(u) using central finite differences in cylindrical coordinates:
          omega_r = -du_theta / dz
          omega_theta = du_r / dz - du_z / dr
          omega_z = (1/r) * d(r * u_theta) / dr
        """
        r_arr = np.asarray(r, dtype=np.float64)
        z_arr = np.asarray(z, dtype=np.float64)

        # Shifted evaluations for gradients
        cyl_c = self.evaluate_cylindrical(r_arr, z_arr, t)
        cyl_rp = self.evaluate_cylindrical(r_arr + dr, z_arr, t)
        cyl_rm = self.evaluate_cylindrical(np.maximum(r_arr - dr, 0.0), z_arr, t)
        cyl_zp = self.evaluate_cylindrical(r_arr, z_arr + dz, t)
        cyl_zm = self.evaluate_cylindrical(r_arr, z_arr - dz, t)

        du_theta_dz = (cyl_zp["u_theta"] - cyl_zm["u_theta"]) / (2.0 * dz)
        du_r_dz = (cyl_zp["u_r"] - cyl_zm["u_r"]) / (2.0 * dz)
        du_z_dr = (cyl_rp["u_z"] - cyl_rm["u_z"]) / (2.0 * dr)

        # d(r * u_theta) / dr:
        flux_p = (r_arr + dr) * cyl_rp["u_theta"]
        flux_m = np.maximum(r_arr - dr, 0.0) * cyl_rm["u_theta"]
        d_flux_dr = (flux_p - flux_m) / (2.0 * dr)

        safe_r = np.where(r_arr > 1e-12, r_arr, 1.0)
        omega_r = -du_theta_dz
        omega_theta = du_r_dz - du_z_dr
        omega_z = np.where(r_arr > 1e-12, d_flux_dr / safe_r, 2.0 * cyl_rp["u_theta"] / dr)

        omega_mag = np.sqrt(omega_r ** 2 + omega_theta ** 2 + omega_z ** 2)

        return {
            "omega_r": omega_r,
            "omega_theta": omega_theta,
            "omega_z": omega_z,
            "omega_mag": omega_mag,
        }
