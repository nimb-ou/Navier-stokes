"""
Similarity Coordinates Solver
Implements the exact similarity transformations from Section 3.1 and Section 4.1
of "Finite Time Blowup for Navier-Stokes" (OpenAI).

Coordinates:
  tau = 1 - t
  A = 1/2 + h
  D = 1/2 - h
  q - z^2 * q^(2h) = tau
  eta = z / q^D
  X = r^2 / (2 * q)
"""

import numpy as np


class SimilarityCoordinates:
    def __init__(self, h: float = 0.005):
        """
        Initialize similarity coordinates.
        h: parameter in (0, 1/100), default 0.005.
        """
        assert 0 < h < 0.01, "h must be in (0, 1/100)"
        self.h = h
        self.A = 0.5 + h
        self.D = 0.5 - h

    def solve_q(self, z: np.ndarray, tau: float, tol: float = 1e-12, max_iter: int = 50) -> np.ndarray:
        """
        Solves q - z^2 * q^(2h) = tau for q > |z|^(1/D) via Newton-Raphson.
        Guaranteed unique root because d/dq(q - z^2 q^(2h)) = 1 - 2h eta^2 >= 1 - 2h > 0.
        """
        z_arr = np.asarray(z, dtype=np.float64)
        tau = float(tau)
        if tau <= 0:
            raise ValueError(f"tau must be strictly positive (tau = {tau})")

        # Initial guess:
        # If z == 0, q = tau.
        # In general, q >= tau and q > |z|^(1/D).
        z_abs = np.abs(z_arr)
        q_min = np.power(np.maximum(z_abs, 1e-15), 1.0 / self.D)
        q = np.maximum(tau + q_min, tau)

        for _ in range(max_iter):
            q_2h = np.power(q, 2.0 * self.h)
            F = q - (z_arr ** 2) * q_2h - tau
            F_prime = 1.0 - 2.0 * self.h * (z_arr ** 2) * (q_2h / q)
            F_prime = np.maximum(F_prime, 1.0 - 2.0 * self.h)
            delta_q = F / F_prime
            q = q - delta_q
            q = np.maximum(q, q_min + 1e-14 * tau)

            if np.all(np.abs(delta_q) < tol * (1.0 + np.abs(q))):
                break

        return q

    def compute_all(self, r: np.ndarray, z: np.ndarray, t: float):
        """
        Computes all similarity variables (q, eta, X, d, L, s) from physical (r, z, t).
        """
        tau = 1.0 - t
        r_arr = np.asarray(r, dtype=np.float64)
        z_arr = np.asarray(z, dtype=np.float64)

        q = self.solve_q(z_arr, tau)
        q_D = np.power(q, self.D)
        eta = z_arr / q_D
        eta = np.clip(eta, -0.999999999999, 0.999999999999)

        s = 0.5 * (r_arr ** 2)
        X = s / q
        d = 1.0 - (eta ** 2)
        L = 1.0 - 2.0 * self.h * (eta ** 2)

        return {
            "tau": tau,
            "q": q,
            "eta": eta,
            "X": X,
            "s": s,
            "d": d,
            "L": L,
            "ell_r": np.sqrt(q),
            "ell_z": q_D,
        }
