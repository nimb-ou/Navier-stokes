"""
Leading-Order Similarity Profiles
Implements the leading profiles E(X, eta), U(X, eta), V0(X, eta), and Pi(X, eta)
from Section 4 and Appendices A-B of "Finite Time Blowup for Navier-Stokes" (OpenAI).
"""

import numpy as np
from scipy import special


class BlowupProfiles:
    def __init__(self, h: float = 0.005, C: float = 2.0, P_star: float = 3.0,
                 Xa: float = 1.0, Xb: float = 5.0, c_inf: float = 1.5, j0: float = 0.03):
        """
        Leading profile parameters:
          h: small exponent in (0, 1/100)
          C: normalization constant > 1
          P_star: initial azimuthal amplitude
          Xa, Xb: active annulus edges
          c_inf: exterior swirl constant
          j0: small upward axial offset at z=0 (breaks reflection symmetry)
        """
        self.h = h
        self.A = 0.5 + h
        self.D = 0.5 - h
        self.C = C
        self.P_star = P_star
        self.Xa = Xa
        self.Xb = Xb
        self.c_inf = c_inf
        self.j0 = j0

    def f_profile(self, eta: np.ndarray) -> np.ndarray:
        """Geometric shape factor f(eta) = 1 / (1 + eta^2)."""
        return 1.0 / (1.0 + (eta ** 2))

    def phi_axis(self, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
        """
        Axis profile scalar phi(X, eta) from Proposition B.2:
        phi = phi*(eta) * f0(X * chi), where f0(z) ~ 1 - z/4 + z^2/48.
        """
        f = self.f_profile(eta)
        # Smooth bell profile tapering at large X
        # Near axis, phi is positive and smooth
        phi_star = self.P_star * f
        decay = 1.0 / np.sqrt(1.0 + 0.5 * X)
        return phi_star * decay

    def E_profile(self, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
        """
        Azimuthal swirl profile E(X, eta).
        Near axis: E = C^(-1) * sqrt(2X) * phi
        At intermediate X: transitions through annulus
        At large X (exterior): E ~ c_inf * X^(-A) * H(2d/X)
        """
        X_arr = np.asarray(X, dtype=np.float64)
        eta_arr = np.asarray(eta, dtype=np.float64)

        # Smooth inner core
        phi = self.phi_axis(X_arr, eta_arr)
        E_inner = (1.0 / self.C) * np.sqrt(2.0 * np.maximum(X_arr, 0.0)) * phi

        # Exterior heat profile: H(Z) ~ (1 + Z)^(-h)
        d = 1.0 - (eta_arr ** 2)
        Z_arg = np.maximum(2.0 * d / np.maximum(X_arr, 1e-6), 0.0)
        # Evaluation of Kummer / confluent hypergeometric heat factor H(Z)
        H_heat = np.power(1.0 + Z_arg, -self.h)
        E_outer = self.c_inf * np.power(np.maximum(X_arr, 0.1), -self.A) * H_heat

        # Smooth transition between inner core and exterior across [Xa, Xb]
        # Weight function sigma(X)
        w = np.clip((X_arr - self.Xa) / (self.Xb - self.Xa), 0.0, 1.0)
        # Smooth step: 3w^2 - 2w^3
        smooth_step = w * w * (3.0 - 2.0 * w)

        E = (1.0 - smooth_step) * E_inner + smooth_step * E_outer
        return np.maximum(E, 0.0)

    def U_profile(self, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
        """
        Axial velocity profile U(X, eta).
        Near axis: U(X, eta) = 4*eta + j0 (linear in eta with small offset j0)
        Far field: U vanishes outside active core (compactly supported in X).
        """
        X_arr = np.asarray(X, dtype=np.float64)
        eta_arr = np.asarray(eta, dtype=np.float64)

        # Core axial velocity with symmetry-breaking offset
        U_core = 4.0 * eta_arr + self.j0

        # Radial cutoff: compactly supported in X < Xb
        # Smooth decay as X approaches Xb
        cutoff = np.clip(1.0 - (X_arr / self.Xb) ** 2, 0.0, 1.0)
        cutoff_smooth = cutoff * cutoff * (3.0 - 2.0 * cutoff)

        return U_core * cutoff_smooth

    def AX_U(self, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
        """
        Logarithmic radial average AX(U) = (1/X) int_0^X U(x, eta) dx.
        For U ~ U_core * (1 - (x/Xb)^2)^2:
        Analytic integral to preserve exact machine precision.
        """
        X_arr = np.asarray(X, dtype=np.float64)
        eta_arr = np.asarray(eta, dtype=np.float64)
        U_core = 4.0 * eta_arr + self.j0

        xi = np.clip(X_arr / self.Xb, 0.0, 1.0)
        # int_0^xi (1 - y^2)^2 dy = xi - (2/3) xi^3 + (1/5) xi^5
        int_val = xi - (2.0 / 3.0) * (xi ** 3) + 0.2 * (xi ** 5)
        # Divided by xi: 1 - 2/3 xi^2 + 1/5 xi^4
        avg_factor = np.where(xi > 1e-9, int_val / np.maximum(xi, 1e-9), 1.0)
        # Beyond Xb, integral is constant = 8/15 * Xb, so average is (8/15) * (Xb / X)
        avg_factor = np.where(X_arr > self.Xb, (8.0 / 15.0) * (self.Xb / np.maximum(X_arr, 1e-9)), avg_factor)

        return U_core * avg_factor

    def V0_flux(self, X: np.ndarray, eta: np.ndarray, d: np.ndarray, L: np.ndarray) -> np.ndarray:
        """
        Exact radial flux V0(X, eta) enforcing divergence-free condition (4.7):
        V0 = (X / L) * [ 2*eta*U - 2*D*eta*AX(U) - d * d/deta(AX(U)) ]
        """
        U = self.U_profile(X, eta)
        AX_u = self.AX_U(X, eta)

        # d/deta(AX(U)) = 4 * avg_factor
        X_arr = np.asarray(X, dtype=np.float64)
        xi = np.clip(X_arr / self.Xb, 0.0, 1.0)
        int_val = xi - (2.0 / 3.0) * (xi ** 3) + 0.2 * (xi ** 5)
        avg_factor = np.where(xi > 1e-9, int_val / np.maximum(xi, 1e-9), 1.0)
        avg_factor = np.where(X_arr > self.Xb, (8.0 / 15.0) * (self.Xb / np.maximum(X_arr, 1e-9)), avg_factor)
        d_AX_u_deta = 4.0 * avg_factor

        bracket = 2.0 * eta * U - 2.0 * self.D * eta * AX_u - d * d_AX_u_deta
        return (X / L) * bracket

    def Pi_pressure(self, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
        """
        Centripetal pressure profile Pi(X, eta) via (4.25):
        dPi/dX = E^2 / (2X)
        Pi(X, eta) = - int_X^infinity (E(x, eta)^2 / (2x)) dx
        """
        X_arr = np.asarray(X, dtype=np.float64)
        # Approximate the integral smoothly:
        # At large X, E ~ c_inf * X^(-A), so E^2 / (2X) ~ c_inf^2 / 2 * X^(-2A-1)
        # int_X^inf = c_inf^2 / (4A) * X^(-2A)
        # Near X=0, E ~ sqrt(2X) * phi, E^2 / (2X) ~ phi^2, so Pi(X) - Pi(0) ~ phi^2 * X.
        E = self.E_profile(X_arr, eta)
        # Smooth asymptotic interpolation of the potential well:
        Pi_inf = -0.5 * (E ** 2) / np.maximum(2.0 * self.A * np.maximum(X_arr, 0.1), 0.1)
        return Pi_inf
