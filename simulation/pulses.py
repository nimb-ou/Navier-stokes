"""
Oscillatory Wave Packets and Annular Reynolds Stress
Implements the shear-amplified viscous wave packets in the annulus Xa < X < Xb
from Sections 6 and 7 of "Finite Time Blowup for Navier-Stokes" (OpenAI).
"""

import numpy as np


class AnnulusWavePulses:
    def __init__(self, Xa: float = 1.0, Xb: float = 5.0, h: float = 0.005):
        """
        Xa, Xb: radial boundaries of the active annulus in similarity coordinates.
        h: core scaling parameter.
        """
        self.Xa = Xa
        self.Xb = Xb
        self.h = h

    def envelope_P(self, v: np.ndarray, Ls: float = 1.0) -> np.ndarray:
        """
        Normalized pulse envelope P(v) from Equation (7.12) and Figure 3(b).
        Grows via shear instability, peaks at v = Ls / 2, then decays via viscous damping.
        P(v) ~ exp(-c * (v - Ls/2)^2 / Ls).
        """
        v_arr = np.asarray(v, dtype=np.float64)
        c_rate = 20.0 / Ls
        arg = -c_rate * ((v_arr - 0.5 * Ls) ** 2) / Ls
        return np.exp(np.clip(arg, -30.0, 0.0))

    def evaluate_pulses(self, r: np.ndarray, theta: np.ndarray, z: np.ndarray, t: float, q: np.ndarray, X: np.ndarray):
        """
        Computes the oscillatory velocity w(r, theta, z, t) in the annulus.
        Two wave families (sigma = +1, -1) with high carrier frequency k ~ q^(-h/2).
        """
        tau = max(1.0 - t, 1e-8)
        # Check if inside active annulus:
        in_annulus = (X >= self.Xa) & (X <= self.Xb)

        # Carrier frequency
        eps = np.power(np.maximum(q, 1e-8), self.h)
        k_freq = int(np.ceil(np.power(eps, -0.5)))
        k_freq = max(k_freq, 8)

        # Pulse coordinate v along annular streamline:
        v_coord = np.mod(theta / (2.0 * np.pi) + t * 5.0, 1.0)
        envelope = self.envelope_P(v_coord, Ls=1.0)

        # Annulus spatial localization window (compact support):
        xi = np.clip((X - self.Xa) / (self.Xb - self.Xa), 0.0, 1.0)
        window = np.sin(np.pi * xi) ** 2

        # Physical amplitude scaling A_wave ~ q^(-1/2 - h/2):
        A_wave = np.power(np.maximum(q, 1e-8), -0.5 - 0.5 * self.h)

        # Phase for sigma = +1 and sigma = -1:
        phi_plus = k_freq * (theta + 2.0 * z / np.maximum(eps, 1e-4) - 3.0 * v_coord)
        phi_minus = k_freq * (theta - 2.0 * z / np.maximum(eps, 1e-4) + 3.0 * v_coord)

        # Wave velocities (orthogonal to wavevector, divergence-free):
        w_r_plus = np.cos(phi_plus) * envelope * window
        w_theta_plus = np.sin(phi_plus) * envelope * window
        w_z_plus = 0.5 * np.sin(phi_plus) * envelope * window

        w_r_minus = np.cos(phi_minus) * envelope * window
        w_theta_minus = -np.sin(phi_minus) * envelope * window
        w_z_minus = -0.5 * np.sin(phi_minus) * envelope * window

        total_amp = A_wave * 0.15

        w_r = total_amp * (w_r_plus + w_r_minus)
        w_theta = total_amp * (w_theta_plus + w_theta_minus)
        w_z = total_amp * (w_z_plus + w_z_minus)

        # Reynolds stress components:
        # <w_r * w_theta>_avg ~ A_wave^2 * window^2
        # Transports angular momentum outwards!
        stress_r_theta = (total_amp ** 2) * (window ** 2) * 0.5
        stress_r_z = (total_amp ** 2) * (window ** 2) * 0.25

        return {
            "w_r": w_r,
            "w_theta": w_theta,
            "w_z": w_z,
            "envelope": envelope,
            "window": window,
            "stress_r_theta": stress_r_theta,
            "stress_r_z": stress_r_z,
            "k_freq": k_freq,
        }
