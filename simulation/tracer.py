"""
Particle Tracer & Streamline Integrator
Simulates 3D parcel trajectories dx/dt = u(x, t) in the Navier-Stokes blowup field
using high-precision 4th-order Runge-Kutta (RK4) integration.
"""

import numpy as np
from .fields import NavierStokesBlowupField


class ParticleTracer:
    def __init__(self, field: NavierStokesBlowupField):
        self.field = field

    def velocity_fn(self, pos: np.ndarray, t: float) -> np.ndarray:
        """
        pos: shape (N, 3) or (3,) containing (x, y, z)
        returns: velocity (vx, vy, vz) with same shape
        """
        pos_arr = np.asarray(pos, dtype=np.float64)
        single = (pos_arr.ndim == 1)
        if single:
            pos_arr = pos_arr.reshape(1, 3)

        x = pos_arr[:, 0]
        y = pos_arr[:, 1]
        z = pos_arr[:, 2]

        eval_res = self.field.evaluate_cartesian(x, y, z, t)
        vel = np.column_stack([eval_res["u_x"], eval_res["u_y"], eval_res["u_z"]])

        if single:
            return vel[0]
        return vel

    def integrate_trajectories(self, initial_positions: np.ndarray, t_start: float = 0.8,
                               t_end: float = 0.999, dt: float = 0.001) -> dict:
        """
        Integrates a swarm of particles forward in time using RK4:
          initial_positions: (N, 3) array of (x, y, z)
          t_start, t_end: time window in [0, 1)
          dt: time step
        """
        positions = np.array(initial_positions, dtype=np.float64)
        if positions.ndim == 1:
            positions = positions.reshape(1, 3)

        n_particles = positions.shape[0]
        t = t_start

        times = [t]
        trajectories = [positions.copy()]
        speeds = [self.field.evaluate_cartesian(positions[:, 0], positions[:, 1], positions[:, 2], t)["speed"]]

        while t < t_end:
            step = min(dt, t_end - t)
            # Adaptive step near singularity as velocities grow:
            v_curr = self.velocity_fn(positions, t)
            max_speed = np.max(np.linalg.norm(v_curr, axis=1))
            if max_speed > 50.0:
                # Subdivide step for numerical stability
                step = min(step, 0.2 / max_speed)

            # RK4 stages:
            k1 = self.velocity_fn(positions, t)
            k2 = self.velocity_fn(positions + 0.5 * step * k1, t + 0.5 * step)
            k3 = self.velocity_fn(positions + 0.5 * step * k2, t + 0.5 * step)
            k4 = self.velocity_fn(positions + step * k3, t + step)

            positions = positions + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            t += step

            times.append(t)
            trajectories.append(positions.copy())
            speeds.append(self.field.evaluate_cartesian(positions[:, 0], positions[:, 1], positions[:, 2], t)["speed"])

            if t >= t_end or (1.0 - t) < 1e-6:
                break

        # Convert to arrays: shape (n_steps, n_particles, 3)
        return {
            "times": np.array(times),
            "trajectories": np.array(trajectories),
            "speeds": np.array(speeds),
        }
