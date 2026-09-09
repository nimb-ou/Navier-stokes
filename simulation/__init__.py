"""
Navier-Stokes Finite Time Blowup Simulation Package
Based on OpenAI's formalization of "Finite Time Blowup for Navier-Stokes".
"""

from .similarity import SimilarityCoordinates
from .profiles import BlowupProfiles
from .fields import NavierStokesBlowupField
from .pulses import AnnulusWavePulses
from .tracer import ParticleTracer
from .diagnostics import BlowupDiagnostics

__all__ = [
    "SimilarityCoordinates",
    "BlowupProfiles",
    "NavierStokesBlowupField",
    "AnnulusWavePulses",
    "ParticleTracer",
    "BlowupDiagnostics",
]
