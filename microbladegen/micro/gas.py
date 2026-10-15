"""Ideal-gas helpers shared by the micro-turbomachinery modules (air by default)."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Gas:
    """Calorically perfect gas with Sutherland viscosity."""
    gamma: float = 1.4
    R: float = 287.05          # J/(kg K)
    mu_ref: float = 1.716e-5   # Pa s at T_ref (air)
    T_ref: float = 273.15      # K
    S: float = 110.4           # K, Sutherland constant (air)

    @property
    def cp(self) -> float:
        return self.R * self.gamma / (self.gamma - 1.0)

    def viscosity(self, T):
        """Sutherland dynamic viscosity (Pa s)."""
        T = np.asarray(T, dtype=float)
        return self.mu_ref * (T / self.T_ref) ** 1.5 * (self.T_ref + self.S) / (T + self.S)

    def speed_of_sound(self, T):
        return np.sqrt(self.gamma * self.R * np.asarray(T, dtype=float))

    def static_from_total(self, Tt: float, Pt: float, velocity):
        """Static T, p, rho from total conditions and a velocity magnitude (isentropic)."""
        v = np.asarray(velocity, dtype=float)
        T = Tt - v**2 / (2.0 * self.cp)
        p = Pt * (T / Tt) ** (self.gamma / (self.gamma - 1.0))
        rho = p / (self.R * T)
        return T, p, rho


AIR = Gas()
