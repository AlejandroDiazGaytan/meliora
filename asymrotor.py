"""
Root-level compatibility module for Dr. Mario Ulises Delgado-Jaime's `asymrotor.py`
and the extended `AsymRotor_Watson` (Option B: Watson + Coriolis + Laser Mode-Hop).

Allows direct execution of:
    from asymrotor import asym_rotor, asym_rotor_watson, asym_rotor_lines, rotor_levels, rotor_levels_watson, planar_top_c, ray_kappa
"""

from specwatson.asymrotor import (
    asym_rotor,
    asym_rotor_watson,
    asym_rotor_lines,
    rotor_levels,
    rotor_levels_watson,
    planar_top_c,
    ray_kappa,
    convergent_grid_fit,
)

__all__ = [
    "asym_rotor",
    "asym_rotor_watson",
    "asym_rotor_lines",
    "rotor_levels",
    "rotor_levels_watson",
    "planar_top_c",
    "ray_kappa",
    "convergent_grid_fit",
]
