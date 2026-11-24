"""ABI L1b radiance to brightness temperature.

The Planck coefficients come from the same file as the radiance. The
conversion does not remove reflected sunlight.
"""

from __future__ import annotations

import numpy as np

def radiance_to_bt(rad, fk1, fk2, bc1, bc2):
    """ABI L1b radiance to brightness temperature (K).

    T = [fk2 / ln(fk1 / L + 1) - bc1] / bc2

    The coefficients live on the radiance variable and already match that
    file's radiance units. This is a Planck inversion of the observed
    radiance. It does not remove reflected sunlight, which band 7
    (3.9 µm) still sees while the sun is up.
    """
    rad = np.asarray(rad, dtype=np.float64)
    bt = np.full(rad.shape, np.nan, dtype=np.float64)
    ok = np.isfinite(rad) & (rad > 0)
    bt[ok] = (fk2 / np.log(fk1 / rad[ok] + 1.0) - bc1) / bc2
    return bt

def bt_ceiling(rad, fk1, fk2, bc1, bc2):
    """Brightness temperature of the largest scaled radiance, or None."""
    enc = rad.encoding
    valid = enc.get("valid_range", rad.attrs.get("valid_range"))
    if valid is None:
        return None
    scale = float(enc.get("scale_factor", rad.attrs.get("scale_factor", 1.0)))
    offset = float(enc.get("add_offset", rad.attrs.get("add_offset", 0.0)))
    rad_max = float(valid[1]) * scale + offset
    if rad_max <= 0:
        return None
    return float(radiance_to_bt(np.array([rad_max]), fk1, fk2, bc1, bc2)[0])
