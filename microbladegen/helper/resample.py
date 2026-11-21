
from typing import Tuple
from scipy.interpolate import splprep, splev, PchipInterpolator
import numpy as np 
import numpy.typing as npt




def clean_curve_points(points):
    """Remove consecutive zero-length segments without sorting or opening loops."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] not in (2, 3) or len(points) < 2:
        raise ValueError("Expected at least two 2D or 3D points")
    if not np.isfinite(points).all():
        raise ValueError("Curve coordinates must be finite")
    keep = np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1) > 1e-12]
    points = points[keep]
    if len(points) < 2:
        raise ValueError("Curve must have nonzero length")
    return points


def validate_count(count):
    if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 2:
        raise ValueError("Sample count must be an integer >= 2")


def resample_curve(curve: npt.NDArray, M: int = 100):
    validate_count(M)
    curve = clean_curve_points(curve)
    tck, _ = splprep(curve.T, s=0, per=0, k=min(3, len(curve) - 1))
    return np.stack(splev(np.linspace(0, 1, M), tck), axis=1)


def arc_length_param_strict(points: npt.NDArray) -> Tuple[npt.NDArray, npt.NDArray]:
    points = clean_curve_points(points)
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    arc /= arc[-1]
    return arc, points

def curvature_2d(x, y):
    dx = np.gradient(x)
    dy = np.gradient(y)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    numerator = np.abs(dx * ddy - dy * ddx)
    denominator = (dx**2 + dy**2) ** 1.5 + 1e-8
    return numerator / denominator

def curvature_3d(x, y, z):
    dx = np.gradient(x)
    dy = np.gradient(y)
    dz = np.gradient(z)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    ddz = np.gradient(dz)
    cross = np.cross(np.stack([dx, dy, dz], axis=1),
                     np.stack([ddx, ddy, ddz], axis=1))
    num = np.linalg.norm(cross, axis=1)
    denom = (dx**2 + dy**2 + dz**2) ** 1.5 + 1e-8
    return num / denom

def resample_by_curvature(points: npt.NDArray, N: int) -> npt.NDArray:
    """
    Resample 2D or 3D curve using curvature-based importance sampling and PCHIP interpolation.
    """
    validate_count(N)
    points = clean_curve_points(points)

    dim = points.shape[1]

    # Interpolate each dimension with PCHIP
    arc, points = arc_length_param_strict(points)
    pchips = [PchipInterpolator(arc, points[:, i]) for i in range(dim)]
    
    u_fine = np.linspace(0, 1, 1000)
    coords_fine = np.stack([p(u_fine) for p in pchips], axis=1)

    # Compute curvature
    if dim == 2:
        kappa = curvature_2d(coords_fine[:, 0], coords_fine[:, 1])
    else:
        kappa = curvature_3d(coords_fine[:, 0], coords_fine[:, 1], coords_fine[:, 2])

    # Normalize curvature into a sampling density
    density = kappa + 1e-3
    density /= np.sum(density)
    cdf = np.cumsum(density)
    cdf /= cdf[-1]

    # Importance sample along CDF
    u_resampled = np.interp(np.linspace(0, 1, N), cdf, u_fine)
    resampled = np.stack([p(u_resampled) for p in pchips], axis=1)
    return resampled