"""Geometric postprocessing of (104)-(106); no new dispersion coefficients.
A horizontal section at any instant is a disk. Its swept envelope is the
zone where the specified concentration was reached during the simulated time.
All maxima/envelopes are sampled in time; clipping and time bounds are explicit.
"""
from math import pi
import numpy as np
from scipy.optimize import brentq


def validate_threshold(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not np.isfinite(value) or value <= 0:
        raise ValueError('Порог концентрации должен быть конечным положительным числом, кг/м³')


def section(model, t, threshold, height_m=0.):
    """Disk with C >= threshold. None = empty section, not a disk of radius zero."""
    validate_threshold(threshold)
    if not np.isfinite(height_m) or height_m < 0:
        raise ValueError('Высота сечения должна быть конечной и неотрицательной')
    if t < 0 or t > model.times[-1]:
        raise ValueError('Время вне рассчитанного интервала')
    s = model.state(model.solution.sol(t))
    # Logarithmic form avoids underflow at large section heights.
    log_ratio = np.log(s['centre_concentration_kg_m3'] / threshold) - (height_m / s['sz_m'])**model.beta
    if log_ratio < 0:
        return None
    radius2 = s['core_radius_m']**2 + s['sy2_m2'] * log_ratio
    radius = float(np.sqrt(max(radius2, 0.)))
    return dict(time_s=float(t), centre_x_m=float(s['centre_x_m']), radius_m=radius,
                x_min_m=float(s['centre_x_m']-radius), x_max_m=float(s['centre_x_m']+radius),
                area_m2=pi*radius2)


def envelope(model, threshold, height_m=0., spatial_points=2001):
    """Union of the sampled disks in a horizontal plane.
    Area integrates their upper/lower envelope; NOT pi*R_eff**2.
    An envelope is returned only for the computed interval, never extrapolated.
    """
    validate_threshold(threshold)
    if isinstance(spatial_points, bool) or not isinstance(spatial_points, int) or not 101 <= spatial_points <= 100001:
        raise ValueError('spatial_points: целое число от 101 до 100001')
    disks = [disk for t in model.times if (disk := section(model, t, threshold, height_m)) is not None]
    base = dict(threshold_kg_m3=float(threshold), section_height_m=float(height_m),
                time_start_s=0., time_end_s=float(model.times[-1]),
                temporal_method='union_of_sampled_disks',
                active_at_end=section(model, model.times[-1], threshold, height_m) is not None)
    if not disks:
        return dict(base, status='not_reached_in_computed_interval', x_min_m=None,
                    x_max_m=None, max_width_m=0., footprint_area_m2=0.,
                    x_m=[], half_width_m=[], disks=[])
    xmin = min(d['x_min_m'] for d in disks)
    xmax = max(d['x_max_m'] for d in disks)
    x = np.linspace(xmin, xmax, spatial_points)
    half2 = np.zeros_like(x)
    for d in disks:
        half2 = np.maximum(half2, d['radius_m']**2 - (x-d['centre_x_m'])**2)
    half = np.sqrt(half2)
    # np.trapz was removed in newer NumPy. np.trapezoid appeared in NumPy 2.
    integrate = getattr(np, 'trapezoid', None)
    if integrate is None: integrate = np.trapz
    return dict(base, status='reached_in_computed_interval', x_min_m=xmin, x_max_m=xmax,
                max_width_m=2*max(d['radius_m'] for d in disks),
                footprint_area_m2=float(integrate(2*half, x)),
                x_m=x.tolist(), half_width_m=half.tolist(), disks=disks)


def receptor(model, x_m, y_m=0., z_m=0., thresholds=()):
    """Concentration history and threshold exposure in a fixed point.
    Bracketed threshold crossings are refined on the dense ODE solution.
    Extremely brief crossings between all sample points can still be missed.
    """
    for value in (x_m, y_m, z_m):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
            raise ValueError('Координаты точки должны быть конечными числами')
    if z_m < 0: raise ValueError('z >= 0')
    c = np.array([float(model.concentration(t, x_m, y_m, z_m)) for t in model.times])
    peak_index = int(c.argmax())
    bands = []
    for threshold in thresholds:
        validate_threshold(threshold)
        above = c >= threshold
        intervals = []
        start = 0. if above[0] else None
        f = lambda t: float(model.concentration(t, x_m, y_m, z_m)) - threshold
        for i in range(len(c)-1):
            if above[i] == above[i+1]: continue
            # The calibrated early cloud has a steep edge; refine the crossing
            # enough to avoid a visible concentration residual at entry/exit.
            crossing = float(brentq(f, model.times[i], model.times[i+1], xtol=1e-12))
            if above[i+1]:
                start = crossing
            else:
                intervals.append([float(start), crossing]); start = None
        if start is not None: intervals.append([float(start), float(model.times[-1])])
        bands.append(dict(threshold_kg_m3=threshold, intervals_s=intervals,
                          arrival_s=intervals[0][0] if intervals else None,
                          last_exit_s=intervals[-1][1] if intervals and not above[-1] else None,
                          duration_in_computed_interval_s=sum(b-a for a,b in intervals),
                          active_at_end=bool(above[-1])))
    return dict(x_m=x_m, y_m=y_m, z_m=z_m, time_s=model.times.tolist(),
                concentration_kg_m3=c.tolist(), max_concentration_kg_m3=float(c[peak_index]),
                time_of_sampled_max_s=float(model.times[peak_index]), threshold_exposure=bands,
                method='sampled_history_with_bracketed_crossings',
                time_end_s=float(model.times[-1]))
