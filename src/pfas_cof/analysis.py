"""Estimators with explicit units; no inference of rejection from an energy barrier."""
import numpy as np

R = 0.008314462618  # kJ mol^-1 K^-1


def msd(positions, dt_ps, max_lag, axes=(0, 1, 2)):
    """Time-origin averaged MSD of unwrapped COMs, shape (frames, particles, 3).

    Remove system drift and align to the membrane before calling. A global z MSD
    in a confining membrane is generally not an asymptotic diffusion estimator.
    """
    x = np.asarray(positions, dtype=float)
    if x.ndim != 3 or x.shape[2] != 3 or not np.isfinite(x).all():
        raise ValueError("Expected finite (frames, particles, 3) positions")
    if x.shape[1] < 1 or dt_ps <= 0 or not 1 <= max_lag < len(x):
        raise ValueError("Invalid particle count, timestep or lag")
    if not axes or len(set(axes)) != len(axes) or any(a not in (0, 1, 2) for a in axes):
        raise ValueError("axes must be distinct Cartesian indices")
    x = x[:, :, axes]
    values = [0.0] + [float(np.mean(np.sum((x[k:] - x[:-k])**2, axis=-1)))
                      for k in range(1, max_lag + 1)]
    return np.arange(max_lag + 1) * dt_ps, np.asarray(values)


def diffusion(time_ps, msd_nm2, fit_min_ps, fit_max_ps, dimensions):
    t, y = np.asarray(time_ps), np.asarray(msd_nm2)
    if t.shape != y.shape or dimensions not in (1, 2, 3) or fit_min_ps >= fit_max_ps:
        raise ValueError("Invalid fit definition")
    if not np.isfinite(t).all() or not np.isfinite(y).all() or np.any(np.diff(t) <= 0):
        raise ValueError("Time must increase and observations must be finite")
    mask = (t >= fit_min_ps) & (t <= fit_max_ps)
    if mask.sum() < 3:
        raise ValueError("At least three fit observations required")
    slope, intercept = np.polyfit(t[mask], y[mask], 1)
    if slope <= 0:
        raise ValueError("Nonpositive slope: no diffusive fit")
    residual = np.sum((y[mask] - (slope*t[mask] + intercept))**2)
    total = np.sum((y[mask] - y[mask].mean())**2)
    return {"D_nm2_ps": float(slope / (2 * dimensions)),
            "D_m2_s": float(slope / (2 * dimensions) * 1e-6),
            "intercept_nm2": float(intercept), "r_squared": float(1-residual/total),
            "fit_min_ps": fit_min_ps, "fit_max_ps": fit_max_ps,
            "warning": "Correlated lag points: use independent replica estimates for uncertainty."}


def wham(samples, centers_nm, force_constants, edges_nm, temperature_K=298.15,
         tolerance=1e-8, max_iterations=20000):
    """Binned 1D harmonic WHAM. k in kJ/mol/nm²; equilibrated samples in nm.

    Assumes uniform bin width and independent/effectively subsampled observations.
    Empty bins retain NaN; disconnected adjacent histograms raise an error.
    This reference implementation does not estimate statistical uncertainty.
    """
    edges = np.asarray(edges_nm, dtype=float)
    c, k = np.asarray(centers_nm), np.asarray(force_constants)
    if len(samples) < 2 or c.shape != (len(samples),) or k.shape != c.shape:
        raise ValueError("Provide matching samples, centers and force constants")
    if temperature_K <= 0 or not np.isfinite(c).all() or not np.isfinite(k).all() or np.any(k <= 0):
        raise ValueError("Invalid restraints or temperature")
    if len(edges) < 3 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0) or not np.allclose(np.diff(edges), np.diff(edges)[0]):
        raise ValueError("Uniform increasing finite edges required")
    arrays = [np.asarray(s, dtype=float) for s in samples]
    if any(s.ndim != 1 or not len(s) or not np.isfinite(s).all() or np.any(s < edges[0]) or np.any(s > edges[-1]) for s in arrays):
        raise ValueError("Samples must be finite, nonempty and within histogram range")
    h = np.array([np.histogram(s, edges)[0] for s in arrays])
    normalized = h / h.sum(axis=1)[:, None]
    order = np.argsort(c)
    overlap = [float(np.minimum(normalized[a], normalized[b]).sum())
               for a, b in zip(order[:-1], order[1:])]
    if min(overlap) <= 0:
        raise ValueError("Disconnected windows: acquire bridging samples")
    z = (edges[1:] + edges[:-1]) / 2
    bias = k[:, None] * (z[None, :] - c[:, None])**2 / (2*R*temperature_K)
    n = h.sum(axis=1)
    occupied = h.sum(axis=0) > 0
    logh = np.log(h.sum(axis=0)[occupied])
    f = np.zeros(len(samples))
    def lse(a, axis):
        m = np.max(a, axis=axis, keepdims=True)
        return np.squeeze(m + np.log(np.exp(a-m).sum(axis=axis, keepdims=True)), axis=axis)
    for iteration in range(max_iterations):
        logp = logh - lse(np.log(n)[:, None] + f[:, None] - bias[:, occupied], 0)
        logp -= lse(logp, 0)
        new = -lse(logp[None, :] - bias[:, occupied], 1)
        new -= new[0]
        if np.max(np.abs(new-f)) < tolerance:
            break
        f = new
    else:
        raise RuntimeError("WHAM did not converge")
    pmf = np.full(len(z), np.nan)
    pmf[occupied] = -R*temperature_K*logp
    pmf[occupied] -= np.min(pmf[occupied])
    return z, pmf, {"iterations": iteration+1, "adjacent_overlap": overlap,
                    "minimum_overlap": min(overlap), "unoccupied_bins": int((~occupied).sum())}


def residence_events(inside, dt_ps):
    """Continuous residence episodes with explicit boundary censoring.

    Boolean (frames, particles). Duration = observed occupied frames * dt.
    First-frame episodes are left-censored, last-frame episodes right-censored.
    Fit survival models to events; averaging completed events alone is biased.
    """
    x = np.asarray(inside)
    if x.ndim != 2 or x.dtype != np.bool_ or dt_ps <= 0:
        raise ValueError("Expected boolean (frames, particles) and positive dt")
    rows = []
    for j in range(x.shape[1]):
        padded = np.r_[False, x[:, j], False].astype(int)
        starts, ends = np.where(np.diff(padded)==1)[0], np.where(np.diff(padded)==-1)[0]
        rows.extend({"particle": j, "duration_ps": float((b-a)*dt_ps),
                     "left_censored": bool(a==0), "right_censored": bool(b==len(x))}
                    for a, b in zip(starts, ends))
    return rows


def water_permeance(n_forward, n_reverse, time_ps, area_nm2, pressure_bar,
                     molecular_volume_nm3=0.0299):
    """Steady net crossing volume / area / time / pressure; counts full passages.

    Pressure must be hydraulic minus osmotic pressure. Count with two reservoir
    boundaries and molecule identities, not every crossing of a single plane.
    """
    if min(time_ps, area_nm2, pressure_bar, molecular_volume_nm3) <= 0 or min(n_forward, n_reverse) < 0:
        raise ValueError("Invalid counts or physical inputs")
    flux_m_s = (n_forward-n_reverse)*molecular_volume_nm3/(time_ps*area_nm2)*1e3
    return {"flux_m_s": flux_m_s, "permeance_m_Pa_s": flux_m_s/(pressure_bar*1e5),
            "permeance_L_m2_h_bar": flux_m_s*3.6e6/pressure_bar}


def hydration_counts(pfas_head_nm, water_oxygen_nm, box_nm, cutoff_nm=0.35):
    """Per-head oxygen coordination for one orthorhombic frame using minimum image."""
    head, water, box = map(np.asarray, (pfas_head_nm, water_oxygen_nm, box_nm))
    if head.ndim != 2 or water.ndim != 2 or head.shape[1] != 3 or water.shape[1] != 3 or box.shape != (3,) or not np.isfinite(box).all() or min(box) <= 2*cutoff_nm or cutoff_nm <= 0:
        raise ValueError("Invalid geometry/cutoff; only orthorhombic boxes supported")
    if not np.isfinite(head).all() or not np.isfinite(water).all():
        raise ValueError("Nonfinite coordinates")
    d = head[:, None, :] - water[None, :, :]
    d -= box * np.rint(d/box)
    return (np.sum(d*d, axis=-1) < cutoff_nm**2).sum(axis=1)
