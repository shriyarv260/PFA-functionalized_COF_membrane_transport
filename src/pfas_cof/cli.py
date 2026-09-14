import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .analysis import msd, diffusion, wham, residence_events, water_permeance


def demo(out):
    """Synthetic Brownian motion and analytically biased harmonic umbrellas."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20250901)
    expected_D = 0.002
    x = np.cumsum(rng.normal(0, np.sqrt(2*expected_D), (2000, 128, 3)), axis=0)
    np.savez_compressed(out/'synthetic_trajectory.npz', positions_nm=x, dt_ps=1.0)
    t, m = msd(x, 1, 100)
    fit = diffusion(t, m, 10, 80, 3)
    np.savetxt(out/'msd.csv', np.c_[t, m], delimiter=',', header='lag_ps,msd_nm2', comments='')
    centers = np.linspace(-1.5, 1.5, 16)
    spring, underlying, rt = 30., 4., 0.008314462618*298.15
    samples = [rng.normal(spring*c/(spring+underlying), np.sqrt(rt/(spring+underlying)), 6000) for c in centers]
    z, pmf, diagnostics = wham(samples, centers, np.full(16, spring), np.linspace(-2.5, 2.5, 101))
    analytic = 0.5*underlying*z*z
    analytic -= analytic.min()
    valid = np.isfinite(pmf) & (np.abs(z) < 1.3)
    np.savetxt(out/'pmf.csv', np.c_[z, pmf, analytic], delimiter=',', header='z_nm,synthetic_pmf_kJ_mol,analytic_kJ_mol', comments='')
    result = {"data_status": "SYNTHETIC VALIDATION ONLY — no PFAS/COF scientific results",
              "seed": 20250901, "diffusion": fit, "expected_D_nm2_ps": expected_D,
              "harmonic_pmf_rmse_kJ_mol": float(np.sqrt(np.mean((pmf[valid]-analytic[valid])**2))),
              "wham": diagnostics,
              "residence_example": residence_events(np.array([[1,0],[1,1],[0,1],[1,0]], dtype=bool), 1),
              "water_example": water_permeance(120,20,1000,25,100)}
    (out/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    (out/'README.md').write_text('# Synthetic verification\n\nNo atomistic PFAS simulations were run. These seeded fixtures validate estimators only.\n\n'+f'Brownian D: {fit["D_nm2_ps"]:.6f} nm²/ps (expected {expected_D}). Harmonic PMF RMSE: {result["harmonic_pmf_rmse_kJ_mol"]:.3f} kJ/mol.\n')
    print(json.dumps(result, indent=2))


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='command', required=True)
    d = sub.add_parser('demo'); d.add_argument('--out', default='results/synthetic')
    m = sub.add_parser('msd')
    m.add_argument('trajectory', help='NPZ: positions_nm (frames,molecules,3), dt_ps')
    m.add_argument('--axes', choices=['xy','xyz','z'], default='xy')
    m.add_argument('--max-lag', type=int, required=True)
    m.add_argument('--fit', nargs=2, type=float, required=True, metavar=('START_PS','END_PS'))
    m.add_argument('--out', required=True)
    w = sub.add_parser('wham')
    w.add_argument('manifest', help='CSV columns: path,center_nm,k_kJ_mol_nm2; paths relative to CSV')
    w.add_argument('--range', nargs=2, type=float, required=True)
    w.add_argument('--bins', type=int, default=100)
    w.add_argument('--temperature', type=float, default=298.15)
    w.add_argument('--out', required=True)
    args = p.parse_args()
    if args.command == 'demo':
        demo(args.out)
    elif args.command == 'msd':
        data = np.load(args.trajectory)
        axes = {'xy':(0,1),'xyz':(0,1,2),'z':(2,)}[args.axes]
        t, y = msd(data['positions_nm'], float(data['dt_ps']), args.max_lag, axes)
        result = diffusion(t, y, *args.fit, len(axes))
        Path(args.out).mkdir(parents=True, exist_ok=True)
        np.savetxt(Path(args.out)/'msd.csv', np.c_[t,y], delimiter=',', header='lag_ps,msd_nm2', comments='')
        (Path(args.out)/'diffusion.json').write_text(json.dumps(result, indent=2))
    else:
        manifest = Path(args.manifest)
        with manifest.open() as f: rows = list(csv.DictReader(f))
        z, pmf, diag = wham([np.loadtxt(manifest.parent/r['path'], ndmin=1) for r in rows],
                           [float(r['center_nm']) for r in rows], [float(r['k_kJ_mol_nm2']) for r in rows],
                           np.linspace(*args.range, args.bins+1), args.temperature)
        Path(args.out).mkdir(parents=True, exist_ok=True)
        np.savetxt(Path(args.out)/'pmf.csv', np.c_[z,pmf], delimiter=',', header='z_nm,pmf_kJ_mol', comments='')
        (Path(args.out)/'diagnostics.json').write_text(json.dumps(diag, indent=2))


if __name__ == '__main__': main()
