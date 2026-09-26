"""Reproduce single-case fitting; writes candidates, never changes production constants.

python fit_primary_calibration.py
Control rows are omitted from the optimizer, not independent accident validation.
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud

ROOT = Path(__file__).resolve().parent
CASE = ROOT / 'validation' / 'primary_chlorine'


def main():
    data = json.loads((CASE / 'input.json').read_text(encoding='utf-8'))
    with (CASE / 'reference.csv').open() as stream:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    times = np.array([r['time_s'] for r in rows])
    ref_x = np.array([r['centre_x_m'] for r in rows])
    ref_sy = np.array([r['sy_m'] for r in rows])
    usable = np.flatnonzero(times >= .01)
    train, control = usable[::2], usable[1::2]
    model = PrimaryCloud(Gas(**data['gas']), Vessel(**data['vessel']),
                         Atmosphere(**data['atmosphere']), Options(**data['options']))
    # Thermodynamic/radius equations do not depend on x or Sy. Freeze Sy and
    # integrate nominal speed for fitting, without introducing a user-facing mode.
    # The data's final row (1820.93 s) is included only in this fitting calculation;
    # production exposure remains limited to 1800 s and its physical stop events.
    def base_rhs(t, y):
        derivative = model.rhs(t, y)
        derivative[2] = 0.
        derivative[4] = model.state(y)['speed_m_s']
        return derivative
    y0 = [model.q, model.initial['radius_m'], 0.,
          model.q * model.cv_g * model.initial['temperature_k'], 0.]
    baseline = solve_ivp(base_rhs, (0., 1821.), y0, method='DOP853',
                         dense_output=True, max_step=.5, rtol=1e-8, atol=1e-10)
    if not baseline.success:
        raise RuntimeError(baseline.message)
    gain = float(np.exp(np.mean(np.log(ref_x[train]/baseline.sol(times[train])[4]))))

    def predict(parameters):
        a, p, b, q = parameters
        def lateral(t, y):
            state = model.state(baseline.sol(t))
            sy = np.sqrt(max(y[0], 0.))
            core = np.sqrt(max(state['radius_m']**2-y[0], 0.))
            x = gain*state['centre_x_m']
            sigma_prime = model.delta*(1+.5e-4*x)/(1+1e-4*x)**1.5
            early, late = min(t/600., 1.), max(t/600., 1.)
            factor = a*early**p/(1+b*early)*late**q
            return [4*np.sqrt(2/np.pi)*state['speed_m_s'] *
                    (core+.5*np.sqrt(np.pi)*sy)*sigma_prime*factor]
        sol = solve_ivp(lateral, (0., 1821.), [ref_sy[0]**2], rtol=2e-8,
                        atol=1e-12, max_step=4., dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        return np.sqrt(np.maximum(sol.sol(times)[0], 0.))

    fit = least_squares(lambda p: np.log(predict(p)[train]/ref_sy[train]),
                        [1.2, 1., .2, -.1], bounds=([.01, .2, 0., -2.], [4., 2., 5., 2.]),
                        diff_step=1e-4, xtol=1e-8, ftol=1e-8, gtol=1e-8)
    if not fit.success:
        raise RuntimeError(fit.message)
    predicted = predict(fit.x)
    result = dict(drift_gain=gain, **dict(zip(
        ['lateral_gain', 'lateral_power', 'lateral_denominator', 'lateral_late_power'],
        map(float, fit.x))), train_rows=(train+1).tolist(), holdout_rows=(control+1).tolist(),
        control_sy_max_error_percent=float(100*np.max(abs(predicted[control]/ref_sy[control]-1))),
        source='one supplied chlorine accident; not independent validation')
    out = ROOT / 'results_chlorine_validation'
    out.mkdir(exist_ok=True)
    (out / 'fit_candidate.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
