import json
from dataclasses import replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import unittest
from math import isclose
from scipy.integrate import quad
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud, R_GAS

DATA = json.loads((Path(__file__).resolve().parents[1] / 'examples/chloromethane.json').read_text(encoding='utf-8'))

def model(duration=30, **options):
    return PrimaryCloud(Gas(**DATA['gas']), Vessel(**DATA['vessel']), Atmosphere(**DATA['atmosphere']), Options(duration_s=duration, **options))


def test_initial_mass_and_geometry():
    m = model()
    # Independent ideal-gas calculation: pV = nRT, n * molecular mass.
    expected = 101325 * 2000 / (8.3144 * 291.15) * .051
    assert isclose(m.initial['mass_kg'], expected, rel_tol=1e-12)
    assert isclose(m.initial['radius_m'], 8.6025401383, rel_tol=1e-10)
    assert isclose(np.pi * m.initial['radius_m']**2 * m.initial['height_m'], 2000, rel_tol=1e-10)
    # Preserve the published discrepancy; do not fit the implementation to 4227.81.
    assert abs(m.initial['mass_kg'] / 4227.81 - 1) < .011
    assert abs(m.initial['mass_kg'] / 4227.81 - 1) > .009


def test_concentration_integral_conserves_released_mass():
    m = model()
    m.run()
    s = m.state(m.solution.sol(30))
    r = s['core_radius_m']
    upper = r + 12 * np.sqrt(s['sy2_m2'])
    # Numerical integration of the actual concentration profile, not a copy of (106).
    radial = quad(lambda q: float(m.concentration(30, s['centre_x_m'] + q)) * 2*np.pi*q,
                  0, upper, points=[r], epsabs=1e-7)[0]
    vertical = quad(lambda z: np.exp(-(z/s['sz_m'])**m.beta), 0, np.inf)[0]
    assert isclose(radial * vertical, m.q, rel_tol=1e-7)


def test_isothermal_invariant_and_monotone_air_mass():
    m = model(120)
    result = m.run()
    assert result['stop_reason'] == 'requested_duration'
    assert np.max(np.abs(np.array([s['temperature_k'] for s in m.states]) - 291.15)) < 1e-8
    assert np.all(np.diff([s['total_mass_kg'] for s in m.states]) > 0)
    assert np.all(np.diff([s['centre_x_m'] for s in m.states]) > 0)


def test_solver_convergence():
    a, b = model(120, max_step_s=1.), model(120, max_step_s=.25, rtol=1e-9, atol=1e-11)
    a.run(); b.run()
    np.testing.assert_allclose(a.solution.sol(120), b.solution.sol(120), rtol=2e-6)
    grid = np.array([30., 100., 200., 300.])
    np.testing.assert_allclose(a.maximum_concentration(grid), b.maximum_concentration(grid), rtol=2e-6)


def test_pressure_expansion_and_phase_guard():
    base = model()
    vessel = replace(base.vessel, pressure_abs_pa=2*base.met.pressure_pa, temperature_k=400.)
    m = PrimaryCloud(base.gas, vessel, base.met)
    assert isclose(m.initial['temperature_k'], 400. * 2**(-.2), rel_tol=1e-10)
    assert isclose(m.initial['volume_m3'], 2000. * 2**.8, rel_tol=1e-10)
    m.run()
    assert m.states[-1]['temperature_k'] < m.initial['temperature_k']
    with unittest.TestCase().assertRaisesRegex(ValueError, 'фазовая модель'):
        PrimaryCloud(base.gas, replace(base.vessel, pressure_abs_pa=10*base.met.pressure_pa), base.met)


def test_bad_meteorology_rejected():
    m = model()
    for change in [dict(wind_speed_m_s=0.), dict(wind_speed_m_s=float('nan')), dict(stability='Z'), dict(alpha_basis='')]:
        with unittest.TestCase().assertRaises(ValueError):
            PrimaryCloud(m.gas, m.vessel, replace(m.met, **change))


def test_limits_are_explicit():
    m = model()
    with unittest.TestCase().assertRaisesRegex(ValueError, '600'):
        PrimaryCloud(m.gas, m.vessel, m.met, Options(duration_s=601))
    with unittest.TestCase().assertRaisesRegex(ValueError, 'подтвердить'):
        PrimaryCloud(replace(m.gas, phase_and_ideal_gas_assumptions_confirmed=False), m.vessel, m.met)
    limited = PrimaryCloud(m.gas, m.vessel, replace(m.met, alpha_height_limit_m=9.), Options(duration_s=600))
    result = limited.run()
    assert result['stop_reason'] == 'alpha_height_exceeded'
    assert limited.times[-1] < 600
    with unittest.TestCase().assertRaises(ValueError): limited.concentration(limited.times[-1]+1, 100)


def load_tests(loader, tests, pattern):
    return unittest.TestSuite(unittest.FunctionTestCase(func) for name, func in globals().items() if name.startswith('test_'))

if __name__ == '__main__':
    unittest.main()
