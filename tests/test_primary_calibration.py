import csv
import json
from pathlib import Path
import unittest

import numpy as np
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud


class PrimaryCalibrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        case = Path(__file__).resolve().parents[1] / 'validation' / 'primary_chlorine'
        data = json.loads((case / 'input.json').read_text(encoding='utf-8'))
        cls.fit = json.loads((case / 'fit_parameters.json').read_text(encoding='utf-8'))
        with (case / 'reference.csv').open() as stream:
            cls.rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
        cls.model = PrimaryCloud(Gas(**data['gas']), Vessel(**data['vessel']),
                                 Atmosphere(**data['atmosphere']), Options(**data['options']))
        cls.result = cls.model.run()

    def test_control_rows_not_used_for_fitting(self):
        self.assertFalse(set(self.fit['train_rows']) & set(self.fit['holdout_rows']))
        checked = 0
        for row in self.rows:
            if int(row['row']) not in self.fit['holdout_rows']:
                continue
            s = self.model.state(self.model.solution.sol(row['time_s']))
            self.assertLess(abs(s['centre_x_m']/row['centre_x_m']-1), .004)
            self.assertLess(abs(np.sqrt(s['sy2_m2'])/row['sy_m']-1), .004)
            self.assertLess(abs(s['centre_concentration_kg_m3']/row['centre_concentration_kg_m3']-1), .01)
            checked += 1
        self.assertEqual(checked, 21)

    def test_full_exposure_and_calibration_provenance(self):
        self.assertEqual(self.model.times[-1], 1800.)
        self.assertEqual(self.result['stop_reason'], 'requested_duration')
        self.assertFalse(self.result['primary_calibration']['independently_validated'])
        self.assertTrue(self.result['primary_calibration']['default_for_primary'])

    def test_dose_boundaries_at_reference_axes(self):
        # Independent time integration at the reported downwind boundaries.
        # It checks the physical field, not values copied into report output.
        from scipy.integrate import trapezoid
        from scipy.optimize import brentq
        t = self.model.times
        states = self.model.states
        cc = np.array([s['centre_concentration_kg_m3'] for s in states])
        xc = np.array([s['centre_x_m'] for s in states])
        rc2 = np.array([s['core_radius_m']**2 for s in states])
        sy2 = np.array([s['sy2_m2'] for s in states])
        def dose(x):
            return trapezoid(cc*np.exp(-np.maximum((x-xc)**2-rc2, 0)/sy2), t)
        # Zone distances were not used to fit trajectory/dispersions.
        for threshold, expected in [(.036, 1396.48), (.36, 356.69)]:
            right = brentq(lambda x: dose(x)-threshold, 0, 3000)
            self.assertLess(abs(right/expected-1), .01)


if __name__ == '__main__':
    unittest.main()
