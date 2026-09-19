import copy
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
from combined import CombinedCloud, validate_pair

ROOT = Path(__file__).resolve().parents[1]

class CombinedTests(unittest.TestCase):
    def test_sum_at_same_time_not_separate_peaks(self):
        p = SimpleNamespace(solution=SimpleNamespace(t=[0., 2.]), times=np.array([0., 1., 2.]),
                            concentration=lambda t, x, y, z: np.asarray(x)*0 + (10. if t == 0 else 0.))
        s = SimpleNamespace(options=SimpleNamespace(observation_s=2.),
                            concentration=lambda x, y, z, t: np.asarray(x)*0 + (8. if t == 2 else 0.))
        m = CombinedCloud(p, s)
        np.testing.assert_allclose(m.maximum([0., 1.]), [10., 10.])
        self.assertEqual(float(m.concentration(2., 0.)[2]), 8.)
        with self.assertRaises(ValueError): m.concentration(3., 0.)

    def test_real_models_sum_domain_and_feed_events(self):
        a = json.loads((ROOT/'examples/chloromethane.json').read_text())
        b = json.loads((ROOT/'examples/secondary_gas.json').read_text())
        a['options']['duration_s'] = 10.
        b['plume_options']['observation_s'] = 8.
        b['feed']['duration_s'] = 2.
        m = CombinedCloud.from_inputs(a, b)
        self.assertEqual(m.time_end, 8.)
        for t in [0., 1., 4., 8.]:
            c1,c2,total=m.concentration(t, np.array([-1., 0., 5.]))
            np.testing.assert_allclose(total, c1+c2)
        self.assertTrue(np.isnan(m.concentration(1., m.secondary.x[-1]+1)[2]))
        self.assertIn(2., m.receptor_times(0.))
        before = np.nextafter(2., -np.inf)
        self.assertGreater(float(m.concentration(before, 0.)[1]), 0.)
        self.assertEqual(float(m.concentration(2., 0.)[1]), 0.)
        altered=copy.deepcopy(b);altered['atmosphere']['wind_speed_m_s']+=1
        with self.assertRaises(ValueError):validate_pair(a,altered)
        altered=copy.deepcopy(b);altered['gas']['name']='Другое вещество'
        with self.assertRaises(ValueError):validate_pair(a,altered)
