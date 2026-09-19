import json
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from exposure import dose, validate
from substance_catalog import flammable_thresholds
from combined import CombinedCloud
import input_data

class ExposureTests(unittest.TestCase):
    def test_units_and_linear_history(self):
        p=SimpleNamespace(times=np.array([0.,30.,60.]),
            concentration=lambda t,x,y,z:np.asarray(x)*0+.001*t/60)
        self.assertAlmostEqual(float(dose(p,0.)),.5)

    def test_secondary_pulse_and_combined_common_window(self):
        p=SimpleNamespace(times=np.array([0.,3.,10.]),
            concentration=lambda t,x,y,z:np.asarray(x)*0+.001)
        s=SimpleNamespace(feed=SimpleNamespace(duration_s=4.),options=SimpleNamespace(observation_s=10.),
            x=np.array([0.,100.]),solution=SimpleNamespace(sol=lambda x:x),
            state=lambda x,v:{'travel_time_s':np.asarray(x)*0+2.},
            concentration=lambda x,y,z,t:np.asarray(x)*0+.002)
        self.assertAlmostEqual(float(dose(s,1.)),.002*4*1000/60)
        m=SimpleNamespace(primary=p,secondary=s,time_end=3.)
        self.assertAlmostEqual(float(dose(m,1.)),(.001*3+.002*1)*1000/60)

    def test_limits_recompute_with_temperature(self):
        cold=flammable_thresholds('Пропан',.044,300.,100000.)
        hot=flammable_thresholds('Пропан',.044,600.,100000.)
        np.testing.assert_allclose(cold,np.array(hot)*2)
        self.assertAlmostEqual(cold[1]/cold[0],9.5/2)
        self.assertEqual(flammable_thresholds('Хлор',.0709,300.,100000.),[])

    def test_dose_only_input_and_validation(self):
        path=Path(__file__).resolve().parents[1]/'examples/chloromethane.json'
        d=json.loads(path.read_text());d['toxicity']={'pct50_mg_min_l':.6,'lct50_mg_min_l':6.}
        out=input_data.collect(d,input_data.display_values(d),True,'','','')
        self.assertEqual(out['thresholds_kg_m3'],[])
        with self.assertRaises(ValueError):validate({'toxicity':{'pct50_mg_min_l':-1}})
