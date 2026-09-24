import json
import tempfile
import unittest
from math import pi
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc
from combined import CombinedCloud, validate_pair
from flammable_mass import MassIntegrator, column_mass, limits, validate, save_flammable_mass

ROOT=Path(__file__).resolve().parents[1]


def gaussian_primary():
    return NS(times=np.array([0.,10.]),beta=2.,solution=NS(sol=lambda t:t),
        state=lambda t:dict(centre_concentration_kg_m3=4.,core_radius_m=0.,
                            sy2_m2=9.,centre_x_m=t,sz_m=2.))


class FlammableMassTests(unittest.TestCase):
    def test_vertical_shell_excludes_rich_core(self):
        for beta in (1.,1.4,2.):
            a=np.array([4.,2.]);s=np.array([.8,3.])
            def c(z):return float(np.sum(a*np.exp(-(z/s)**beta)))
            from scipy.optimize import brentq
            z0=brentq(lambda z:c(z)-3.,0.,100.)
            z1=brentq(lambda z:c(z)-.1,0.,100.)
            expected=quad(c,z0,z1,epsabs=1e-10)[0]
            self.assertAlmostEqual(float(column_mass(a,s,beta,.1,3.)),expected,places=9)
        self.assertEqual(float(column_mass([.01],[1.],2.,.1,3.)),0.)

    def test_overlapping_subthreshold_components_are_counted(self):
        self.assertEqual(float(column_mass([.7],[1.],2.,1.,3.)),0.)
        self.assertGreater(float(column_mass([.7,.7],[1.,1.],2.,1.,3.)),0.)
        self.assertAlmostEqual(float(column_mass([.7,.7],[1.,1.],2.,1.,3.)),
                               float(column_mass([1.4],[1.],2.,1.,3.)),places=10)

    def test_full_3d_gaussian_against_analytic_shell(self):
        m=MassIntegrator(gaussian_primary(),.1,1.)
        expected=4*pi**1.5*9*2/2*(gammainc(1.5,np.log(4/.1))-gammainc(1.5,np.log(4/1.)))
        value,clipped=m.mass(5.,192,128)
        self.assertFalse(clipped)
        self.assertLess(abs(value/expected-1),.0003)

    def test_secondary_finite_pulse_and_spatial_truncation(self):
        def state(x,v):
            x=np.asarray(x);one=np.ones_like(x)
            return dict(travel_time_s=x/2,centre_concentration_kg_m3=one*4,
                        core_half_width_m=one*0,sy_m=one*3,sz_m=one*2)
        s=NS(feed=NS(duration_s=2.),options=NS(observation_s=20.),beta=2.,
             x=np.array([0.,20.]),solution=NS(sol=lambda x:x),state=state)
        m=MassIntegrator(s,.1,1.)
        per_length=pi*3*2/2*(1-.1)
        for t,length in [(0.,0.),(1.,2.),(2.,4.),(5.,4.),(11.,2.),(13.,0.)]:
            value,clipped=m.mass(t,96,128)
            self.assertLess(abs(value-per_length*length),.005*max(1.,per_length*length))
            self.assertEqual(clipped,t>10.)

    def test_real_model_mass_balance(self):
        a=json.loads((ROOT/'examples/chloromethane.json').read_text())
        b=json.loads((ROOT/'examples/secondary_gas.json').read_text())
        a['options']['duration_s']=10.;b['plume_options']['observation_s']=10.
        b['feed']['duration_s']=2.
        m=CombinedCloud.from_inputs(a,b)
        for t in (0.,1.,5.):
            expected=m.primary.q+b['feed']['rate_kg_s']*min(t,2.)
            mass,truncated=MassIntegrator(m,1e-12,1e6).mass(t,192,128)
            self.assertFalse(truncated)
            self.assertLess(abs(mass/expected-1),.008)

    def test_limits_validation_and_conflicting_pair(self):
        self.assertIsNone(limits({'threshold_labels':['НКПР'],'thresholds_kg_m3':[.1]}))
        with self.assertRaises(ValueError):validate({'threshold_labels':['НКПР','ВКПР'],'thresholds_kg_m3':[1.,.1]})
        for cfg in ({'nx':15},{'ny':16.5},{'time_step_s':float('nan')}):
            with self.assertRaises(ValueError):validate({'flammable_mass':cfg})
        a=json.loads((ROOT/'examples/chloromethane.json').read_text())
        b=json.loads((ROOT/'examples/secondary_gas.json').read_text())
        for d in (a,b):d.update(threshold_labels=['НКПР','ВКПР'],thresholds_kg_m3=[.1,.2])
        validate_pair(a,b)
        b['thresholds_kg_m3'][1]=.3
        with self.assertRaises(ValueError):validate_pair(a,b)

    def test_output_and_missing_limits(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=dict(threshold_labels=['НКПР','ВКПР'],thresholds_kg_m3=[.1,1.],
                   flammable_mass=dict(time_step_s=3.,nx=32,ny=32))
            save_flammable_mass(gaussian_primary(),d,tmp)
            result=json.loads((Path(tmp)/'flammable_mass.json').read_text())
            self.assertEqual([r['time_s'] for r in result['history']],[0.,3.,6.,9.,10.])
            self.assertGreater(result['sampled_max_mass_kg'],0.)
            self.assertTrue((Path(tmp)/'flammable_mass.png').is_file())
        with tempfile.TemporaryDirectory() as tmp:
            save_flammable_mass(gaussian_primary(),{},tmp)
            self.assertEqual(json.loads((Path(tmp)/'flammable_mass.json').read_text())['status'],'not_calculated')
