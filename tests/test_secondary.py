import json,unittest
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from secondary import SecondaryCloud,GasFeed,PlumeOptions
from cloud import Gas,Atmosphere
import secondary_input
D=json.loads((Path(__file__).resolve().parents[1]/'examples/secondary_gas.json').read_text(encoding='utf-8'))

def make(**changes):
    opt=dict(D['plume_options']);opt.update(changes)
    m=SecondaryCloud(Gas(**D['gas']),GasFeed(**D['feed']),Atmosphere(**D['atmosphere']),PlumeOptions(**opt));m.run();return m

class SecondaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=make()

    def test_initial_cross_section_matches_flow(self):
        m=self.m;b=m.initial_b
        self.assertAlmostEqual(2*b*b*m.rho_gas*m.speed(b),m.q,places=10)

    def test_profile_flux_integral(self):
        m=self.m;x=30.;s=m.state(x,m.solution.sol(x));b=float(s['core_half_width_m']);sy=float(s['sy_m'])
        lateral=2*quad(lambda y:float(m.concentration(x,y,0.)),0,b+12*sy,points=[b])[0]
        vertical=quad(lambda z:np.exp(-(z/float(s['sz_m']))**m.beta),0,np.inf)[0]
        self.assertAlmostEqual(lateral*vertical*float(s['speed_m_s']),m.q,places=7)

    def test_source_shutoff_and_travel(self):
        m=self.m;x=50.;tau=float(m.state(x,m.solution.sol(x))['travel_time_s'])
        self.assertEqual(float(m.concentration(x,time_s=tau-1)),0.)
        self.assertGreater(float(m.concentration(x,time_s=tau+1)),0.)
        self.assertEqual(float(m.concentration(x,time_s=tau+m.feed.duration_s+1)),0.)
        self.assertEqual(float(m.concentration(-1,time_s=10)),0.)
        self.assertTrue(np.isnan(m.concentration(m.x[-1]+1)))

    def test_material_mass_after_shutoff(self):
        m=self.m;t=100.;tau=lambda x:float(m.state(x,m.solution.sol(x))['travel_time_s'])
        back=brentq(lambda x:tau(x)-(t-m.feed.duration_s),0,m.x[-1])
        front=brentq(lambda x:tau(x)-t,0,m.x[-1])
        # Integrated material per unit length = mass rate / effective drift velocity.
        mass=quad(lambda x:m.q/float(m.state(x,m.solution.sol(x))['speed_m_s']),back,front,epsabs=1e-7)[0]
        self.assertAlmostEqual(mass,m.q*m.feed.duration_s,places=5)

    def test_convergence_and_xy_orientation(self):
        fine=make(max_step_m=.25,rtol=1e-9,atol=1e-11)
        x=np.array([10.,50.,100.])
        np.testing.assert_allclose(self.m.concentration(x),fine.concentration(x),rtol=1e-5)
        xx,yy=np.meshgrid(x,[-3.,0.,3.]);field=self.m.concentration(xx,yy,time_s=80.)
        self.assertEqual(field.shape,(3,3));np.testing.assert_allclose(field[0],field[2])

    def test_secondary_form_roundtrip_and_source_change(self):
        values=secondary_input.display_values(D);values[('feed','rate_kg_s')]='2,0'
        d=secondary_input.collect(D,values,True,'Порог; 0,1','Точка; 40; 0; 0','30; 90')
        self.assertEqual(d['feed']['rate_kg_s'],2.)
        self.assertEqual(d['atmosphere']['temperature_k'],d['atmosphere']['surface_temperature_k'])
        self.assertEqual(secondary_input.build_model(d).q,2.)

    def test_invalid_modes(self):
        m=self.m
        with self.assertRaises(ValueError):SecondaryCloud(m.gas,GasFeed(-1,60),m.met)
        with self.assertRaises(ValueError):SecondaryCloud(m.gas,m.feed,replace(m.met,surface_temperature_k=300.))
        with self.assertRaises(ValueError):SecondaryCloud(replace(m.gas,molar_mass_kg_mol=.02),m.feed,m.met)
        with self.assertRaises(ValueError):SecondaryCloud(m.gas,m.feed,m.met,PlumeOptions(distance_m=10001))

if __name__=='__main__':unittest.main()
