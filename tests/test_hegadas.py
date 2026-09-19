"""Physics/numerics checks independent of end-to-end empirical validation."""
import unittest
from dataclasses import replace
import numpy as np
from scipy.integrate import quad
from secondary import SecondaryCloud,GasFeed,PlumeOptions,A_EDGE
from cloud import Gas,Atmosphere
from test_secondary import D,make

class HegadasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=make()
        # Synthetic numerical transition case; enlarged alpha range is NOT a validated atmosphere.
        cls.long=SecondaryCloud(Gas(**D['gas']),GasFeed(**D['feed']),
            replace(Atmosphere(**D['atmosphere']),alpha_height_limit_m=200),
            PlumeOptions(distance_m=2500,observation_s=2000,max_step_m=5,gaussian_core_fraction=.001))
        cls.long.run()

    def test_inverse_briggs(self):
        x=np.r_[0.,np.geomspace(1e-8,1e5,100)]
        np.testing.assert_allclose(self.m.inverse_sigma(self.m.sigma(x)),x,rtol=2e-14)

    def test_cross_section_integrates_actual_wind_flux(self):
        m=self.long
        for x in [1.,30.,2000.]:
            s=m.state(x,m.solution.sol(x));core=max(0.,float(s['core_half_width_m']));sy=float(s['sy_m'])
            side=2*quad(lambda y:float(m.concentration(x,y)),0,core+12*sy,points=[core])[0]
            vertical=quad(lambda z:np.exp(-(z/s['sz_m'])**m.beta)*m.met.wind_speed_m_s*(z/m.met.wind_height_m)**m.met.alpha_wind,0,np.inf)[0]
            self.assertAlmostEqual(side*vertical,m.q,places=6)

    def test_piecewise_continuity_and_gaussian_solution(self):
        m=self.long
        self.assertEqual([v['type'] for v in m.transitions],['gravity_collapse','gaussian_transition'])
        for a,b in zip(m.segments,m.segments[1:]):
            left=a[3](a[1]);right=b[3](b[0])
            np.testing.assert_allclose(left[:3],right[:3],rtol=1e-12)
            self.assertLessEqual(abs(np.sqrt(left[3]/right[3])-1),.00100001)
        tr=m.transitions[-1];x=np.linspace(tr['x_m'],m.x[-1],40)
        p=m.state(x,m.solution.sol(x))
        np.testing.assert_allclose(p['sy_m'],np.sqrt(2)*m.sigma(x+tr['virtual_origin_m']),rtol=2e-7)
        np.testing.assert_allclose(p['half_width_m'],A_EDGE*p['sy_m'],rtol=2e-7)

    def test_molar_balance_changes_at_collapse(self):
        m=self.m
        for x in [5.,30.]:
            eps=1e-3;xs=np.array([x-eps,x,x+eps]);s=m.state(xs,m.solution.sol(xs))
            total=s['total_rate_kg_s'];N=m.q/m.gas.molar_mass_kg_mol+(total-m.q)/m.met.air_molar_mass_kg_mol
            B=s['half_width_m'];_,rist=m.richardson(s)
            source=.41*m.u_star*m.beta/np.sqrt(1+.8*rist[1])/.0224
            derivative=(N[2]-N[0])/(2*eps)/(2*B[1]) if x==5 else ((N/(2*B))[2]-(N/(2*B))[0])/(2*eps)
            self.assertAlmostEqual(derivative/source,1.,places=6)

    def test_zero_initial_edge_is_regular_and_parameters_change_results(self):
        m=self.m
        self.assertEqual(m.profile['sy_m'][0],0.)
        self.assertTrue(np.all(np.isfinite(m.profile['sy_m'])))
        self.assertTrue(np.all(m.profile['core_half_width_m']>=-1e-9))
        other=SecondaryCloud(m.gas,replace(m.feed,rate_kg_s=2.),m.met,m.options);other.run()
        self.assertGreater(float(other.concentration(50)),float(m.concentration(50)))

    def test_long_run_step_convergence(self):
        m=self.long
        fine=SecondaryCloud(m.gas,m.feed,m.met,replace(m.options,max_step_m=1,rtol=1e-9,atol=1e-11));fine.run()
        x=np.array([5.,100.,1200.,2400.])
        np.testing.assert_allclose(m.concentration(x),fine.concentration(x),rtol=2e-6)
        self.assertAlmostEqual(m.transitions[-1]['x_m']/fine.transitions[-1]['x_m'],1.,places=6)

if __name__=='__main__':unittest.main()
