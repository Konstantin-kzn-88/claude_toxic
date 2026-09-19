import unittest
import numpy as np
from test_cloud import model
from zones import section, envelope, receptor

class ZoneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = model(120, output_step_s=.1)
        cls.m.run()

    def test_section_boundary_matches_concentration(self):
        m = self.m
        centre_x = m.state(m.solution.sol(20))['centre_x_m']
        threshold = float(m.concentration(20, centre_x, 0., 1.)) / 2
        d = section(m, 20, threshold, 1.)
        self.assertIsNotNone(d)
        self.assertAlmostEqual(float(m.concentration(20, d['x_max_m'], 0., 1.)), threshold, places=10)
        self.assertGreater(float(m.concentration(20, d['centre_x_m'], d['radius_m']*.99, 1.)), threshold)
        self.assertLess(float(m.concentration(20, d['centre_x_m'], d['radius_m']*1.01, 1.)), threshold)

    def test_initial_geometry_and_empty_zone(self):
        d = section(self.m, 0., .2)
        self.assertAlmostEqual(d['radius_m'], self.m.initial['radius_m'])
        e = envelope(self.m, 10.)
        self.assertIsNone(e['x_max_m'])
        self.assertEqual(e['footprint_area_m2'], 0.)

    def test_envelope_nesting_and_grid_convergence(self):
        large = envelope(self.m, .1, spatial_points=2001)
        small = envelope(self.m, .2, spatial_points=2001)
        fine = envelope(self.m, .1, spatial_points=4001)
        self.assertGreaterEqual(large['x_max_m'], small['x_max_m'])
        self.assertGreaterEqual(large['max_width_m'], small['max_width_m'])
        self.assertGreaterEqual(large['footprint_area_m2'], small['footprint_area_m2'])
        self.assertLess(abs(large['footprint_area_m2']/fine['footprint_area_m2'] - 1), .001)
        # Every instantaneous circle's area is bounded by the union area.
        self.assertGreaterEqual(large['footprint_area_m2'], max(d['area_m2'] for d in large['disks'])*.999)

    def test_receptor_crossings(self):
        obs = receptor(self.m, 50., thresholds=[.2])
        band = obs['threshold_exposure'][0]
        self.assertIsNotNone(band['arrival_s'])
        self.assertIsNotNone(band['last_exit_s'])
        self.assertGreater(band['duration_in_computed_interval_s'], 0)
        for a,b in band['intervals_s']:
            self.assertAlmostEqual(float(self.m.concentration(a, 50.)), .2, places=7)
            self.assertAlmostEqual(float(self.m.concentration(b, 50.)), .2, places=7)

    def test_open_exposure_not_reported_as_completed(self):
        obs = receptor(self.m, 50., thresholds=[1e-10])
        self.assertTrue(obs['threshold_exposure'][0]['active_at_end'])
        self.assertIsNone(obs['threshold_exposure'][0]['last_exit_s'])

    def test_numeric_booleans_and_thresholds_rejected(self):
        from dataclasses import replace
        from cloud import PrimaryCloud
        with self.assertRaises(ValueError):
            PrimaryCloud(self.m.gas, replace(self.m.vessel, volume_m3=True), self.m.met)
        for value in [0., -1., float('nan'), True]:
            with self.assertRaises(ValueError): envelope(self.m, value)

if __name__ == '__main__': unittest.main()
