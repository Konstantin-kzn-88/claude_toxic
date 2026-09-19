import copy
import json
from pathlib import Path
import unittest
from input_data import display_values, collect
from cloud import Gas,Vessel,Atmosphere,Options,PrimaryCloud

D=json.loads((Path(__file__).resolve().parents[1]/'examples/chloromethane.json').read_text(encoding='utf-8'))

def gather(values=None, confirmed=True, thresholds=None, times='10; 30; 60'):
    return collect(D,values or display_values(D),confirmed,thresholds if thresholds is not None else '\n'.join(f'{name}; {c:.16g}' for name,c in zip(D['threshold_labels'],D['thresholds_kg_m3'])), 'Точка; 50; -10; 1,5',times)

class InputTests(unittest.TestCase):
    def test_display_units_and_roundtrip(self):
        v=display_values(D)
        self.assertAlmostEqual(float(v[('vessel','temperature_k')]),18.)
        self.assertAlmostEqual(float(v[('vessel','pressure_abs_pa')]),.101325)
        self.assertAlmostEqual(float(v[('gas','molar_mass_kg_mol')]),51.)
        result=gather(v)
        for section in ('gas','vessel','atmosphere','options'):
            for key,value in D[section].items():
                if isinstance(value,float):self.assertAlmostEqual(result[section][key],value)
                else:self.assertEqual(result[section][key],value)
        self.assertEqual(result['receptors'][0]['z_m'],1.5)
        self.assertEqual(result['receptors'][0]['y_m'],-10.)

    def test_modified_volume_changes_mass_without_mutating_source(self):
        baseline=copy.deepcopy(D);v=display_values(D)
        v[('vessel','volume_m3')]='1000,0'
        d=gather(v)
        m=PrimaryCloud(Gas(**d['gas']),Vessel(**d['vessel']),Atmosphere(**d['atmosphere']),Options(**d['options']))
        self.assertAlmostEqual(m.initial['mass_kg'],4269.421919491348/2)
        self.assertEqual(D,baseline)

    def test_invalid_input_rejected(self):
        for value in ['nan','inf','abc']:
            v=display_values(D);v[('vessel','volume_m3')]=value
            with self.assertRaises(ValueError):gather(v)
        with self.assertRaises(ValueError):gather(confirmed=False)
        with self.assertRaises(ValueError):gather(thresholds='Порог; -1')
        with self.assertRaises(ValueError):gather(times='601')
        v=display_values(D);v[('maps_xy','nx')]='30,5'
        with self.assertRaises(ValueError):gather(v)

if __name__=='__main__':unittest.main()
