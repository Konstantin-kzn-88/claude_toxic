import unittest
import numpy as np
from test_cloud import model
from maps_xy import concentration_grid

class MapsXYTests(unittest.TestCase):
    def test_xy_orientation_symmetry_and_maximum(self):
        m=model(10,output_step_s=.2);m.run()
        x=np.array([-10.,0.,20.,40.]);y=np.array([-5.,0.,5.])
        field=concentration_grid(m,x,y,0.,5.)
        self.assertEqual(field.shape,(3,4))
        self.assertAlmostEqual(field[2,3],float(m.concentration(5.,40.,5.,0.)))
        np.testing.assert_allclose(field[0],field[2])
        maximum=concentration_grid(m,x,y,0.)
        for j in range(3):
            for i in range(4):
                expected=max(float(m.concentration(t,x[i],y[j],0.)) for t in m.times)
                self.assertAlmostEqual(maximum[j,i],expected,places=12)

if __name__=='__main__':unittest.main()
