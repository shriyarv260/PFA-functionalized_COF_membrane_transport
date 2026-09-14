import unittest
import numpy as np
from pfas_cof.analysis import msd, diffusion, wham, residence_events, water_permeance, hydration_counts


class AnalysisTests(unittest.TestCase):
    def test_diffusion_units(self):
        t = np.arange(30.)
        self.assertAlmostEqual(diffusion(t, 6*.002*t, 2, 20, 3)['D_m2_s'], 2e-9)

    def test_msd_stationary(self):
        _, y = msd(np.zeros((10, 2, 3)), 1, 5)
        np.testing.assert_array_equal(y, 0)

    def test_brownian(self):
        rng = np.random.default_rng(7)
        x = np.cumsum(rng.normal(0,np.sqrt(.004),(1200,100,3)),axis=0)
        t,y = msd(x,1,70)
        self.assertLess(abs(diffusion(t,y,5,60,3)['D_nm2_ps']/.002-1), .12)

    def test_harmonic_wham(self):
        rng = np.random.default_rng(8)
        c = np.linspace(-1.5,1.5,16); k = 30.; a = 4.
        s = [rng.normal(k*z/(k+a),np.sqrt(.008314462618*298.15/(k+a)),8000) for z in c]
        z,f,d = wham(s,c,np.full(16,k),np.linspace(-3,3,121))
        exact = .5*a*z*z; exact -= exact.min()
        mask = np.abs(z)<1.2
        self.assertLess(np.sqrt(np.mean((f[mask]-exact[mask])**2)),.25)
        self.assertGreater(d['minimum_overlap'], 0)

    def test_disconnected_rejected(self):
        with self.assertRaises(ValueError): wham([[-1.], [1.]], [-1,1], [30,30], np.linspace(-2,2,10))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): msd(np.zeros((3,1,3)),1,3)
        with self.assertRaises(ValueError): diffusion(np.arange(5),-np.arange(5),0,4,1)
        with self.assertRaises(ValueError): wham([[0],[1]],[0,1],[1,1],[0,1,3])

    def test_censoring(self):
        events = residence_events(np.array([[1],[1],[0],[1]],dtype=bool),2)
        self.assertEqual(events[0], {'particle':0,'duration_ps':4.,'left_censored':True,'right_censored':False})
        self.assertTrue(events[1]['right_censored'])

    def test_permeance(self):
        self.assertAlmostEqual(water_permeance(10,10,100,20,1)['flux_m_s'],0)
        self.assertAlmostEqual(water_permeance(1,0,1,1,1,1)['permeance_m_Pa_s'],.01)

    def test_periodic_hydration(self):
        np.testing.assert_array_equal(hydration_counts([[.1,0,0]],[[2.9,0,0],[1,1,1]],[3,3,3]),[1])


if __name__ == '__main__': unittest.main()
