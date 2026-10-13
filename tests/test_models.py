"""Independent mathematical checks, conservation and data reproducibility."""
import csv
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from logos_pkpd import finite, tail, synthetic, static
from logos_pkpd.periodic import cycle, constant

class ModelTests(unittest.TestCase):
    def test_finite_and_conservation(self):
        result=finite.run()
        self.assertEqual(result['timing_starts'],24)
        objectives=[r['J'] for r in result['rows']]
        self.assertTrue(np.all(np.diff(objectives)>0))
        self.assertAlmostEqual(objectives[-1],.679552690590195,places=10)
        for r in result['rows']:
            dose=sum(r['doses'])
            if r['regimen']=='expanded':dose+=r['concentration']*r['maintenance']
            self.assertAlmostEqual(dose,1.,places=8)
        for k,V,E,D in [(1,1,1,1),(.4,2,.7,3),(2,.5,4,.8)]:
            threshold=np.log(2+E*V/D)/k
            for T in [.5*threshold,threshold,2*threshold]:
                p=finite.plateau(D,T,k,V,E)
                self.assertAlmostEqual(V*p['concentration']*(1+k*p['maintenance']),D,places=10)
                self.assertGreaterEqual(p['maintenance'],0.)

    def test_tail_against_ode(self):
        for rho in [0,.2,2.,3.]:
            for z0 in [0,.3]:
                t=np.linspace(0,20,401)
                sol=solve_ivp(lambda t,y:[.7*np.exp(-2*t)-rho*y[0]],(0,20),[z0],method="DOP853",t_eval=t,rtol=1e-11,atol=1e-13)
                np.testing.assert_allclose(tail.injury(t,c0=1,z0=z0,alpha=.7,k=1,rho=rho),sol.y[0],rtol=0,atol=2e-10)
        p=tail.peak();self.assertAlmostEqual(p['time'],1.2792139405522476,places=10)
        self.assertGreater(p['value'],.1);self.assertLess(float(tail.injury(.05)),.1)
        self.assertFalse(tail.peak(rho=0)['attained'])
        self.assertTrue(np.all(np.isfinite(tail.injury([0,1,1000],rho=3))))

    def test_periodic_null(self):
        reference=constant(False,False)['J']
        for T in [.08,1.,3.,6.,12.]:
            result=cycle(T,False,False)
            self.assertLessEqual(result['J'],reference+1e-8)
            self.assertLess(result['closure'],1e-10)

    def test_synthetic_determinism(self):
        with tempfile.TemporaryDirectory() as directory:
            a,b=Path(directory)/'a.csv',Path(directory)/'b.csv'
            synthetic.generate(a);synthetic.generate(b)
            self.assertEqual(a.read_bytes(),b.read_bytes())
            with a.open() as f:
                rows=list(csv.DictReader(f))
            self.assertEqual(len(rows),600)
            self.assertTrue(all(r['is_synthetic']=='true' for r in rows))

    def test_hill_envelope(self):
        for h in [2.,3.,4.]:
            for cap in [.5,3.]:
                mean=np.linspace(0,cap,501)
                env=static.envelope(mean,cap,h)
                self.assertTrue(np.all(env>=static.hill(mean,h)-1e-14))
                self.assertLessEqual(np.diff(env,n=2).max(),1e-12)
                support=min(cap,(h-1)**(1/h))
                slope=static.hill(support,h)/support
                self.assertTrue(np.all(static.hill(mean,h)<=slope*mean+1e-12))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):finite.plateau(dose=0)
        with self.assertRaises(ValueError):finite.objective([1,1],[1,0])
        with self.assertRaises(ValueError):tail.injury([-1])

if __name__=='__main__':unittest.main()
