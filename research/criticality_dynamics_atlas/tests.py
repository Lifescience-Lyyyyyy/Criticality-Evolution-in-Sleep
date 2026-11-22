import unittest

import numpy as np

from metrics import ar_modes, lag_cov, mr, pr, shift, unit_metrics
from theory import ar, branch
from analyse import fold, selected


class ScientificChecks(unittest.TestCase):
    def test_ar_and_mr(self):
        rng=np.random.default_rng(42)
        for m in (.55,.85,.95):
            y=ar(m,80000,rng.normal(size=80000))[1000:]
            rho=lag_cov([y],4)
            np.testing.assert_allclose(rho[1:]/rho[0],m**np.arange(1,5),atol=.018)
            x=branch(m,np.full(80000,(1-m)*10),rng)[1000:]
            for f in (1.,.3):
                z=rng.binomial(x.astype(int),f).astype(float)
                q,_,_=mr([z],.04)
                self.assertTrue(q["valid"])
                self.assertLess(abs(q["m"]-m),.04)
            self.assertAlmostEqual(ar_modes([m],.04)[0]["tau"],-.04/np.log(m))

    def test_segment_and_covariance(self):
        rng=np.random.default_rng(51)
        x=rng.poisson(2,(600,5)).astype(float)
        y=rng.poisson(8,(650,5)).astype(float)
        c=lag_cov([x,y],1)
        z=[a-a.mean(0) for a in (x,y)]
        np.testing.assert_allclose(c[1],sum((a[1:]*a[:-1]).sum(0) for a in z)/1248)
        r=unit_metrics([x,y],.04,99,True)
        np.testing.assert_allclose(np.add(r["cov_auto"],r["cov_cross"]),r["cov_total"],atol=1e-12)
        self.assertAlmostEqual(pr(np.eye(7)),7)
        self.assertAlmostEqual(pr(np.ones((7,7))),1)

    def test_separate_windows_and_time(self):
        a=np.arange(1510)
        b=np.arange(1515)+100000
        self.assertEqual([len(x) for x in fold([a,b],5)],[302,303])
        self.assertEqual(fold([a,b],5)[1][0],500010)
        windows=list(selected([a,b],.04))
        self.assertEqual([(s,k,start) for s,k,start,_ in windows],[(60,0,0),(60,1,0)])
        np.testing.assert_equal(windows[1][-1],b[:1500])

    def test_animal_is_unit(self):
        from summarize import aggregate_deltas
        ds=[dict(animal="a",contrast="NREM_Wake",d2=i) for i in (2,4,6)]
        ds.append(dict(animal="b",contrast="NREM_Wake",d2=-2))
        g=aggregate_deltas(ds,["d2"])
        self.assertEqual([(r["session_n"],r["d2"]) for r in g],[(3,4),(1,-2)])
        self.assertEqual(np.mean([r["d2"] for r in g]),1)

    def test_shift(self):
        rng=np.random.default_rng(84)
        x=rng.integers(0,8,size=(1500,8))
        y=shift(x,rng)
        np.testing.assert_equal(np.sort(x,axis=0),np.sort(y,axis=0))
        np.testing.assert_allclose(abs(np.fft.rfft(x,axis=0)),abs(np.fft.rfft(y,axis=0)),atol=1e-10)


if __name__=="__main__":
    unittest.main()
