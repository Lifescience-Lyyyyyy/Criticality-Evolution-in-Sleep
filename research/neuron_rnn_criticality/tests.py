"""Small scientific checks for the spike adapter, covariance and surrogate."""
import unittest
import numpy as np

from physiology import hmm,merge_short,periods
from analyse_snn import components, pure_episodes, shifted


class ScientificChecks(unittest.TestCase):
    def test_detector_and_intervals(self):
        rng=np.random.default_rng(4)
        actual=np.repeat([False,True,False,True,False],100)
        path,rate=hmm(rng.poisson(np.where(actual,4.5,.05)))
        self.assertGreater(np.mean(path==actual),.98)
        self.assertLess(rate[0],rate[1])
        actual[125:128]=False
        self.assertEqual(len(periods(merge_short(actual,5))),5)

    def test_covariance_identity_direct(self):
        rng=np.random.default_rng(5)
        x=rng.integers(0,7,size=(180,11))
        _,a,b,t=components(x)
        np.testing.assert_allclose(a+b,t,atol=1e-12)
        c=x-x.mean(0)
        for k in (0,1,4,33):
            cc=c[:len(c)-k].T @ c[k:] / (len(c)-k)
            np.testing.assert_allclose(np.trace(cc),a[k])
            np.testing.assert_allclose(cc.sum()-np.trace(cc),b[k])

    def test_boundary_and_shift(self):
        state=np.r_[np.zeros(45),np.ones(47),np.zeros(108)].astype(bool)
        boundaries=pure_episodes(state)
        self.assertEqual(boundaries,[(0,0,11),(2,11,12),(1,12,23),(0,23,50)])
        x=np.random.default_rng(6).poisson(1,size=(50,7))
        y,info=shifted(x,boundaries,111)
        self.assertEqual(info["max_cyclic_autocorrelation_error"],0)
        np.testing.assert_equal(x[11],y[11])
        for _,a,b in boundaries:
            np.testing.assert_equal(np.sort(x[a:b],axis=0),np.sort(y[a:b],axis=0))
        self.assertFalse(np.array_equal(x,y))


if __name__=="__main__":
    unittest.main()
