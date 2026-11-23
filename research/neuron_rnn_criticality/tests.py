"""Small scientific checks for the spike adapter, covariance and surrogate."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np

from analyse_bugnon import archive_wave_stats, paired_timescales, parameter_evidence, phase_spike_counts
from join_fits import verify_fit_input
from physiology import hmm,merge_short,periods
from analyse_snn import components, pure_episodes, shifted


class ScientificChecks(unittest.TestCase):
    def test_archived_slow_waves_threshold_and_boundaries(self):
        with TemporaryDirectory() as temp:
            root=Path(temp)
            (root/"df_sw.csv").write_text(
                "subnet,startTime,endTime,maxPeak_Amp,duration,halfSlopetoFirstPeakTime,halfSlopefromLastPeakTime\n"
                "full,100,400,1,300,170,310\n"
                "full,410,810,2,400,510,710\n"
                "full,850,1200,3,350,920,1100\n"
                "local,100,400,4,300,170,310\n"
            )
            result=archive_wave_stats(root,100,1000,{"strict":1,"lower_threshold":0})
        self.assertEqual(result["candidate_waves"],2)
        strict=result["threshold_sensitivity"]["strict"]
        self.assertEqual(strict["above_threshold"],1)
        self.assertEqual(strict["median_down_half_slope_s"],.2)
        self.assertEqual(result["threshold_sensitivity"]["lower_threshold"]["largest_quartile"],1)

    def test_paired_timescales_require_a_valid_estimate_in_both_windows(self):
        a=np.array([[.2,.4],[np.nan,.2],[.1,np.nan]])
        b=np.array([[.4,.6],[.3,.1],[.2,.4]])
        results=paired_timescales(a,b,np.array(["PC","PC","BC"]))
        pc=results["PC"]
        self.assertEqual(pc["exponential"]["n_paired"],1)
        np.testing.assert_allclose(pc["exponential"]["change_s_p10_p50_p90"],[.2,.2,.2])
        self.assertEqual(pc["integrated"]["n_paired"],2)
        self.assertAlmostEqual(pc["integrated"]["change_s_p10_p50_p90"][1],.05)
        self.assertIsNone(results["BC"]["integrated"]["change_s_p10_p50_p90"])

    def test_bugnon_recording_flags_come_from_verified_parameters(self):
        tree=dict(simulation={"params":{"sessions":["transition","sleep"],"output_dir":"original"}},
                  session_models={"transition":{"params":{"record":False}},"sleep":{"params":{}}},
                  network={k:{} for k in ("layers","projection_models","synapse_models","topology")},kernel={})
        with TemporaryDirectory() as temp:
            frozen,executed=(Path(temp)/f"{name}.yml" for name in ("frozen","executed"))
            frozen.write_text(json.dumps(tree))
            tree["simulation"]["params"]["output_dir"]="local"
            executed.write_text(json.dumps(tree))
            result=parameter_evidence(frozen,executed)
            self.assertEqual(result["session_recorded"],{"00_transition":False,"01_sleep":True})
            tree["session_models"]["transition"]["params"]["record"]=True
            executed.write_text(json.dumps(tree))
            with self.assertRaises(ValueError):
                parameter_evidence(frozen,executed)

    def test_bugnon_sessions_exclude_the_transition_from_stable_counts(self):
        times=np.array([2500,2999.8,3000,4999.8,5000,5000.2,125000])
        sessions={"wake":(1000,3000),"transition":(3000,5000),"sleep":(5000,125000)}
        self.assertEqual(phase_spike_counts(times,sessions),{"wake":2,"transition":2,"sleep":2})

    def test_fit_receipt_preserves_counts_and_segments(self):
        rng=np.random.default_rng(7)
        blocks=[rng.poisson(4,60).tolist(),rng.poisson(8,70).tolist()]
        population=np.concatenate(blocks)
        centered=[b-population.mean() for b in map(np.asarray,blocks)]
        acf=[sum(x[:-k]@x[k:] for x in centered)/((len(population)-2*k)*population.var())
             for k in range(1,21)]
        fit=dict(segment_bins=[60,70],duration_s=5.2,rate_hz=population.mean()/.04,
                 count_variance=population.var(),observed_acf=acf)
        row={"record":{"id":"segments"},"fit":dict(fit,ar_order=5),"long_reach":dict(fit,ar_order=10)}
        verify_fit_input(row,blocks)
        row["long_reach"]["segment_bins"]=[130]
        with self.assertRaises(AssertionError):
            verify_fit_input(row,blocks)

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
