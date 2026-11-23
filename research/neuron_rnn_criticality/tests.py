"""Small scientific checks for the spike adapter, covariance and surrogate."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np

from analyse_bugnon import (archive_wave_stats, fixed_subsamples, load_prepared_condition,
                            local_transition_synchrony, paired_timescales, paper_physiology_gates,
                            paper_run_evidence, parameter_evidence, phase_spike_counts,
                            save_prepared_condition, transition_crossings, verify_prepared_spikes)
from run_bugnon_paper import COMMIT as PAPER_COMMIT, session_overrides
from join_fits import matching_directions, verify_fit_input
from physiology import digest,hmm,merge_short,periods
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

    def test_prepared_condition_roundtrip_requires_its_checksum_and_condition(self):
        counts=np.tile(np.arange(1500)%5,(7,1))
        summary={"receipt":{"source_commit":PAPER_COMMIT,"condition":"wake",
                            "duration_seconds":120,"neuron_count_recorded":2},
                 "decomposition":{"total":[counts[0].var()]}}
        prepared=dict(summary=summary,counts=counts,neuron_ids=np.array([11,12]),
                      neuron_types=np.array(["PC","BC"]),neuron_layers=np.array(["L4","L4"]),
                      rates_hz=np.array([.5,2.]),timescales=np.full((2,2),np.nan))
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/"wake.npz"
            save_prepared_condition(path,prepared)
            out=load_prepared_condition(path,"wake")
            np.testing.assert_array_equal(out["counts"],counts)
            np.testing.assert_array_equal(out["timescales"],prepared["timescales"])
            legacy=Path(tmp)/"legacy-wake.npz"
            save_prepared_condition(legacy,dict(prepared,counts=counts.astype(float)))
            converted=load_prepared_condition(legacy,"wake")["counts"]
            np.testing.assert_array_equal(converted,counts)
            self.assertEqual(converted.dtype,np.int64)
            fractional=Path(tmp)/"fractional-wake.npz"
            save_prepared_condition(fractional,dict(prepared,counts=counts+.5))
            with self.assertRaisesRegex(ValueError,"not exact nonnegative integers"):
                load_prepared_condition(fractional,"wake")
            with self.assertRaisesRegex(ValueError,"completed 120 s paper-linked sleep"):
                load_prepared_condition(path,"sleep")
            path.with_suffix(".sha256").write_text("0"*64)
            with self.assertRaisesRegex(ValueError,"Checksum mismatch"):
                load_prepared_condition(path,"wake")

    def test_native_archive_rebins_every_transferred_series(self):
        types=np.repeat(["PC","BC","MC"],[395,79,26])
        ids=np.arange(len(types))
        layers=np.repeat("L4",len(ids))
        offsets=np.arange(len(ids)+1)
        times=np.arange(len(ids))*.04
        counts=np.zeros((7,1500),dtype=np.int64)
        counts[0,ids]=1
        for i,(_,_,cols) in enumerate(fixed_subsamples(types),1):
            counts[i,cols]=1
        with TemporaryDirectory() as temp:
            path=Path(temp)/"spikes.npz"
            np.savez_compressed(path,neuron_ids=ids,neuron_types=types,neuron_layers=layers,
                                spike_times_offsets=offsets,spike_times_seconds=times)
            data=dict(summary={"receipt":{"spike_record_sha256":digest(path),"duration_seconds":120}},
                      neuron_ids=ids,neuron_types=types,neuron_layers=layers,rates_hz=np.full(len(ids),1/120),
                      counts=counts)
            self.assertEqual(verify_prepared_spikes(path,data)["compared_series"],7)
            data["counts"][-1,0]+=1
            with self.assertRaises(AssertionError):
                verify_prepared_spikes(path,data)

    def test_paper_fit_gate_requires_sleep_specific_physiology(self):
        types=np.array(["PC","PC","PC","BC","BC","MC"])
        rates=(np.array([0,0,2,10,12,2]),np.array([0,1,3,11,13,3]))
        corr={k:{"pearson":.8} for k in ("PC","BC")}
        episode=dict(up_n=2,down_n=2,E=dict(up_hz=2,down_hz=.2),I=dict(up_hz=10,down_hz=1))
        rows={"wake":{"lfp_proxy":{"statistics":{"band_power":1}}},
              "sleep":{"lfp_proxy":{"statistics":{"band_power":3}},"episode_proxy":episode,
                       "voltage":{"Vp_L4_PC":{"up_minus_down_mv":3}}}}
        self.assertTrue(all(paper_physiology_gates(rows,rates,types,corr).values()))
        episode["E"]["down_hz"]=episode["E"]["up_hz"]
        gates=paper_physiology_gates(rows,rates,types,corr)
        self.assertFalse(gates["sleep_bouts"])
        self.assertTrue(gates["slow_modulation"])
        for state,stage in (("wake","noise"),("sleep","sleep")):
            plan=session_overrides(state)
            self.assertEqual(plan["session_models"][stage]["params"]["simulation_time"],120000)
            self.assertNotIn("grating",plan["simulation"]["params"]["sessions"])
        with self.assertRaises(ValueError):
            paper_run_evidence({"condition":"network-check"},"sleep",None,None)

    def test_voltage_halfway_crossings_pair_the_same_cells(self):
        traces=np.array([[4,4,3,2,0,0,1,3,4],
                         [4,3,1,0,0,1,3,4,4],
                         [0,0,0,0,0,0,0,0,0],
                         [4,4,3,2,np.nan,1,2,3,4]]).T
        down,up=transition_crossings(traces,4)
        np.testing.assert_array_equal(down,[4,2])
        np.testing.assert_array_equal(up,[7,6])
        self.assertEqual(transition_crossings(traces[:,2:],4)[0].size,0)

    def test_bugnon_synchrony_inputs_survive_archival(self):
        timeline=10*np.arange(400)
        wave=np.cos(2*np.pi*timeline/1000)
        voltage=[(-65-wave)[:,None],(-65-np.roll(wave,1))[:,None]]
        currents=np.vstack((np.zeros_like(wave),2*wave))
        coverage=np.full(len(wave),2.)
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/"local_cortical_voltage_10ms.npz"
            result=local_transition_synchrony(currents,coverage,voltage,[11,12],8000,path)
            self.assertEqual(result["input_archive_sha256"],digest(path))
            self.assertGreater(result["n_valid_paired_waves"],0)
            with np.load(path,allow_pickle=False) as saved:
                np.testing.assert_array_equal(saved["voltage_mv"],np.column_stack(voltage))
                np.testing.assert_array_equal(saved["neuron_ids"],[11,12])
                np.testing.assert_array_equal(saved["time_ms"],timeline+8000)
                np.testing.assert_array_equal(saved["ampa_current_sum"],currents[0])
                np.testing.assert_array_equal(saved["gaba_current_sum"],currents[1])
                np.testing.assert_array_equal(saved["sample_counts"],coverage)

    def test_bugnon_sessions_exclude_the_transition_from_stable_counts(self):
        times=np.array([2500,2999.8,3000,4999.8,5000,5000.2,125000])
        sessions={"wake":(1000,3000),"transition":(3000,5000),"sleep":(5000,125000)}
        self.assertEqual(phase_spike_counts(times,sessions),{"wake":2,"transition":2,"sleep":2})

    def test_matched_directions_keep_unavailable_fits_and_technical_replicates(self):
        pairs=[dict(n_units=10,sample_seed=seed,delta_sleep_minus_wake={"fit":value,"ar10":None})
               for seed,value in zip((2401,2402,2403),(.2,-.1,None),strict=True)]
        empirical={f"N10_w{i}_60s":{"animal_deltas":{"a":[-.2,0],"b":[.1,None]}} for i in range(1,4)}
        result=matching_directions(pairs,empirical)
        self.assertEqual(result["model_networks"],1)
        self.assertEqual(result["model_technical_subsamples"]["ar5"]["unavailable"],1)
        self.assertEqual(result["model_technical_subsamples"]["ar10"]["unavailable"],3)
        self.assertEqual(result["reference_animals"]["w1_60s"]["ar10"]["equal"],1)

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
