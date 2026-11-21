"""Deterministic window and legacy-compatibility check for the shared preparer."""
import json
from unittest import TestCase, main

import numpy as np

from build_qc_fixture import ROOT
from prepare_crcns_sleep import prepare_group
from crcns_adapter import SleepRecording


class PairwiseDurationTest(TestCase):
    def test_pairwise_does_not_depend_on_unrelated_states(self):
        cfg = json.loads((ROOT/"config/sleep_dandi.json").read_text())
        ranges = {"Wake":[(0,123),(140,202)],"NREM":[(220,411)],"REM":[(430,499),(510,575)]}
        spikes = tuple(np.arange(i+.004,800,.053+i*.01) for i in range(6))
        args = ("mouse-a","session-1","000939","PoSub",tuple(map(str,range(6))),spikes)
        rest = ([(0,800)], {"seconds":[0,800]}, {"path":"test","sha256":"test","asset_id":"test"}, {}, "REST")
        pair = prepare_group([SleepRecording(*args,ranges,*rest)],cfg,"pairwise")
        self.assertEqual([(m["variant"],m["state"]) for m,_ in pair],
                         [("equal_duration_NW","Wake"),("equal_duration_NW","NREM")])
        self.assertEqual(pair[0][0]["intervals_s"],[[0,60],[60,120],[140,200]])
        self.assertEqual(pair[1][0]["intervals_s"],[[220,280],[280,340],[340,400]])
        for m,counts in pair:
            self.assertEqual(m["unit_ids"],list(args[4]))
            self.assertEqual([len(x[0]) for x in counts],[1500,1500,1500])
        legacy = prepare_group([SleepRecording(*args,ranges,*rest)],cfg,"states")
        self.assertFalse(any(m["variant"]=="equal_duration" for m,_ in legacy))
        longer={**ranges,"REM":[(430,671)]}
        legacy = prepare_group([SleepRecording(*args,longer,*rest)],cfg,"states")
        self.assertEqual({m["state"]:m["intervals_s"] for m,_ in legacy if m["variant"]=="equal_duration"},
            {"Wake":[[0,60],[60,120],[140,200]],"NREM":[[220,280],[280,340],[340,400]],
             "REM":[[430,490],[490,550],[610,670]]})
        only_rw={k:v for k,v in longer.items() if k != "NREM"}
        rw=prepare_group([SleepRecording(*args,only_rw,*rest)],cfg,"pairwise")
        self.assertEqual({m["variant"] for m,_ in rw},{"equal_duration_RW"})


if __name__=="__main__":
    main()
