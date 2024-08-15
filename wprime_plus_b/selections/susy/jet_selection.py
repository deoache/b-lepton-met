import json
import numpy as np
import awkward as ak
import importlib.resources


def select_good_jets(
    jets,
    year: str = "2017",
    jet_pt_threshold: int = 20,
    jet_abs_eta: float = 2.4,
    jet_pileup_id: str = "T",
) -> ak.highlevel.Array:
    """
    Selects and filters 'good' b-jets from a collection of jets based on specified criteria

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    year: {'2016', '2017', '2018'}
        Year for which the data is being analyzed. Default is '2017'.

    btag_working_point: {'L', 'M', 'T'}
        Working point for b-tagging. Default is 'M'.

    jet_id: https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Run_II
        Jet ID flags {1, 2, 3, 6, 7}
        For 2016 samples:
            1 means: pass loose ID, fail tight, fail tightLepVeto
            3 means: pass loose and tight ID, fail tightLepVeto
            7 means: pass loose, tight, tightLepVeto ID.
        For 2017 and 2018 samples:
            2 means: pass tight ID, fail tightLepVeto
            6 means: pass tight and tightLepVeto ID.

    jet_pileup_id: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        Pileup ID flags for pre-UL trainings {0, 4, 6, 7}. Should be applied only to AK4 CHS jets with pT < 50 GeV
        
        2016:
        0 means 000: fail all PU ID;
        1 means 001: pass loose ID, fail medium, fail tight;
        3 means 011: pass loose and medium ID, fail tight;
        7 means 111: pass loose, medium, tight ID.
        
        2017 and 2018:
        0 means 000: fail all PU ID;
        4 means 100: pass loose ID, fail medium, fail tight;
        6 means 110: pass loose and medium ID, fail tight;
        7 means 111: pass loose, medium, tight ID.

    Returns:
    --------
        An Awkward Array mask containing the selected "good" b-jets that satisfy the specified criteria.
    """
    puid_wps = {
        "2016APV": {
            "L": 1,
            "M": 3,
            "T": 7,
        },
        "2016": {
            "L": 1,
            "M": 3,
            "T": 7,
        },
        "2017": {
            "L": 4,
            "M": 6,
            "T": 7,
        },
        "2018": {
            "L": 4,
            "M": 6,
            "T": 7,
        }
    }
    jet_id = {
        "2016APV": 1,
        "2016": 1,
        "2017": 6,
        "2018": 6
    }

    # break up selection for low and high pT jets
    low_pt_jets_mask = (
        (jets.pt > jet_pt_threshold)
        & (jets.pt < 50)
        & (np.abs(jets.eta) < jet_abs_eta)
        & (jets.jetId == jet_id[year])
        & (jets.puId == puid_wps[year][jet_pileup_id])
    )
    high_pt_jets_mask = (
        (jets.pt >= 50)
        & (np.abs(jets.eta) < jet_abs_eta)
        & (jets.jetId == jet_id[year])
    )

    return ak.where(
        (jets.pt > jet_pt_threshold) & (jets.pt < 50),
        low_pt_jets_mask,
        high_pt_jets_mask,
    )
