# selections are only define for the muon channel

import json
import numpy as np
import awkward as ak
import importlib.resources
from coffea.nanoevents.methods.base import NanoEventsArray


def select_good_electrons(
    events: NanoEventsArray,
    region: str,
) -> ak.highlevel.Array:
    """
    Selects and filters "good" electrons from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    region:

    Returns:
    --------
        An Awkward Array mask containing the selected "good" electrons that satisfy the specified criteria.
    """
    # electron pT threshold
    electron_pt_mask = events.Electron.pt >= 50

    # electron pseudorapidity mask
    electron_eta_mask = (np.abs(events.Electron.eta) < 2.4) & (
        (np.abs(events.Electron.eta) < 1.44) | (np.abs(events.Electron.eta) > 1.57)
    )

    # electron ID and Iso mask
    id_wps = {
        # mva ID working points https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
        "wp80iso": events.Electron.mvaFall17V2Iso_WP80,
        "wp90iso": events.Electron.mvaFall17V2Iso_WP90,
        "wp80noiso": events.Electron.mvaFall17V2noIso_WP80,
        "wp90noiso": events.Electron.mvaFall17V2noIso_WP90,
        # cutbased ID working points https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
        "loose": events.Electron.cutBased == 2,
        "medium": events.Electron.cutBased == 3,
        "tight": events.Electron.cutBased == 4,
    }
    iso_wps = {
        # https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection
        "loose": events.Electron.pfRelIso04_all < 0.25
        if hasattr(events.Electron, "pfRelIso04_all")
        else events.Electron.pfRelIso03_all < 0.25,
        "medium": events.Electron.pfRelIso04_all < 0.20
        if hasattr(events.Electron, "pfRelIso04_all")
        else events.Electron.pfRelIso03_all < 0.20,
        "tight": events.Electron.pfRelIso04_all < 0.15
        if hasattr(events.Electron, "pfRelIso04_all")
        else events.Electron.pfRelIso03_all < 0.15,
    }

    return electron_pt_mask & electron_eta_mask & id_wps["wp90iso"]


def select_good_muons(
    events: NanoEventsArray, region: str
) -> ak.highlevel.Array:
    """
    Selects and filters "good" muons from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    region:

    Returns:
    --------
        An Awkward Array mask containing the selected "good" muons that satisfy the specified criteria.
    """
    # muon pT threshold
    muon_pt_mask = events.Muon.pt >= 35

    # electron pseudorapidity mask
    muon_eta_mask = np.abs(events.Muon.eta) < 2.4

    # muon ID mask https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    id_wps = {
        # cutbased ID working points 
        "loose": events.Muon.looseId,
        "medium": events.Muon.mediumId,
        "tight": events.Muon.tightId,
    }
    if region in ["A", "D"]:
        muon_id_mask = id_wps["tight"]
    else:
        muon_id_mask = (id_wps["medium"]) & (~id_wps["tight"])

    # muon ID and Iso mask https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection
    iso_wps = {
        "loose": events.Muon.pfRelIso04_all < 0.25
        if hasattr(events.Muon, "pfRelIso04_all")
        else events.Muon.pfRelIso03_all < 0.25,
        "medium": events.Muon.pfRelIso04_all < 0.20
        if hasattr(events.Muon, "pfRelIso04_all")
        else events.Muon.pfRelIso03_all < 0.20,
        "tight": events.Muon.pfRelIso04_all < 0.15
        if hasattr(events.Muon, "pfRelIso04_all")
        else events.Muon.pfRelIso03_all < 0.15,
    }

    return muon_pt_mask & muon_eta_mask & muon_id_mask & iso_wps["tight"]


def select_good_taus(events: NanoEventsArray) -> ak.highlevel.Array:
    """
    Selects and filters "good" taus from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    Returns:
    --------
        An Awkward Array mask containing the selected "good" taus that satisfy the specified criteria.
    """
    return (
        (events.Tau.pt > 20)
        & (np.abs(events.Tau.eta) < 2.3)
        & (events.Tau.dz < 0.2)
        & (events.Tau.idDeepTau2017v2p1VSjet > 8)
        & (events.Tau.idDeepTau2017v2p1VSe > 8)
        & (events.Tau.idDeepTau2017v2p1VSmu > 1)
    )