import json
import numpy as np
import awkward as ak
import importlib.resources
from pathlib import Path
from coffea.nanoevents.methods.base import NanoEventsArray


def select_good_electrons(
    events: NanoEventsArray,
    electron_pt_threshold: int,
    electron_id_wp: str,
    electron_iso_wp: str = None,
) -> ak.highlevel.Array:
    """
    Selects and filters "good" electrons from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    electron_pt_threshold:
        Electron transverse momentum threshold

    electron_id_wp:
        Electron ID working point. Available working point for the CutBased and the MVA IDs.
        MVA: {'wp80iso', 'wp90iso', 'wp80noiso', 'wp90noiso'}
        CutBased: {'loose', 'medium', 'tight'}

    electron_iso_wp:
        Electron ISO working point {'loose', 'medium', 'tight'}. Only used for CutBased IDs or noIso MVA IDs

    Returns:
    --------
        An Awkward Array mask containing the selected "good" electrons that satisfy the specified criteria.
    """
    # electron pT threshold
    electron_pt_mask = events.Electron.pt >= electron_pt_threshold
    # electron pseudorapidity mask
    electron_eta_mask = (np.abs(events.Electron.eta) < 2.4) & (
        (np.abs(events.Electron.eta) < 1.44) | (np.abs(events.Electron.eta) > 1.57)
    )
    # electron ID and Iso masks
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
    if electron_id_wp in ["wp80iso", "wp90iso"]:
        electron_id_iso_mask = id_wps[electron_id_wp]
    else:
        electron_id_iso_mask = (id_wps[electron_id_wp]) & (iso_wps[electron_iso_wp])

    return electron_pt_mask & electron_eta_mask & electron_id_iso_mask


def select_good_muons(
    muons: ak.Array, muon_pt_threshold: int, muon_id_wp: str, muon_iso_wp: str
) -> ak.highlevel.Array:
    """
    Selects and filters "good" muons from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    muon_pt_threshold:
        Muon transverse momentum threshold

    muon_id_wp:
        Muon ID working point. Available working points for the CutBased ID {'loose', 'medium', 'tight'}

    muon_iso_wp:
        Muon ISO working point {'loose', 'medium', 'tight'}

    Returns:
    --------
        An Awkward Array mask containing the selected "good" muons that satisfy the specified criteria.
    """
    # muon pT threshold
    muon_pt_mask = muons.pt >= muon_pt_threshold

    # electron pseudorapidity mask
    muon_eta_mask = np.abs(muons.eta) < 2.4

    # muon ID mask https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    id_wps = {
        # cutbased ID working points 
        "loose": muons.looseId,
        "medium": muons.mediumId,
        "tight": muons.tightId,
    }
    muon_id_mask = id_wps[muon_id_wp]

    # muon ID and Iso mask https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection
    iso_wps = {
        "loose": muons.pfRelIso04_all < 0.25
        if hasattr(muons, "pfRelIso04_all")
        else muons.pfRelIso03_all < 0.25,
        "medium": muons.pfRelIso04_all < 0.20
        if hasattr(muons, "pfRelIso04_all")
        else muons.pfRelIso03_all < 0.20,
        "tight": muons.pfRelIso04_all < 0.15
        if hasattr(muons, "pfRelIso04_all")
        else muons.pfRelIso03_all < 0.15,
    }
    muon_iso_mask = iso_wps[muon_iso_wp]

    return (muon_pt_mask) & (muon_eta_mask) & (muon_id_mask) & (muon_iso_mask)


def select_good_taus(
    taus: ak.Array,
    tau_pt_threshold: float,
    tau_eta_threshold: float,
    tau_dz_threshold: float,
    tau_vs_jet: str,
    tau_vs_ele: str,
    tau_vs_mu: str,
    prong: int,
) -> ak.highlevel.Array:
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
    with importlib.resources.open_text("wprime_plus_b.data", "tau_wps.json") as file:
        taus_wps = json.load(file)
    # The DM is defined using the "new" decay mode reconstruction, binned in DMs 0, 1, 10, and 11.
    # https://github.com/uhh-cms/hh2bbtautau/blob/78fb359ea275e4eb2bc4cafbe238efa052d6355f/hbt/production/tau.py#L83
    prong_to_modes = {
        1: [0, 1, 2],
        2: [5, 6, 7],
        3: [10, 11],
        13: [0, 1, 2, 10, 11],
        12: [0, 1, 2, 5, 6, 7],
        23: [5, 6, 7, 10, 11],
    }
    if prong not in prong_to_modes:
        raise ValueError(
            "Invalid prong value. Please specify 1, 2, 3, 12, 13 or 23 for the prong parameter."
        )
    tau_dm = taus.decayMode
    decay_mode_mask = ak.zeros_like(tau_dm)
    for mode in prong_to_modes[prong]:
        decay_mode_mask = np.logical_or(decay_mode_mask, tau_dm == mode)

    good_taus = (
        (taus.pt > tau_pt_threshold)
        & (np.abs(taus.eta) < tau_eta_threshold)
        & (np.abs(taus.dz) < tau_dz_threshold)
        & (
            taus.idDeepTau2017v2p1VSjet
            > taus_wps["DeepTau2017"]["deep_tau_jet"][tau_vs_jet]
        )
        & (
            taus.idDeepTau2017v2p1VSe
            > taus_wps["DeepTau2017"]["deep_tau_electron"][tau_vs_ele]
        )
        & (
            taus.idDeepTau2017v2p1VSmu
            > taus_wps["DeepTau2017"]["deep_tau_muon"][tau_vs_mu]
        )
        & (decay_mode_mask)
    )
    return good_taus