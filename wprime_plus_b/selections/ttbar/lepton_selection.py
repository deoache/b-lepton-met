import numpy as np
import awkward as ak
from coffea.nanoevents.methods.base import NanoEventsArray


def select_good_electrons(
    events: NanoEventsArray, channel: str, lepton_flavor: str
) -> ak.highlevel.Array:
    """
    Selects and filters "good" electrons from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    channel: {'1b1l', '1b1e1mu', '2b1l'}
        The channel of the event.

    lepton_flavor: {'ele', 'mu'}
        The flavor of the lepton

    Returns:
    --------
        An Awkward Array mask containing the selected "good" electrons that satisfy the specified criteria.
    """
    if channel != "1b1e1mu":
        good_electron_pt = 55 if lepton_flavor == "ele" else 30
    else:
        good_electron_pt = 55 if lepton_flavor == "mu" else 30
    return (
        (events.Electron.pt >= good_electron_pt)
        & (np.abs(events.Electron.eta) < 2.4)
        & ((np.abs(events.Electron.eta) < 1.44) | (np.abs(events.Electron.eta) > 1.57))
        & (
            events.Electron.mvaFall17V2Iso_WP80
            if lepton_flavor == "ele"
            else events.Electron.mvaFall17V2Iso_WP90
        )
        & (
            events.Electron.pfRelIso04_all < 0.25
            if hasattr(events.Electron, "pfRelIso04_all")
            else events.Electron.pfRelIso03_all < 0.25
        )
    )


def select_good_muons(events: NanoEventsArray) -> ak.highlevel.Array:
    """
    Selects and filters "good" muons from a collection of events based on specified criteria.

    Parameters:
    -----------
    events:
        A collection of events represented using the NanoEventsArray class.

    Returns:
    --------
        An Awkward Array mask containing the selected "good" muons that satisfy the specified criteria.
    """
    return (
        (events.Muon.pt >= 35)
        & (np.abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId)
        & (
            events.Muon.pfRelIso04_all < 0.15
            if hasattr(events.Muon, "pfRelIso04_all")
            else events.Muon.pfRelIso03_all < 0.15
        )
    )


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