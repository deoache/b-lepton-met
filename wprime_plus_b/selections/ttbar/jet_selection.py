import json
import numpy as np
import awkward as ak
import importlib.resources
from coffea.nanoevents.methods.nanoaod import JetArray


def select_good_bjets(
    jets: JetArray, year: str = "2017", working_point: str = "M"
) -> ak.highlevel.Array:
    """
    Selects and filters 'good' b-jets from a collection of jets based on specified criteria

    Parameters:
    -----------
    jets:
        A collection of jets
    year: {'2016', '2017', '2018'}
        Year for which the data is being analyzed. Default is '2017'.
    working_point: {'L', 'M', 'T'}
        Working point for b-tagging. Default is 'M'.

    Returns:
    --------
        An Awkward Array mask containing the selected "good" b-jets that satisfy the specified criteria.
    """
    # open and load btagDeepFlavB working point
    with importlib.resources.open_text("wprime_plus_b.data", "btagWPs.json") as file:
        btag_threshold = json.load(file)["deepJet"][year][working_point]
        
    return (
        (jets.pt >= 20)
        & (np.abs(jets.eta) < 2.4)
        & (jets.jetId == 6)
        & (jets.puId == 7)
        & (jets.btagDeepFlavB > btag_threshold)
    )
