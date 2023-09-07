import json
import numpy as np
import awkward as ak
import importlib.resources


def select_good_bjets(jets, year: str, btag_working_point: str, jet_pt_threshold: int) -> ak.highlevel.Array:
    """
    Selects and filters 'good' b-jets from a collection of jets based on specified criteria

    Parameters:
    -----------
    jets:
        A collection of jets

    year: {'2016', '2017', '2018'}
        Year for which the data is being analyzed. Default is '2017'.

    btag_working_point: {'L', 'M', 'T'}
        Working point for b-tagging. Default is 'M'.

    Returns:
    --------
        An Awkward Array mask containing the selected "good" b-jets that satisfy the specified criteria.
    """
    # open and load btagDeepFlavB working point
    with importlib.resources.open_text("wprime_plus_b.data", "btagWPs.json") as file:
        btag_threshold = json.load(file)["deepJet"][year][btag_working_point]

    return (
        (jets.pt >= jet_pt_threshold)
        & (np.abs(jets.eta) < 2.4)
        & (jets.btagDeepFlavB > btag_threshold)
    )
    