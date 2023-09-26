import correctionlib
import numpy as np
import awkward as ak
from typing import Type
from coffea.analysis_tools import Weights
from wprime_plus_b.corrections.utils import get_pog_json


def add_pujetid_weight(
    jets: ak.Array,
    weights: Type[Weights],
    year: str = "2017",
    year_mod: str = "",
    working_point: str = "T",
    variation: str = "nominal",
):
    """
    add jet pileup ID scale factor

    Parameters:
    -----------
        jets:
            Jet collection
        weights:
            Weights object from coffea.analysis_tools
        year:
            dataset year {'2016', '2017', '2018'}
        year_mod:
            year modifier {"", "APV"}
        working_point:
            pujetId working point {'L', 'M', 'T'}
        variation:
            if 'nominal' (default) add 'nominal', 'up' and 'down'
            variations to weights container. else, add only 'nominal' weights.
    """
    # define correction set
    cset = correctionlib.CorrectionSet.from_file(
        get_pog_json("pujetid", year + year_mod)
    )

    # jet transverse momentum and pseudorapidity
    jet_pt = np.clip(jets.pt, 12.5, 50)
    jet_eta = np.array(jets.eta)

    # get scale factors
    scale_factors = {}
    scale_factors["nominal"] = cset["PUJetID_eff"].evaluate(
        jet_eta, jet_pt, "nom", working_point
    )

    if variation == "nominal":
        scale_factors["up"] = cset["PUJetID_eff"].evaluate(
            jet_eta, jet_pt, "up", working_point
        )
        scale_factors["down"] = cset["PUJetID_eff"].evaluate(
            jet_eta, jet_pt, "down", working_point
        )
        # add scale factors to weights container
        weights.add(
            name="pujetid",
            weight=scale_factors["nominal"]
            if jet_pt.ndim == 1
            else np.prod(scale_factors["nominal"], axis=1),
            weightUp=scale_factors["up"]
            if jet_pt.ndim == 1
            else np.prod(scale_factors["up"], axis=1),
            weightDown=scale_factors["down"]
            if jet_pt.ndim == 1
            else np.prod(scale_factors["down"], axis=1),
        )
    else:
        weights.add(
            name="pujetid",
            weight=scale_factors["nominal"]
            if jet_pt.ndim == 1
            else np.prod(scale_factors["nominal"], axis=1),
        )
