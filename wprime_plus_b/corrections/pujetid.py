import correctionlib
import numpy as np
import awkward as ak
from typing import Type
from coffea.analysis_tools import Weights
from wprime_plus_b.corrections.utils import get_pog_json
from .utils import unflat_sf

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
    # flat jets array since correction function works only on flat arrays
    j, n = ak.flatten(jets), ak.num(jets)
    
    # get jet transverse momentum and pseudorapidity
    jets_pt = np.clip(j.pt, 12.5, 50)
    jets_eta = np.clip(j.eta, -5.8, 5.8)

    # define correction set
    cset = correctionlib.CorrectionSet.from_file(
        get_pog_json("pujetid", year + year_mod)
    )
    # get nominal scale factors
    nominal_sf = unflat_sf(
        cset["PUJetID_eff"].evaluate(jets_eta, jets_pt, "nom", working_point), n
    )

    if variation == "nominal":
        # get 'up' and 'down' variations
        up_sf = unflat_sf(
            cset["PUJetID_eff"].evaluate(jets_eta, jets_pt, "up", working_point), n
        )
        down_sf = unflat_sf(
            cset["PUJetID_eff"].evaluate(jets_eta, jets_pt, "down", working_point), n
        )

        # add nominal, up and down scale factors to weights container
        weights.add(
            name="pujetid",
            weight=nominal_sf,
            weightUp=up_sf,
            weightDown=down_sf,
        )
    else:
        # add nominal scale factors to weights container
        weights.add(
            name="pujetid",
            weight=nominal_sf
        )