import json
import copy
import correctionlib
import numpy as np
import awkward as ak
import importlib.resources
from typing import Type
from pathlib import Path
from .utils import unflat_sf
from coffea.analysis_tools import Weights
from wprime_plus_b.corrections.utils import pog_years, get_pog_json


# ----------------------------------
# lepton scale factors
# -----------------------------------
#
# Electron
#    - ID: wp80noiso?
#    - Recon: RecoAbove20?
#    - Trigger: ?
#
# working points: (Loose, Medium, RecoAbove20, RecoBelow20, Tight, Veto, wp80iso, wp80noiso, wp90iso, wp90noiso)
#
class ElectronCorrector:
    """
    Electron corrector class

    Parameters:
    -----------
    electrons:
        electron collection
    hlt:
        high level trigger branch
    weights:
        Weights object from coffea.analysis_tools
    year:
        Year of the dataset {'2016', '2017', '2018'}
    year_mod:
        Year modifier {'', 'APV'}
    variation:
        if 'nominal' (default) add 'nominal', 'up' and 'down'
        variations to weights container. else, add only 'nominal' weights.
    """

    def __init__(
        self,
        electrons: ak.Array,
        weights: Type[Weights],
        year: str = "2017",
        year_mod: str = "",
        variation: str = "nominal",
    ) -> None:
        self.electrons = electrons
        self.variation = variation
        self.nevents = len(electrons)

        # flat electrons array
        self.e, self.n = ak.flatten(electrons), ak.num(electrons)

        # weights container
        self.weights = weights

        # define correction set
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron", year=year + year_mod)
        )
        self.year = year
        self.year_mod = year_mod  # 2018
        self.pog_year = pog_years[year + year_mod]

    def add_trigger_weight(self, trigger_mask):
        """trigger weights are computed from only leading Electrons"""
        # get leading electrons
        le = ak.firsts(self.electrons)
        
        # get 'in-limits' electrons
        electron_pt_mask = (le.pt > 10.0) & (le.pt < 499.999)
        electron_eta_mask = np.abs(le.eta) < 2.4
        in_electron_mask = electron_pt_mask & electron_eta_mask & trigger_mask
        in_electrons = le.mask[in_electron_mask]

        # get electrons transverse momentum and pseudorapidity (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 10.0)
        electron_eta = ak.fill_none(in_electrons.eta, 0.0)

        # get eletron trigger correction
        cset = correctionlib.CorrectionSet.from_file(
            f"wprime_plus_b/data/correction_electron_trigger_{self.year + self.year_mod}.json.gz"
        )
        sf = cset["trigger_eff"].evaluate(electron_pt, electron_eta)
        nominal_sf = ak.where(in_electron_mask, sf, 1.0)

        # replace zero-value SF for 1
        # zero_sf_mask = nominal_sf == 0.0
        # nominal_sf = ak.where(zero_sf_mask, 1., nominal_sf)

        self.weights.add(
            name=f"ele_trigger",
            weight=nominal_sf,
        )

    def add_id_weight(self, id_working_point: str) -> None:
        """
        add electron identification scale factors to weights container

        Parameters:
        -----------
            id_working_point:
                Working point {'Loose', 'Medium', 'Tight', 'wp80iso', 'wp80noiso', 'wp90iso', 'wp90noiso'}
        """
        id_wps = {
            # mva ID working points https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
            "wp80iso": self.e.mvaFall17V2Iso_WP80,
            "wp90iso": self.e.mvaFall17V2Iso_WP90,
            "wp80noiso": self.e.mvaFall17V2noIso_WP80,
            "wp90noiso": self.e.mvaFall17V2noIso_WP90,
            # cutbased ID working points https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            "loose": self.e.cutBased == 2,
            "medium": self.e.cutBased == 3,
            "tight": self.e.cutBased == 4,
        }
        # get 'in-limits' electrons
        electron_pt_mask = (self.e.pt > 10.0) & (
            self.e.pt < 499.999
        )  # potential problems with pt > 500 GeV
        electron_id_mask = id_wps[id_working_point]
        in_electron_mask = electron_pt_mask & electron_id_mask
        in_electrons = self.e.mask[in_electron_mask]
        
        # get electrons transverse momentum and pseudorapidity (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 10.0)
        electron_eta = ak.fill_none(in_electrons.eta, 0.0)
        
        # remove '_UL' from year
        year = self.pog_year.replace("_UL", "")
        
        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset["UL-Electron-ID-SF"].evaluate(
                year, "sf", id_working_point, electron_eta, electron_pt
            ),
            in_electron_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset["UL-Electron-ID-SF"].evaluate(
                    year, "sfup", id_working_point, electron_eta, electron_pt
                ),
                in_electron_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset["UL-Electron-ID-SF"].evaluate(
                    year, "sfdown", id_working_point, electron_eta, electron_pt
                ),
                in_electron_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"electron_id",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"electron_id",
                weight=nominal_sf,
            )

    def add_reco_weight(self) -> None:
        """add electron reconstruction scale factors to weights container"""
        # get 'in-limits' electrons
        electron_pt_mask = (self.e.pt > 20.1) & (
            self.e.pt < 499.999
        )  # potential problems with pt > 500 GeV
        in_electron_mask = electron_pt_mask
        in_electrons = self.e.mask[in_electron_mask]
        
        # get electrons transverse momentum and pseudorapidity (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 20.1)
        electron_eta = ak.fill_none(in_electrons.eta, 0.0)
        
        # remove _UL from year
        year = self.pog_year.replace("_UL", "")
        
        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset["UL-Electron-ID-SF"].evaluate(
                year, "sf", "RecoAbove20", electron_eta, electron_pt
            ),
            in_electron_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset["UL-Electron-ID-SF"].evaluate(
                    year, "sfup", "RecoAbove20", electron_eta, electron_pt
                ),
                in_electron_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset["UL-Electron-ID-SF"].evaluate(
                    year, "sfdown", "RecoAbove20", electron_eta, electron_pt
                ),
                in_electron_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"electron_reco",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"electron_reco",
                weight=nominal_sf,
            )

