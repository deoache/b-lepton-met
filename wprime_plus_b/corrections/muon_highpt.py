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


# Muon
#
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018
#
#    - ID: medium prompt ID NUM_MediumPromptID_DEN_TrackerMuon?
#    - Iso: LooseRelIso with mediumID (NUM_LooseRelIso_DEN_MediumID)?
#    - Trigger iso:
#          2016: for IsoMu24 (and IsoTkMu24?) NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight?
#          2017: for isoMu27 NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight?
#          2018: for IsoMu24 NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight?
#
def get_iso_wps(muons):
    return {
        "loose": (
            muons.pfRelIso04_all < 0.25
            if hasattr(muons, "pfRelIso04_all")
            else muons.pfRelIso03_all < 0.25
        ),
        "medium": (
            muons.pfRelIso04_all < 0.20
            if hasattr(muons, "pfRelIso04_all")
            else muons.pfRelIso03_all < 0.20
        ),
        "tight": (
            muons.pfRelIso04_all < 0.15
            if hasattr(muons, "pfRelIso04_all")
            else muons.pfRelIso03_all < 0.15
        ),
    }


class MuonHighPtCorrector:
    """
    Muon corrector class

    Parameters:
    -----------
    muons:
        muons collection
    weights:
        Weights object from coffea.analysis_tools
    year:
        Year of the dataset {'2016', '2017', '2018'}
    year_mod:
        Year modifier {'', 'APV'}
    variation:
        syst variation
    id_wp:
        ID working point {'loose', 'medium', 'tight'}
    iso_wp:
        Iso working point {'loose', 'medium', 'tight'}
    """

    def __init__(
        self,
        muons: ak.Array,
        weights: Type[Weights],
        year: str = "2017",
        year_mod: str = "",
        variation: str = "nominal",
        iso_wp: str = "tight",
    ) -> None:
        self.muons = muons
        self.variation = variation
        self.iso_wp = iso_wp
        self.id_wp = "highpt"
        
        # muon array
        self.muons = muons
        
        # flat muon array
        self.m, self.n = ak.flatten(muons), ak.num(muons)
        
        # weights container
        self.weights = weights
        
        # define correction set
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="muon_highpt", year=year + year_mod)
        )

        self.year = year
        self.year_mod = year_mod
        self.pog_year = pog_years[year + year_mod]

    def add_id_weight(self):
        """
        add muon ID scale factors to weights container
        """
        # get muons that pass the id wp, and within SF binning
        muon_pt_mask = self.m.pt > 50.0
        muon_eta_mask = np.abs(self.m.eta) < 2.39
        muon_id_mask = self.m.highPtId == 2
        in_muon_mask = muon_pt_mask & muon_eta_mask & muon_id_mask
        in_muons = self.m.mask[in_muon_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 50.0)
        muon_eta = np.abs(ak.fill_none(in_muons.eta, 0.0))

        # 'id' scale factors names
        id_corrections = {
            "2016APV": {},
            "2016": {},
            "2017": {
                "highpt": "NUM_HighPtID_DEN_GlobalMuonProbes"
            },
            "2018": {},
        }

        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset[id_corrections[self.year + self.year_mod][self.id_wp]].evaluate(
                muon_eta, muon_pt, "nominal"
            ),
            in_muon_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset[
                    id_corrections[self.year + self.year_mod][self.id_wp]
                ].evaluate(muon_eta, muon_pt, "systup"),
                in_muon_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset[
                    id_corrections[self.year + self.year_mod][self.id_wp]
                ].evaluate(muon_eta, muon_pt, "systdown"),
                in_muon_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"muon_id",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"muon_id",
                weight=nominal_sf,
            )

    def add_iso_weight(self):
        """
        add muon Iso (LooseRelIso with mediumID) scale factors to weights container
        """
        # get 'in-limits' muons
        muon_pt_mask = self.m.pt > 50.0
        muon_eta_mask = np.abs(self.m.eta) < 2.39
        muon_id_mask = self.m.highPtId == 2
        muon_iso_mask = get_iso_wps(self.m)[self.iso_wp]
        in_muon_mask = muon_pt_mask & muon_eta_mask & muon_id_mask & muon_iso_mask
        in_muons = self.m.mask[in_muon_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 50.0)
        muon_eta = np.abs(ak.fill_none(in_muons.eta, 0.0))

        iso_corrections = {
            "2016APV": {},
            "2016": {},
            "2017": {
                "loose": "NUM_probe_LooseRelTkIso_DEN_HighPtProbes",
                "medium": None,
                "tight": "NUM_probe_TightRelTkIso_DEN_HighPtProbes",
            },
            "2018": {},
        }

        correction_name = iso_corrections[self.year + self.year_mod][self.iso_wp]
        assert correction_name, "No Iso SF's available"

        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset[correction_name].evaluate(muon_eta, muon_pt, "nominal"),
            in_muon_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset[correction_name].evaluate(
                    muon_eta, muon_pt, "systup"
                ),
                in_muon_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset[correction_name].evaluate(
                    muon_eta, muon_pt, "systdown"
                ),
                in_muon_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"muon_iso",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"muon_iso",
                weight=nominal_sf,
            )

    def add_triggeriso_weight(self, trigger_mask) -> None:
        """
        add muon Trigger Iso (IsoMu24 or IsoMu27) scale factors.
        trigger weights are computed from only leading Muons
        """

        # get leading muons
        lm = ak.firsts(self.muons)

        # get 'in-limits' muons
        muon_pt_mask = lm.pt > 50.0
        muon_eta_mask = np.abs(lm.eta) < 2.399
        muon_id_mask = lm.highPtId == 2
        muon_iso_mask = get_iso_wps(lm)[self.iso_wp]
        in_muon_mask = (
            muon_pt_mask & muon_eta_mask & muon_id_mask & muon_iso_mask & trigger_mask
        )
        in_muons = lm.mask[in_muon_mask]

        # get muons transverse momentum and abs pseudorapidity (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 50.0)
        muon_eta = np.abs(ak.fill_none(in_muons.eta, 0.0))

        # scale factors keys
        sfs_keys = {
            "2016APV": None,
            "2016": None,
            "2017": "NUM_HLT_DEN_HighPtTightRelIsoProbes",
            "2018": None,
        }
        # get nominal scale factors
        sf = self.cset[sfs_keys[self.year + self.year_mod]].evaluate(
            muon_eta, muon_pt, "nominal"
        )
        nominal_sf = ak.where(in_muon_mask, sf, 1.0)
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = self.cset[sfs_keys[self.year]].evaluate(
                muon_eta, muon_pt, "systup"
            )
            up_sf = ak.where(in_muon_mask, up_sf, 1.0)

            down_sf = self.cset[sfs_keys[self.year]].evaluate(
                muon_eta, muon_pt, "systdown"
            )
            down_sf = ak.where(in_muon_mask, down_sf, 1.0)
            # add scale factors to weights container
            self.weights.add(
                name=f"muon_trigger",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"muon_trigger",
                weight=nominal_sf,
            )