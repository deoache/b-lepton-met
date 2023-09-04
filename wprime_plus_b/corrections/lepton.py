import correctionlib
import numpy as np
import awkward as ak
import importlib.resources
from typing import Type
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
    weights:
        Weights object from coffea.analysis_tools
    year:
        Year of the dataset {'2016', '2017', '2018'}
    year_mod:
        Year modifier {'', 'APV'}
    tag:
        label to include in the weight name
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
        tag: str = "electron",
        variation: str = "nominal"
    ) -> None:
        self.variation = variation
        
        # electron array
        self.electrons = electrons

        # electron transverse momentum and pseudorapidity
        self.electron_pt = np.array(ak.fill_none(self.electrons.pt, 0.0))
        self.electron_eta = np.array(ak.fill_none(self.electrons.eta, 0.0))

        # weights container
        self.weights = weights

        # define correction set
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron", year=year + year_mod)
        )
        self.year = year
        self.year_mod = year_mod
        self.pog_year = pog_years[year + year_mod]

        self.tag = tag

    def add_id_weight(self, id_working_point: str) -> None:
        """
        add electron identification scale factors to weights container

        Parameters:
        -----------
            id_working_point:
                Working point {'Loose', 'Medium', 'Tight', 'wp80iso', 'wp80noiso', 'wp90iso', 'wp90noiso'}
        """
        # electron pseudorapidity range: (-inf, inf)
        electron_eta = self.electron_eta

        # electron pt range: [10, inf)
        electron_pt = np.clip(
            self.electron_pt.copy(), 10.0, 499.999
        )  # potential problems with pt > 500 GeV

        # remove '_UL' from year
        year = self.pog_year.replace("_UL", "")

        # get scale factors
        values = {}
        values["nominal"] = self.cset["UL-Electron-ID-SF"].evaluate(
            year, "sf", id_working_point, electron_eta, electron_pt
        )
        if self.variation == "nominal":
            values["up"] = self.cset["UL-Electron-ID-SF"].evaluate(
                year, "sfup", id_working_point, electron_eta, electron_pt
            )
            values["down"] = self.cset["UL-Electron-ID-SF"].evaluate(
                year, "sfdown", id_working_point, electron_eta, electron_pt
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"],
            )
        else:
            self.weights.add(
                name=f"{self.tag}_id",
                weight=values["nominal"],
            )

    def add_reco_weight(self) -> None:
        """add electron reconstruction scale factors to weights container"""
        # electron pseudorapidity range: (-inf, inf)
        electron_eta = self.electron_eta

        # electron pt range: (20, inf)
        electron_pt = np.clip(
            self.electron_pt.copy(), 20.1, 499.999
        )  # potential problems with pt > 500 GeV

        # remove _UL from year
        year = self.pog_year.replace("_UL", "")

        # get scale factors
        values = {}
        values["nominal"] = self.cset["UL-Electron-ID-SF"].evaluate(
            year, "sf", "RecoAbove20", electron_eta, electron_pt
        )
        if self.variation == "nominal":
            values["up"] = self.cset["UL-Electron-ID-SF"].evaluate(
                year, "sfup", "RecoAbove20", electron_eta, electron_pt
            )
            values["down"] = self.cset["UL-Electron-ID-SF"].evaluate(
                year, "sfdown", "RecoAbove20", electron_eta, electron_pt
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_reco",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"],
            )
        else:
            self.weights.add(
                name=f"{self.tag}_reco",
                weight=values["nominal"],
            )

    def add_trigger_weight(self) -> None:
        """add electron trigger scale factors to weights container"""
        # corrections still not provided by POG
        with importlib.resources.path(
            "wprime_plus_b.data", f"electron_trigger_{self.pog_year}.json"
        ) as filename:
            # correction set
            cset = correctionlib.CorrectionSet.from_file(str(filename))

            # electron pt range: (10, 500)
            electron_pt = np.clip(self.electron_pt.copy(), 10, 499.999)

            # electron pseudorapidity range: (-2.5, 2.5)
            electron_eta = np.clip(self.electron_eta.copy(), -2.499, 2.499)

            # get scale factors (only nominal)
            values = {}
            values["nominal"] = cset["UL-Electron-Trigger-SF"].evaluate(
                electron_eta, electron_pt
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_trigger",
                weight=values["nominal"],
            )


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
class MuonCorrector:
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
    tag:
        label to include in the weight name
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
        tag: str = "muon",
        variation: str = "nominal",
        id_wp: str = "tight",
        iso_wp: str = "tight",
    ) -> None:
        self.variation = variation
        self.id_wp = id_wp
        self.iso_wp = iso_wp
        # muon array
        self.muons = muons

        # muon transverse momentum and pseudorapidity
        self.muon_pt = np.array(ak.fill_none(self.muons.pt, 0.0))
        self.muon_eta = np.array(ak.fill_none(self.muons.eta, 0.0))

        # weights container
        self.weights = weights

        # define correction set
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="muon", year=year + year_mod)
        )

        self.year = year
        self.year_mod = year_mod
        self.pog_year = pog_years[year + year_mod]

        self.tag = tag

    def add_id_weight(self) -> None:
        """
        add muon ID scale factors to weights container
        """
        self.add_weight(sf_type="id")

    def add_iso_weight(self) -> None:
        """
        add muon Iso (LooseRelIso with mediumID) scale factors to weights container
        """
        self.add_weight(sf_type="iso")

    def add_triggeriso_weight(self) -> None:
        assert self.id_wp == "tight" and self.iso_wp == "tight", "there's only available muon trigger SF for 'tight' ID and Iso"
        """add muon Trigger Iso (IsoMu24 or IsoMu27) scale factors"""
        # muon absolute pseudorapidity range: [0, 2.4)
        muon_eta = np.clip(np.abs(self.muon_eta).copy(), 0.0, 2.399)

        # muon pt range: [29, 200)
        muon_pt = np.clip(self.muon_pt.copy(), 29.0, 199.999)

        # scale factors keys
        sfs_keys = {
            "2016": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2017": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "2018": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        }
        # get scale factors
        values = {}
        values["nominal"] = self.cset[sfs_keys[self.year]].evaluate(
            self.pog_year, muon_eta, muon_pt, "sf"
        )
        if self.variation == "nominal":
            values["up"] = self.cset[sfs_keys[self.year]].evaluate(
                self.pog_year, muon_eta, muon_pt, "systup"
            )
            values["down"] = self.cset[sfs_keys[self.year]].evaluate(
                self.pog_year, muon_eta, muon_pt, "systdown"
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_triggeriso",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"],
            )
        else:
            self.weights.add(
                name=f"{self.tag}_triggeriso",
                weight=values["nominal"],
            )

    def add_weight(self, sf_type: str) -> None:
        """
        add muon ID (TightID) or Iso (LooseRelIso with mediumID) scale factors

        Parameters:
        -----------
            sf_type:
                Type of scale factor {'id', 'iso'}
        """
        if self.iso_wp == "tight":
            assert self.id_wp != "loose", "there's no available SFs"
            
        # muon absolute pseudorapidity range: [0, 2.4)
        muon_eta = np.clip(np.abs(self.muon_eta).copy(), 0.0, 2.399)

        # muon pt range: [15, 120)
        muon_pt = np.clip(self.muon_pt.copy(), 15.0, 119.999)

        # 'id' and 'iso' scale factors keys
        id_corrections = {
            "loose": "NUM_LooseID_DEN_TrackerMuons",
            "medium": "NUM_MediumID_DEN_TrackerMuons",
            "tight": "NUM_TightID_DEN_TrackerMuons"
        }
        if self.iso_wp == "loose":
            if self.id_wp == "loose":
                iso_correction = "NUM_LooseRelIso_DEN_LooseID"
            elif self.id_wp == "medium":
                iso_correction = "NUM_LooseRelIso_DEN_MediumID"
            elif self.id_wp == "tight": 
                iso_correction = "NUM_LooseRelIso_DEN_TightIDandIPCut"
        if self.iso_wp == "tight":
            if self.id_wp == "medium":
                iso_correction = "NUM_TightRelIso_DEN_MediumID"
            elif self.id_wp == "tight": 
                iso_correction = "NUM_TightRelIso_DEN_TightIDandIPCut"
            
        sfs_keys = {
            "id": id_corrections[self.id_wp],
            "iso": iso_correction
        }
        # get scale factors
        values = {}
        values["nominal"] = self.cset[sfs_keys[sf_type]].evaluate(
            self.pog_year, muon_eta, muon_pt, "sf"
        )
        if self.variation == "nominal":
            values["up"] = self.cset[sfs_keys[sf_type]].evaluate(
                self.pog_year, muon_eta, muon_pt, "systup"
            )
            values["down"] = self.cset[sfs_keys[sf_type]].evaluate(
                self.pog_year, muon_eta, muon_pt, "systdown"
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"{self.tag}_{sf_type}",
                weight=values["nominal"],
                weightUp=values["up"],
                weightDown=values["down"],
            )
        else:
            self.weights.add(
                name=f"{self.tag}_{sf_type}",
                weight=values["nominal"],
            )