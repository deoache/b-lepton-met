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
        self.year_mod = year_mod # 2018
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
        nominal_sf = ak.where(in_electron_mask, sf, 1.)
        
        # replace zero-value SF for 1
        #zero_sf_mask = nominal_sf == 0.0
        #nominal_sf = ak.where(zero_sf_mask, 1., nominal_sf)
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
def get_id_wps(muons):
    return {
        # cutbased ID working points
        "loose": muons.looseId,
        "medium": muons.mediumId,
        "tight": muons.tightId,
    }
    
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
        id_wp: str = "tight",
        iso_wp: str = "tight",
    ) -> None:
        self.muons = muons
        self.variation = variation
        self.id_wp = id_wp
        self.iso_wp = iso_wp
        # muon array
        self.muons = muons
        # flat muon array
        self.m, self.n = ak.flatten(muons), ak.num(muons)
        # weights container
        self.weights = weights
        # define correction set
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="muon", year=year + year_mod)
        )
        self.year = year
        self.year_mod = year_mod
        self.pog_year = pog_years[year + year_mod]

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

    def add_triggeriso_weight(self, trigger_mask) -> None:
        """
        add muon Trigger Iso (IsoMu24 or IsoMu27) scale factors.
        trigger weights are computed from only leading Muons
        """
        assert (
            self.id_wp == "tight" and self.iso_wp == "tight"
        ), "there's only available muon trigger SF for 'tight' ID and Iso"
        # get leading muons
        lm = ak.firsts(self.muons)
        
        # get 'in-limits' muons
        muon_pt_mask = (lm.pt > 29.0) & (lm.pt < 199.999)
        muon_eta_mask = np.abs(lm.eta) < 2.399
        muon_id_mask = get_id_wps(lm)[self.id_wp]
        muon_iso_mask = get_iso_wps(lm)[self.iso_wp]
        in_muon_mask = muon_pt_mask & muon_eta_mask & muon_id_mask & muon_iso_mask & trigger_mask
        in_muons = lm.mask[in_muon_mask]
        
        # get muons transverse momentum and abs pseudorapidity (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 29.0)
        muon_eta = np.abs(ak.fill_none(in_muons.eta, 0.0))
        
        # scale factors keys
        sfs_keys = {
            "2016": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2017": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "2018": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        }
        # get nominal scale factors        
        sf = self.cset[sfs_keys[self.year]].evaluate(
            self.pog_year, muon_eta, muon_pt, "sf"
        )
        nominal_sf = ak.where(in_muon_mask, sf, 1.)
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = self.cset[sfs_keys[self.year]].evaluate(
                self.pog_year, muon_eta, muon_pt, "systup"
            )
            up_sf = ak.where(in_muon_mask, up_sf, 1.)
            
            down_sf = self.cset[sfs_keys[self.year]].evaluate(
                self.pog_year, muon_eta, muon_pt, "systdown"
            )
            down_sf = ak.where(in_muon_mask, down_sf, 1.)
            # add scale factors to weights container
            self.weights.add(
                name=f"muon_triggeriso",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"muon_triggeriso",
                weight=nominal_sf,
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
        assert self.iso_wp != "medium", "Only LooseRelIso and TightRelIso avaliable"

        # get 'in-limits' muons
        muon_pt_mask = (self.m.pt > 15.0) & (self.m.pt < 199.999)
        muon_eta_mask = np.abs(self.m.eta) < 2.399
        muon_id_mask = get_id_wps(self.m)[self.id_wp]
        muon_iso_mask = get_iso_wps(self.m)[self.iso_wp]
        in_muon_mask = muon_pt_mask & muon_eta_mask & muon_id_mask & muon_iso_mask
        in_muons = self.m.mask[in_muon_mask]
        # get muons transverse momentum and abs pseudorapidity (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 15.0)
        muon_eta = np.abs(ak.fill_none(in_muons.eta, 0.0))
        # 'id' and 'iso' scale factors keys
        id_corrections = {
            "loose": "NUM_LooseID_DEN_TrackerMuons",
            "medium": "NUM_MediumID_DEN_TrackerMuons",
            "tight": "NUM_TightID_DEN_TrackerMuons",
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
        sfs_keys = {"id": id_corrections[self.id_wp], "iso": iso_correction}
        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset[sfs_keys[sf_type]].evaluate(
                self.pog_year, muon_eta, muon_pt, "sf"
            ),
            in_muon_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset[sfs_keys[sf_type]].evaluate(
                    self.pog_year, muon_eta, muon_pt, "systup"
                ),
                in_muon_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset[sfs_keys[sf_type]].evaluate(
                    self.pog_year, muon_eta, muon_pt, "systdown"
                ),
                in_muon_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"muon_{sf_type}",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"muon_{sf_type}",
                weight=nominal_sf,
            )
            
            
# ----------------------------------
# Tau scale factors
# ----------------------------------

"""
TauID corrections

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
https://github.com/cms-tau-pog/TauIDSFs
https://github.com/uhh-cms/hh2bbtautau/blob/7666ed0426c87baa8d143ec26a216c3cadde513b/hbt/calibration/tau.py#L60

Good example:
https://github.com/schaefes/hh2bbtautau/blob/da6d47a7ddb2b1e7ffda06b8a96c6ddead2824b8/hbt/production/tau.py#L108

# genuine taus
* DeepTau2017v2p1VSjet (pt [-inf, inf); dm (0, 1, 2, 10, 11) ; genmatch (0, 1, 2, 3, 4, 5, 6 ); wp (Loose, Medium, Tight, VTight); wp_VSe (Tight, VVLoose );  syst ; flag (dm, pt))


# electrons faking taus
* DeepTau2017v2p1VSe = eta [0.0, 2.3); genmatch (0, 1); wp (Loose, Medium, Tight, VLoose, VTight, VVLoose, VVTight; syst (down, nom, up)

# muons faking taus
* DeepTau2017v2p1VSmu (eta [0.0, 2.3); genmatch (0, 2); wp (Loose, Medium, Tight, VLoose); syst (down, nom, up))

        
"""
class TauCorrector:
    def __init__(
        self,
        taus: ak.Array,
        weights: Type[Weights],
        year: str = "2017",
        year_mod: str = "",
        tau_vs_jet: str = "Tight",
        tau_vs_ele: str = "Tight",
        tau_vs_mu: str = "Tight",
        variation: str = "nominal",
    ) -> None:

        # flat taus array
        self.taus, self.n = ak.flatten(taus), ak.num(taus)

        # tau transverse momentum and pseudorapidity
        self.taus_pt = self.taus.pt
        self.taus_eta = self.taus.eta

        # tau genPartFlav and decayMode
        self.taus_genMatch = self.taus.genPartFlav
        self.taus_dm = self.taus.decayMode
        self.taus_wp_jet = self.taus.idDeepTau2017v2p1VSjet
        self.taus_wp_e = self.taus.idDeepTau2017v2p1VSe
        self.taus_wp_mu = self.taus.idDeepTau2017v2p1VSmu

        self.weights = weights
        self.year = year
        self.year_mod = year_mod

        # DeepTau working points
        self.tau_vs_jet = tau_vs_jet
        self.tau_vs_ele = tau_vs_ele
        self.tau_vs_mu = tau_vs_mu
        with importlib.resources.open_text("wprime_plus_b.data", "tau_wps.json") as file:
            taus_wp = json.load(file)
        self.tau_vs_jet_wp = taus_wp["DeepTau2017"]["deep_tau_jet"][tau_vs_jet]
        self.tau_vs_ele_wp = taus_wp["DeepTau2017"]["deep_tau_electron"][tau_vs_ele]
        self.tau_vs_mu_wp = taus_wp["DeepTau2017"]["deep_tau_muon"][tau_vs_mu]

        self.variation = variation

        # define correction set_id
        self.cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="tau", year=self.year + self.year_mod)
        )
        self.pog_year = pog_years[year + year_mod]
        """
        Check: https://github.com/cms-tau-pog/TauFW/blob/43bc39474b689d9712107d53a953b38c3cd9d43e/PicoProducer/python/analysis/ModuleETau.py#L270 
        """

    # e -> tau_h fake rate SFs for DeepTau2017v2p1VSe
    # eta = (0, 2.3]; genMatch = 0,1; wp = Loose, Medium, Tight, VLoose, VTight, VVLoose, VVTight; syst: down, nom, up
    def add_id_weight_DeepTau2017v2p1VSe(self):
        """
        Sf is called with:

        evaluate(eta (real),  genmatch (int) , wp (string), syst (string))

        """
        # tau pseudorapidity range: [0, 2.3)
        tau_eta_mask = (self.taus_eta >= 0) & (self.taus_eta < 2.3)
        # GenMatch = 0 "unmatched", 1 "electron";
        tau_genMatch_mask = (self.taus_genMatch == 1) | (self.taus_genMatch == 3)
        # Only taus passing the wp stablished
        tau_wp_mask = self.taus_wp_e > self.tau_vs_ele_wp
        in_tau_mask = tau_genMatch_mask & tau_wp_mask  #  & tau_eta_mask
        # get 'in-limits' taus
        in_limit_taus = self.taus.mask[in_tau_mask]
        # get pt and eta
        # fill Nones with some 'in-limit' value
        tau_eta = ak.fill_none(in_limit_taus.eta, 0)
        tau_genMatch = ak.fill_none(in_limit_taus.genPartFlav, 0.0)
        # wp
        wp = self.tau_vs_ele
        # syst
        syst = "nom"
        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset["DeepTau2017v2p1VSe"].evaluate(tau_eta, tau_genMatch, wp, "nom"),
            in_tau_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset["DeepTau2017v2p1VSe"].evaluate(
                    tau_eta, tau_genMatch, wp, "up"
                ),
                in_tau_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset["DeepTau2017v2p1VSe"].evaluate(
                    tau_eta, tau_genMatch, wp, "down"
                ),
                in_tau_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"e -> tau_h fake rate_{self.tau_vs_ele}",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"e -> tau_h fake rate_{self.tau_vs_ele}",
                weight=nominal_sf,
            )
        return nominal_sf

    # mu -> tau_h fake rate SFs for DeepTau2017v2p1VSmu
    # eta = (0, 2.3]; genMatch = 0,2; wp = Loose, Medium, Tight, VLoose ; syst: down, nom, up
    def add_id_weight_DeepTau2017v2p1VSmu(self):
        """
        Sf is called with:

        evaluate(eta (real),  genmatch (int) , wp (string), syst (string))

        """
        # tau pseudorapidity range: [0, 2.3)
        tau_eta_mask = (self.taus_eta >= 0) & (self.taus_eta < 2.3)
        # GenMatch = 0 "unmatched", 2 "muon";
        tau_genMatch_mask = (self.taus_genMatch == 2) | (self.taus_genMatch == 4)
        # Only taus passing the wp stablished
        tau_wp_mask = self.taus_wp_mu > self.tau_vs_mu_wp
        in_tau_mask = tau_genMatch_mask & tau_wp_mask  # & tau_eta_mask
        # get 'in-limits' taus
        in_limit_taus = self.taus.mask[in_tau_mask]
        # get pt and etaF
        # fill Nones with some 'in-limit' value
        tau_eta = ak.fill_none(in_limit_taus.eta, 0)
        tau_genMatch = ak.fill_none(in_limit_taus.genPartFlav, 0.0)
        # wp
        wp = self.tau_vs_mu
        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset["DeepTau2017v2p1VSmu"].evaluate(tau_eta, tau_genMatch, wp, "nom"),
            in_tau_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset["DeepTau2017v2p1VSmu"].evaluate(
                    tau_eta, tau_genMatch, wp, "up"
                ),
                in_tau_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset["DeepTau2017v2p1VSmu"].evaluate(
                    tau_eta, tau_genMatch, wp, "down"
                ),
                in_tau_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"mu -> tau_h fake rate_{self.tau_vs_mu}",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"mu -> tau_h fake rate_{self.tau_vs_mu}",
                weight=nominal_sf,
            )
        return nominal_sf

    # By default, use the pT-dependent SFs with the 'pt' flag
    # pt = (-inf, inf); dm = 0, 1, 2, 10, 11; genmatch = 0, 1, 2, 3, 4, 5, 6; wp = Loose, Medium, Tight, VTight; wp_VSe = Tight, VVLoose; syst = down, nom, up; flag = dm, pt
    def add_id_weight_DeepTau2017v2p1VSjet(self, flag: str = "pt"):
        """
        https://github.com/LEAF-HQ/LEAF/blob/d22cc55594a4b16d061c25dbf7ecdec04eedbc34/Analyzer/src/TauScaleFactorApplicatorJson.cc#L28

        Sf is called with:

        evaluate(pt (real),  dm (int), genmatch (int), wp (string), wp_VSe (string), syst (string), flag (string))

         - dm (decay mode): 0 (tau->pi); 1 (tau->rho->pi+pi0); 2 (tau->a1->pi+2pi0); 10 (tau->a1->3pi); 11 (tau->3pi+pi0)
         - getmatch: 0 or 6 = unmatched or jet, 1 or 3 = electron, 2 or 4 = muon, 5 = real tau
         - flag: We have worked in 'pt' = pT-dependent

        """
        # tau decayMode
        tau_dm_mask = (
            (self.taus_dm == 0)
            | (self.taus_dm == 1)
            | (self.taus_dm == 2)
            | (self.taus_dm == 10)
            | (self.taus_dm == 11)
        )
        # GenMatch = 0 or 6 = unmatched or jet, 1 or 3 = electron, 2 or 4 = muon, 5 = real tau
        tau_genMatch_mask = self.taus_genMatch == 5
        # Only taus passing the wp stablished
        tau_wp_mask = self.taus_wp_jet > self.tau_vs_jet_wp
        in_tau_mask = tau_dm_mask & tau_genMatch_mask & tau_wp_mask
        # get 'in-limits' taus
        in_limit_taus = self.taus.mask[in_tau_mask]
        # get pt and eta
        # fill Nones with some 'in-limit' value
        tau_pt = ak.fill_none(in_limit_taus.pt, 0)
        tau_dm = ak.fill_none(in_limit_taus.decayMode, 0)
        tau_genMatch = ak.fill_none(in_limit_taus.genPartFlav, 0.0)
        # wp
        wp = self.tau_vs_jet
        # wp_VSe
        wp_VSe = self.tau_vs_ele
        # get nominal scale factors
        nominal_sf = unflat_sf(
            self.cset["DeepTau2017v2p1VSjet"].evaluate(
                tau_pt, tau_dm, tau_genMatch, wp, wp_VSe, "default", flag
            ),
            in_tau_mask,
            self.n,
        )
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset["DeepTau2017v2p1VSjet"].evaluate(
                    tau_pt, tau_dm, tau_genMatch, wp, wp_VSe, "up", flag
                ),
                in_tau_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset["DeepTau2017v2p1VSjet"].evaluate(
                    tau_pt, tau_dm, tau_genMatch, wp, wp_VSe, "down", flag
                ),
                in_tau_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"jet -> tau_h fake rate_{self.tau_vs_jet}_{flag}",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"jet -> tau_h fake rate_{self.tau_vs_jet}_{flag}",
                weight=nominal_sf,
            )
        return nominal_sf

        # By default, use the pT-dependent SFs with the 'pt' flag

    # pt = [24.59953, inf); dm = -1, 0, 1, 10; trigtype = 'ditau', 'etau', 'mutau', 'ditauvbf; wp "DeepTauVSjet"= Loose, Medium, Tight, VLoose, VTight, VVLoose, VVTight, VVVLoose; corrtype =  eff_data, eff_mc, sf;  syst = down, nom, up

    def add_id_weight_diTauTrigger(
        self, mask_trigger, trigger: str = "ditau", info: str = "sf", dm: int = -1
    ):
        """
            Tau Trigger SFs and efficiencies for {0} ditau, etau, mutau or ditauvbf triggers. Ditauvbf
        trigger SF is only available for 2017 and 2018. To get the usual DM-specific SF's, specify the
        DM, otherwise set DM to -1 to get the inclusive SFs. Default corrections are set to SF's, if you
        require the input efficiencies, you can specify so in the corrtype input variable

        """
        # tau pt range: [24.59953, inf]
        tau_pt_mask = self.taus_pt >= 40
        # tau decayMode
        tau_dm_mask = (
            (self.taus_dm == -1)
            | (self.taus_dm == 0)
            | (self.taus_dm == 1)
            | (self.taus_dm == 10)
        )
        # Only taus passing the wp stablished
        tau_wp_mask = self.taus_wp_jet > self.tau_vs_jet_wp
        tau_mask = tau_pt_mask & tau_dm_mask & tau_wp_mask
        # get 'in-limits' taus
        in_limit_taus = self.taus.mask[tau_mask]
        # get pt and dm
        # fill Nones with some 'in-limit' value
        tau_pt = ak.fill_none(in_limit_taus.pt, 40.0)
        tau_dm = ak.fill_none(in_limit_taus.decayMode, -1)
        trigtype = trigger
        corrtype = info
        wp = self.tau_vs_jet
        # get nominal scale factors
        sf = unflat_sf(
            self.cset["tau_trigger"].evaluate(
                tau_pt, tau_dm, trigtype, wp, corrtype, "nom"
            ),
            tau_mask,
            self.n,
        )
        nominal_sf = np.where(mask_trigger, sf, 1.0)
        if self.variation == "nominal":
            # get 'up' and 'down' scale factors
            up_sf = unflat_sf(
                self.cset["tau_trigger"].evaluate(
                    tau_pt, tau_dm, trigtype, wp, corrtype, "up"
                ),
                tau_mask,
                self.n,
            )
            down_sf = unflat_sf(
                self.cset["tau_trigger"].evaluate(
                    tau_pt, tau_dm, trigtype, wp, corrtype, "down"
                ),
                tau_mask,
                self.n,
            )
            # add scale factors to weights container
            self.weights.add(
                name=f"trigger_{trigtype}",
                weight=nominal_sf,
                weightUp=up_sf,
                weightDown=down_sf,
            )
        else:
            self.weights.add(
                name=f"trigger_{trigtype}",
                weight=nominal_sf,
            )
        return nominal_sf