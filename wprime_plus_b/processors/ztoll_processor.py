import json
import copy
import pickle
import numpy as np
import awkward as ak
import importlib.resources
from coffea import processor
from wprime_plus_b.processors import utils
from coffea.analysis_tools import Weights, PackedSelection
from wprime_plus_b.corrections.jec import jet_corrections
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.lepton import ElectronCorrector, MuonCorrector
from wprime_plus_b.processors.utils.analysis_utils import delta_r_mask, normalize
from wprime_plus_b.selections.ztoll.jet_selection import select_good_bjets
from wprime_plus_b.selections.ztoll.config import (
    ztoll_electron_selection,
    ztoll_muon_selection,
    ztoll_jet_selection,
)
from wprime_plus_b.selections.ztoll.lepton_selection import (
    select_good_electrons,
    select_good_muons,
    select_good_taus,
)


class ZToLLProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year: str = "2017",
        yearmod: str = "",
        lepton_flavor: str = "ele",
        output_type="hist",
    ):
        self._year = year
        self._yearmod = yearmod
        self._lepton_flavor = lepton_flavor
        self._output_type = output_type

        # initialize output histogram
        self.hist_dict = {
            "dilepton_kin": utils.histograms.dilepton_mass_hist
        }
        # define dictionary to store analysis variables
        self.features = {}
        # initialize dictionary of arrays
        self.array_dict = {}
        
    def add_feature(self, name: str, var: ak.Array) -> None:
        """add a variable array to the out dictionary"""
        self.features = {**self.features, name: var}
    
    def process(self, events):
        # get dataset name
        dataset = events.metadata["dataset"]
        
        # check if sample is MC
        self.is_mc = hasattr(events, "genWeight")
        
        # get number of events before selection
        nevents = len(events)

        # create copies of histogram objects
        hist_dict = copy.deepcopy(self.hist_dict)
        # create copy of array dictionary
        array_dict = copy.deepcopy(self.array_dict)
        
        # ------------------
        # event preselection
        # ------------------
        # select good electrons
        good_electrons = select_good_electrons(
            events=events,
            electron_pt_threshold=ztoll_electron_selection["electron_pt_threshold"],
            electron_id_wp=ztoll_electron_selection["electron_id_wp"],
            electron_iso_wp=ztoll_electron_selection["electron_iso_wp"],
        )
        electrons = events.Electron[good_electrons]
        
        # select good muons
        good_muons = select_good_muons(
            events=events,
            muon_pt_threshold=ztoll_muon_selection["muon_pt_threshold"],
            muon_id_wp=ztoll_muon_selection["muon_id_wp"],
            muon_iso_wp=ztoll_muon_selection["muon_iso_wp"],
        )
        muons = events.Muon[good_muons]
        
        # define leptons collection
        leptons = electrons if self._lepton_flavor == "ele" else muons
        
        # select good taus
        good_taus = (
            select_good_taus(events)
            & (delta_r_mask(events.Tau, leptons, threshold=0.4))
        )
        taus = events.Tau[good_taus]
        
        # apply JEC/JER corrections to MC jets (propagate corrections to MET)
        # in data, the corrections are already applied
        if self.is_mc:
            jets, _ = jet_corrections(events, self._year + self._yearmod)
        else:
            jets, _ = events.Jet, events.MET
            
        # select good bjets
        good_bjets = (
            select_good_bjets(
                jets=jets,
                year=self._year + self._yearmod,
                jet_pt_threshold=ztoll_jet_selection["jet_pt_threshold"],
                btag_working_point=ztoll_jet_selection["btag_working_point"],   
            )
            & (delta_r_mask(jets, leptons, threshold=0.4))
            & (delta_r_mask(jets, taus, threshold=0.4))
        )
        bjets = jets[good_bjets]
        
        # ---------------
        # event selection
        # ---------------
        # make a PackedSelection object to manage selections
        self.selections = PackedSelection()
        
        # add luminosity calibration mask (only to data)
        with importlib.resources.path("wprime_plus_b.data", "lumi_masks.pkl") as path:
            with open(path, "rb") as handle:
                self._lumi_mask = pickle.load(handle)
        if not self.is_mc:
            lumi_mask = self._lumi_mask[self._year](events.run, events.luminosityBlock)
        else:
            lumi_mask = np.ones(len(events), dtype="bool")
        self.selections.add("lumi", lumi_mask)

        # add lepton triggers masks
        with importlib.resources.path("wprime_plus_b.data", "triggers.json") as path:
            with open(path, "r") as handle:
                self._triggers = json.load(handle)[self._year]
        trigger = {}
        for ch in ["ele", "mu"]:
            trigger[ch] = np.zeros(nevents, dtype="bool")
            for t in self._triggers[ch]:
                if t in events.HLT.fields:
                    trigger[ch] = trigger[ch] | events.HLT[t]
        self.selections.add("trigger_ele", trigger["ele"])
        self.selections.add("trigger_mu", trigger["mu"])
        
        # add MET filters mask
        with importlib.resources.path("wprime_plus_b.data", "metfilters.json") as path:
            with open(path, "r") as handle:
                self._metfilters = json.load(handle)[self._year]
        metfilters = np.ones(nevents, dtype="bool")
        metfilterkey = "mc" if self.is_mc else "data"
        for mf in self._metfilters[metfilterkey]:
            if mf in events.Flag.fields:
                metfilters = metfilters & events.Flag[mf]
        self.selections.add("metfilters", metfilters)
        
        
        # compute dilepton mass
        leading_lepton = ak.firsts(leptons)
        subleading_lepton = ak.pad_none(leptons, 2)[:, 1]
        dilepton_mass = (leading_lepton + subleading_lepton).mass
        
        # check that we have 2l events
        self.selections.add("two_leptons", ak.num(leptons) == 2)
        # check that dilepton system is neutral
        self.selections.add("neutral", leading_lepton.charge * subleading_lepton.charge < 0)
        # check that dilepton invariant mass is between 60 and 120 GeV
        self.selections.add("mass_range", (60 < dilepton_mass) & (dilepton_mass < 120))
        # veto bjets
        self.selections.add("bjet_veto", ak.num(bjets) == 0)
        # veto taus
        self.selections.add("tau_veto", ak.num(taus) == 0)
        
        # define selection regions for each lepton_channel
        regions = {
            "ele": [
                "lumi",
                "metfilters",
                "trigger_ele",
                "tau_veto",
                "bjet_veto",
                "two_leptons",
                "neutral",
                "mass_range",
            ],
            "mu": [
                "lumi",
                "metfilters",
                "trigger_mu",
                "tau_veto",
                "bjet_veto",
                "two_leptons",
                "neutral",
                "mass_range",
            ],
        }
        # ---------------
        # event variables
        # ---------------
        for lepton_flavor in regions:
            if lepton_flavor != self._lepton_flavor: continue
            
            region_selection = self.selections.all(*regions[lepton_flavor])
            # if there are no events left after selection cuts continue to the next .root file
            nevents_after = ak.sum(region_selection)
            if nevents_after == 0:
                continue
                
            # select region objects
            region_leptons = leptons[region_selection]
            region_leading_lepton = ak.firsts(region_leptons)
            region_subleading_lepton = ak.pad_none(region_leptons, 2)[:, 1]

            # add dilepton invariant mass to out
            region_dilepton_mass = (region_leading_lepton + region_subleading_lepton).mass
            self.add_feature("dilepton_mass", region_dilepton_mass)
            
            # features for corrections
            if self._output_type == "array":
                self.add_feature("L1PreFiringWeight", events.L1PreFiringWeight.Nom[region_selection])
                self.add_feature("leading_lepton_pt", region_leading_lepton.pt)
                self.add_feature("subleading_lepton_pt", region_subleading_lepton.pt)
                self.add_feature("leading_lepton_eta", region_leading_lepton.eta)
                self.add_feature("subleading_lepton_eta", region_subleading_lepton.eta)
                if self.is_mc:
                    self.add_feature("npu", events.Pileup.nPU[region_selection])
        
            # -------------
            # event weights
            # -------------
            # define weights container
            weights_container = Weights(
                    len(events[region_selection]), storeIndividual=True
                )
            if self.is_mc:
                # add gen weigths
                gen_weight = events.genWeight[region_selection]
                weights_container.add("genweight", gen_weight)
                
                # add L1prefiring weights
                if self._year in ("2016", "2017"):
                    weights_container.add(
                        "L1Prefiring",
                        weight=events.L1PreFiringWeight.Nom[region_selection],
                        weightUp=events.L1PreFiringWeight.Up[region_selection],
                        weightDown=events.L1PreFiringWeight.Dn[region_selection],
                    )
                # add pileup reweighting
                add_pileup_weight(
                    n_true_interactions=ak.to_numpy(
                        events.Pileup.nPU[region_selection]
                    ),
                    weights=weights_container,
                    year=self._year,
                    year_mod=self._yearmod,
                )

                # add lepton weights
                if self._lepton_flavor == "ele":
                    # leading electron corrector
                    leading_electron_corrector = ElectronCorrector(
                        electrons=region_leading_lepton,
                        weights=weights_container,
                        year=self._year,
                        year_mod=self._yearmod,
                        tag="leading_electron",
                    )
                    # add leading electron reco weights
                    leading_electron_corrector.add_reco_weight()
                    # add leading electron trigger weights
                    #leading_electron_corrector.add_trigger_weight()
                    
                    # subleading electron corrector
                    subleading_electron_corrector = ElectronCorrector(
                        electrons=region_subleading_lepton,
                        weights=weights_container,
                        year=self._year,
                        year_mod=self._yearmod,
                        tag="subleading_electron",
                    )
                    # add subleading electron reco weights
                    subleading_electron_corrector.add_reco_weight()
                    # add subleading electron trigger weights
                    #subleading_electron_corrector.add_trigger_weight()
                    
                    if ztoll_electron_selection["electron_id_wp"] == "heep":
                        # HEEP ID scale factor
                        # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#HEEP_ID_Scale_Factor_for_UL
                        leading_electron_heep_idsf = ak.where(
                            np.abs(region_leading_lepton.eta) < 1.4442,
                            ak.full_like(array=region_leading_lepton.pt, fill_value=0.979),
                            ak.full_like(array=region_leading_lepton.pt, fill_value=0.987)
                        )
                        subleading_electron_heep_idsf = ak.where(
                            np.abs(region_subleading_lepton.eta) < 1.4442,
                            ak.full_like(array=region_subleading_lepton.pt, fill_value=0.979),
                            ak.full_like(array=region_subleading_lepton.pt, fill_value=0.987)
                        )
                        weights_container.add("leading_electron_heep_id", ak.fill_none(leading_electron_heep_idsf, 1))
                        weights_container.add("subleading_electron_heep_id", ak.fill_none(subleading_electron_heep_idsf, 1))
                    else:
                        leading_electron_corrector.add_id_weight(ztoll_electron_selection["electron_id_wp"])
                        subleading_electron_corrector.add_id_weight(ztoll_electron_selection["electron_id_wp"])
                        
                if self._lepton_flavor == "mu":
                    # leading muon corrector
                    leading_muon_corrector = MuonCorrector(
                        muons=region_leading_lepton,
                        weights=weights_container,
                        year=self._year,
                        year_mod=self._yearmod,
                        tag="leading_muon",
                    )
                    # add muon ID weights
                    leading_muon_corrector.add_id_weight(working_point="tight")
                    # add muon iso weights
                    leading_muon_corrector.add_iso_weight(working_point="tight")
                    # add muon trigger weights
                    leading_muon_corrector.add_triggeriso_weight()
                    
                    # subleading muon corrector
                    subleading_muon_corrector = MuonCorrector(
                        muons=region_subleading_lepton,
                        weights=weights_container,
                        year=self._year,
                        year_mod=self._yearmod,
                        tag="leading_muon",
                    )
                    # add muon ID weights
                    subleading_muon_corrector.add_id_weight(working_point=ztoll_muon_selection["muon_id_wp"])
                    # add muon iso weights
                    subleading_muon_corrector.add_iso_weight(working_point=ztoll_muon_selection["muon_iso_wp"])
                    # add muon trigger weights
                    subleading_muon_corrector.add_triggeriso_weight()
                
            # get total weight from the weights container
            region_weights = weights_container.weight()
            
            # -----------------------------
            # fill histogram
            # -----------------------------
            if self._output_type == "hist":
                for kin in hist_dict:
                    fill_args = {
                        feature: utils.analysis_utils.normalize(self.features[feature])
                        for feature in hist_dict[kin].axes.name
                        if "dataset" not in feature
                    }
                    hist_dict[kin].fill(
                        **fill_args,
                        dataset=dataset,
                        weight=region_weights,
                    )
            else:
                self.add_feature("weights", region_weights)
                if self.is_mc:
                    self.add_feature("genweights", gen_weight)
                # select variables and put them in column accumulators
                array_dict = {
                    feature_name: processor.column_accumulator(
                        normalize(feature_array)
                    )
                    for feature_name, feature_array in self.features.items()
                }
        # define output
        output = {}
        if self._output_type == "hist":
            output["histograms"] = hist_dict
        else:
            output["arrays"] = array_dict
            
        output["metadata"] = {
            "events_before": nevents,
            "events_after": nevents_after,
        }
        
        # if dataset is montecarlo add sumw to output
        if self.is_mc:
            output["metadata"].update({"sumw": ak.sum(events.genWeight)})
        return output
    
    def postprocess(self, accumulator):
        return accumulator
