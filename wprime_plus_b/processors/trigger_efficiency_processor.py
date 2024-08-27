import json
import hist
import pickle
import numpy as np
import pandas as pd
import awkward as ak
from typing import List
from coffea import processor
import importlib.resources
from coffea.analysis_tools import Weights, PackedSelection
from wprime_plus_b.corrections.btag import BTagCorrector
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.l1prefiring import add_l1prefiring_weight
from wprime_plus_b.corrections.pujetid import add_pujetid_weight
from wprime_plus_b.corrections.electron import ElectronCorrector
from wprime_plus_b.corrections.muon import MuonCorrector
from wprime_plus_b.corrections.muon_highpt import MuonHighPtCorrector
from wprime_plus_b.corrections.tau import TauCorrector

from wprime_plus_b.processors.utils.analysis_utils import delta_r_mask, normalize, trigger_match
from wprime_plus_b.corrections.jec import apply_jet_corrections
from wprime_plus_b.corrections.met import apply_met_phi_corrections
from wprime_plus_b.corrections.rochester import apply_rochester_corrections
from wprime_plus_b.corrections.tau_energy import apply_tau_energy_scale_corrections

class TriggerEfficiencyProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year: str = "2017",
        lepton_flavor: str = "ele",
        output_type: str = "hist",
    ):
        self.year = year
        self.lepton_flavor = lepton_flavor

        # open triggers
        with open("wprime_plus_b/data/triggers.json", "r") as f:
            self._triggers = json.load(f)[self.year]
        # open btagDeepFlavB
        with open("wprime_plus_b/data/btagWPs.json", "r") as f:
            self._btagDeepFlavB = json.load(f)["deepJet"][self.year]["M"]
        # open met filters
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
        with open("wprime_plus_b/data/metfilters.json", "rb") as handle:
            self._metfilters = json.load(handle)[self.year]
        # open lumi masks
        with open("wprime_plus_b/data/lumi_masks.pkl", "rb") as handle:
            self._lumi_mask = pickle.load(handle)
        # output histograms
        self.make_output = lambda: {
            "electron_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(
                    [30, 60, 90, 120, 150, 180, 210, 240, 300, 500],
                    name="electron_pt",
                    label=r"electron $p_T$ [GeV]",
                ),
                hist.axis.Regular(
                    50, -2.4, 2.4, name="electron_eta", label="electron $\eta$"
                ),
                hist.storage.Weight(),
            ),
            "muon_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(
                    [30, 60, 90, 120, 150, 180, 210, 240, 300, 500],
                    name="muon_pt",
                    label=r"muon $p_T$ [GeV]",
                ),
                hist.axis.Regular(50, -2.4, 2.4, name="muon_eta", label="muon $\eta$"),
                hist.storage.Weight(),
            ),
            "jet_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(
                    [30, 60, 90, 120, 150, 180, 210, 240, 300, 500],
                    name="jet_pt",
                    label=r"bJet $p_T$ [GeV]",
                ),
                hist.axis.Regular(50, -2.4, 2.4, name="jet_eta", label="bJet $\eta$"),
                hist.storage.Weight(),
            ),
            "met_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(
                    [50, 75, 100, 125, 150, 175, 200, 300, 500],
                    name="met_pt",
                    label=r"$p_T^{miss}$ [GeV]",
                ),
                hist.axis.Regular(
                    50, -4.0, 4.0, name="met_phi", label=r"$\phi(p_T^{miss})$"
                ),
                hist.storage.Weight(),
            ),
            "electron_bjet_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Regular(
                    30, 0, 5, name="electron_bjet_dr", label="$\Delta R(e, bJet)$"
                ),
                hist.axis.Variable(
                    [40, 75, 100, 125, 150, 175, 200, 300, 500],
                    name="invariant_mass",
                    label=r"$m(e, bJet)$ [GeV]",
                ),
                hist.storage.Weight(),
            ),
            "muon_bjet_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Regular(
                    30, 0, 5, name="muon_bjet_dr", label="$\Delta R(\mu, bJet)$"
                ),
                hist.axis.Variable(
                    [40, 75, 100, 125, 150, 175, 200, 300, 500],
                    name="invariant_mass",
                    label=r"$m(\mu, bJet)$ [GeV]",
                ),
                hist.storage.Weight(),
            ),
            "lep_met_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(
                    [40, 75, 100, 125, 150, 175, 200, 300, 500, 800],
                    name="electron_met_transverse_mass",
                    label=r"$m_T(e, p_T^{miss})$ [GeV]",
                ),
                hist.axis.Variable(
                    [40, 75, 100, 125, 150, 175, 200, 300, 500, 800],
                    name="muon_met_transverse_mass",
                    label=r"$m_T(\mu, p_T^{miss})$ [GeV]",
                ),
                hist.storage.Weight(),
            ),
            "lep_bjet_met_kin": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(
                    [40, 75, 100, 125, 150, 175, 200, 300, 500, 800],
                    name="electron_total_transverse_mass",
                    label=r"$m_T^{tot}(e, bJet, p_T^{miss})$ [GeV]",
                ),
                hist.axis.Variable(
                    [40, 75, 100, 125, 150, 175, 200, 300, 500, 800],
                    name="muon_total_transverse_mass",
                    label=r"$m_T^{tot}(\mu, bJet, p_T^{miss})$ [GeV]",
                ),
                hist.storage.Weight(),
            ),
        }

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"]
        nevents = len(events)
        self.is_mc = hasattr(events, "genWeight")
        self.histograms = self.make_output()

        # dictionary to store output data and metadata
        output = {}

        if self.is_mc:
            apply_jet_corrections(events, self.year)

        # apply energy corrections to taus (only to MC)
        if self.is_mc:
            apply_tau_energy_scale_corrections(
                events=events, 
                year=self.year, 
                variation="nominal"
            )
        # apply rochester corretions to muons
        apply_rochester_corrections(
            events=events, 
            is_mc=self.is_mc, 
            year=self.year,
            variation="nominal"
        )
        # apply MET phi modulation corrections
        apply_met_phi_corrections(
            events=events,
            is_mc=self.is_mc,
            year=self.year,
        )
            
        # get trigger mask
        with importlib.resources.path(
            "wprime_plus_b.data", "triggers.json"
        ) as path:
            with open(path, "r") as handle:
                self._triggers = json.load(handle)[self.year]
    
        trigger_paths = {
            "2016APV": {
                "ele": ["Ele27_WPTight_Gsf", "Photon175"],
                "mu": ["IsoMu24"]
            },
            "2016": {
                "ele": ["Ele27_WPTight_Gsf", "Photon175"],
                "mu": ["IsoMu24"]
            },
            "2017": {
                "ele": ["Ele35_WPTight_Gsf", "Photon200"],
                "mu": ["IsoMu27"]
            },
            "2018": {
                "ele": ["Ele32_WPTight_Gsf"],
                "mu": ["IsoMu24"]
            }
        }
        def get_trigger(trigger_path):
            trigger_mask = np.zeros(nevents, dtype="bool")
            for tp in trigger_paths[self.year][self.lepton_flavor]:
                if tp in events.HLT.fields:
                    trigger_mask = trigger_mask | events.HLT[tp]
            return trigger_mask
    
        trigger_ele = get_trigger(trigger_paths[self.year]["ele"])
        trigger_mu = get_trigger(trigger_paths[self.year]["mu"])
        
        trigger_match_mask = np.zeros(nevents, dtype="bool")
        for trigger_path in trigger_paths[self.year]["mu"]:
            trig_match = trigger_match(
                leptons=events.Muon,
                trigobjs=events.TrigObj,
                trigger_path=trigger_path,
            )
            trigger_match_mask = trigger_match_mask | trig_match
                
        # set weights container
        weights_container = Weights(len(events), storeIndividual=True)
        if self.is_mc:
            # add gen weigths
            weights_container.add("genweight", events.genWeight)
            # add l1prefiring weigths
            add_l1prefiring_weight(events, weights_container, self.year, "nominal")
            # add pileup weigths
            add_pileup_weight(events, weights_container, self.year, "nominal")
            # add pujetid weigths
            add_pujetid_weight(
                jets=events.Jet,
                weights=weights_container,
                year=self.year,
                working_point="M",
                variation="nominal",
            )
            # b-tagging corrector
            btag_corrector = BTagCorrector(
                jets=events.Jet,
                weights=weights_container,
                sf_type="comb",
                worging_point="M",
                tagger="deepJet",
                year=self.year,
                full_run=False,
                variation="nominal",
            )
            # add b-tagging weights
            btag_corrector.add_btag_weights(flavor="bc")
            btag_corrector.add_btag_weights(flavor="light")
            # electron corrector
            electron_corrector = ElectronCorrector(
                electrons=events.Electron,
                weights=weights_container,
                year=self.year,
            )
            # add electron ID weights
            electron_corrector.add_id_weight("wp80iso" if self.lepton_flavor == "ele" else "wp90iso")
            # add electron reco weights
            electron_corrector.add_reco_weight()
            # add trigger weights
            if self.lepton_flavor == "ele":
                """
                electron_corrector.add_trigger_weight(
                    trigger_mask=trigger_mask["ele"],
                    trigger_match_mask=trigger_match_mask
                )
                """
                pass

            # muon corrector
            if self.lepton_flavor == "mu":
                mu_corrector = MuonHighPtCorrector
            else:
                mu_corrector = MuonCorrector
            muon_corrector = mu_corrector(
                muons=events.Muon,
                weights=weights_container,
                year=self.year,
                variation="nominal",
                id_wp="tight",
                iso_wp="tight",
            )
            # add muon ID weights
            muon_corrector.add_id_weight()
            # add muon iso weights
            muon_corrector.add_iso_weight()
            # add trigger weights
            if self.lepton_flavor == "ele":
                muon_corrector.add_triggeriso_weight(
                    trigger_mask=trigger_mu,
                    trigger_match_mask=trigger_match_mask,
                )

            # add tau weights
            tau_corrector = TauCorrector(
                taus=events.Tau,
                weights=weights_container,
                year=self.year,
                tau_vs_jet="Loose",
                tau_vs_ele="VVLoose",
                tau_vs_mu="Loose",
                variation="nominal",
            )
            tau_corrector.add_id_weight_deeptauvse()
            tau_corrector.add_id_weight_deeptauvsmu()
            tau_corrector.add_id_weight_deeptauvsjet()

        # save sum of weights before selections
        output["metadata"] = {"sumw": ak.sum(weights_container.weight())}
                
        # --------------------
        # object selection
        # --------------------
        # select electrons
        good_electrons = (
            (events.Electron.pt >= 30)
            & (np.abs(events.Electron.eta) < 2.4)
            & (
                (np.abs(events.Electron.eta) < 1.44)
                | (np.abs(events.Electron.eta) > 1.57)
            )
            & (
                events.Electron.mvaFall17V2Iso_WP80
                if self.lepton_flavor == "ele"
                else events.Electron.mvaFall17V2Iso_WP90
            )
        )
        electrons = events.Electron[good_electrons]  
        # select muons
        good_muons = (
            (events.Muon.pt >= 30)
            & (np.abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId)
            & (
                events.Muon.pfRelIso04_all < 0.15
                if hasattr(events.Muon, "pfRelIso04_all")
                else events.Muon.pfRelIso03_all < 0.15
            )
        )
        good_muons = (good_muons) & (
            delta_r_mask(events.Muon, electrons, threshold=0.4)
        )
        muons = events.Muon[good_muons]
        
        # correct and select taus
        prong_to_modes = {
            1: [0, 1, 2],
            2: [5, 6, 7],
            3: [10, 11],
            13: [0, 1, 2, 10, 11],
            12: [0, 1, 2, 5, 6, 7],
            23: [5, 6, 7, 10, 11],
        }
        tau_dm = events.Tau.decayMode
        decay_mode_mask = ak.zeros_like(tau_dm)
        for mode in prong_to_modes[13]:
            decay_mode_mask = np.logical_or(decay_mode_mask, tau_dm == mode)
            
        good_taus = (
            (events.Tau.pt > 20)
            & (np.abs(events.Tau.eta) < 2.3)
            & (np.abs(events.Tau.dz) < 0.2)
            & (events.Tau.idDeepTau2017v2p1VSjet > 32)
            & (events.Tau.idDeepTau2017v2p1VSe > 32)
            & (events.Tau.idDeepTau2017v2p1VSmu > 8)
            & (decay_mode_mask)
        )
        good_taus = (
            (good_taus)
            & (delta_r_mask(events.Tau, electrons, threshold=0.4))
            & (delta_r_mask(events.Tau, muons, threshold=0.4))
        )
        taus = events.Tau[good_taus]
        
        # b-jets
        # break up selection for low and high pT jets
        low_pt_jets_mask = (
            (events.Jet.pt > 20)
            & (events.Jet.pt < 50)
            & (np.abs(events.Jet.eta) < 2.4)
            & (events.Jet.jetId == 6)
            & (events.Jet.puId == 7)
            & (events.Jet.btagDeepFlavB > self._btagDeepFlavB)
        )
        high_pt_jets_mask = (
            (events.Jet.pt >= 50)
            & (np.abs(events.Jet.eta) < 2.4)
            & (events.Jet.jetId == 6)
            & (events.Jet.btagDeepFlavB > self._btagDeepFlavB)
        )
        good_bjets = ak.where(
            (events.Jet.pt > 20) & (events.Jet.pt < 50),
            low_pt_jets_mask,
            high_pt_jets_mask,
        )
        good_bjets = (
            good_bjets
            & (delta_r_mask(events.Jet, electrons, threshold=0.4))
            & (delta_r_mask(events.Jet, muons, threshold=0.4))
        )
        bjets = events.Jet[good_bjets]

        met = events.MET
        # ----------------
        # event selection
        # ----------------
        # make a PackedSelection object to store selection masks
        self.selections = PackedSelection()

        # luminosity
        if not self.is_mc:
            lumi_mask = self._lumi_mask[self.year](events.run, events.luminosityBlock)
        else:
            lumi_mask = np.ones(len(events), dtype="bool")
        self.selections.add("lumi", lumi_mask)

        # MET filters
        metfilters = np.ones(nevents, dtype="bool")
        metfilterkey = "mc" if self.is_mc else "data"
        for mf in self._metfilters[metfilterkey]:
            if mf in events.Flag.fields:
                metfilters = metfilters & events.Flag[mf]
        self.selections.add("metfilters", metfilters)

        # triggers 
        self.selections.add("trigger_ele", trigger_ele)
        self.selections.add("trigger_mu", trigger_mu)
        self.selections.add("atleastone_bjet", ak.num(bjets) >= 1)
        self.selections.add("one_electron", ak.num(electrons) == 1)
        self.selections.add("one_muon", ak.num(muons) == 1)
        self.selections.add("muon_veto", ak.num(muons) == 0)
        self.selections.add("electron_veto", ak.num(electrons) == 0)
        self.selections.add("tau_veto", ak.num(taus) == 0)

        # regions
        regions = {
            "ele": {
                "numerator": [
                    "trigger_ele",
                    "trigger_mu",
                    "lumi",
                    "metfilters",
                    "atleastone_bjet",
                    "one_muon",
                    "one_electron",
                    "tau_veto"
                ],
                "denominator": [
                    "trigger_mu",
                    "lumi",
                    "metfilters",
                    "atleastone_bjet",
                    "one_muon",
                    "one_electron",
                    "tau_veto"
                ],
            },
            "mu": {
                "numerator": [
                    "trigger_ele",
                    "trigger_mu",
                    "lumi",
                    "metfilters",
                    "atleastone_bjet",
                    "one_electron",
                    "one_muon",
                    "tau_veto"
                ],
                "denominator": [
                    "trigger_ele",
                    "lumi",
                    "metfilters",
                    "atleastone_bjet",
                    "one_electron",
                    "one_muon",
                    "tau_veto"
                ],
            },
        }
        # ---------------
        # event variables
        # ---------------
        # lepton-bjet delta R and invariant mass
        ele_bjet_dr = ak.firsts(bjets).delta_r(electrons)
        ele_bjet_mass = (electrons + ak.firsts(bjets)).mass
        mu_bjet_dr = ak.firsts(bjets).delta_r(muons)
        mu_bjet_mass = (muons + ak.firsts(bjets)).mass

        # lepton-MET transverse mass
        ele_met_tranverse_mass = np.sqrt(
            2.0
            * electrons.pt
            * met.pt
            * (ak.ones_like(met.pt) - np.cos(electrons.delta_phi(met)))
        )
        mu_met_transverse_mass = np.sqrt(
            2.0
            * muons.pt
            * met.pt
            * (ak.ones_like(met.pt) - np.cos(muons.delta_phi(met)))
        )
        # lepton-bJet-MET total transverse mass
        ele_total_transverse_mass = np.sqrt(
            (electrons.pt + ak.firsts(bjets).pt + met.pt) ** 2
            - (electrons + ak.firsts(bjets) + met).pt ** 2
        )
        mu_total_transverse_mass = np.sqrt(
            (muons.pt + ak.firsts(bjets).pt + met.pt) ** 2
            - (muons + ak.firsts(bjets) + met).pt ** 2
        )
        # filling histograms
        def fill(region: str):
            selections = regions[self.lepton_flavor][region]
            region_cut = self.selections.all(*selections)
            region_weight = weights_container.weight()[region_cut]

            self.histograms["jet_kin"].fill(
                region=region,
                jet_pt=normalize(ak.firsts(bjets).pt[region_cut]),
                jet_eta=normalize(ak.firsts(bjets).eta[region_cut]),
                weight=region_weight,
            )
            self.histograms["met_kin"].fill(
                region=region,
                met_pt=normalize(met.pt[region_cut]),
                met_phi=normalize(met.phi[region_cut]),
                weight=region_weight,
            )
            self.histograms["electron_kin"].fill(
                region=region,
                electron_pt=normalize(electrons.pt[region_cut]),
                electron_eta=normalize(electrons.eta[region_cut]),
                weight=region_weight,
            )
            self.histograms["muon_kin"].fill(
                region=region,
                muon_pt=normalize(muons.pt[region_cut]),
                muon_eta=normalize(muons.eta[region_cut]),
                weight=region_weight,
            )
            self.histograms["electron_bjet_kin"].fill(
                region=region,
                electron_bjet_dr=normalize(ele_bjet_dr[region_cut]),
                invariant_mass=normalize(ele_bjet_mass[region_cut]),
                weight=region_weight,
            )
            self.histograms["muon_bjet_kin"].fill(
                region=region,
                muon_bjet_dr=normalize(mu_bjet_dr[region_cut]),
                invariant_mass=normalize(mu_bjet_mass[region_cut]),
                weight=region_weight,
            )
            self.histograms["lep_met_kin"].fill(
                region=region,
                electron_met_transverse_mass=normalize(
                    ele_met_tranverse_mass[region_cut]
                ),
                muon_met_transverse_mass=normalize(mu_met_transverse_mass[region_cut]),
                weight=region_weight,
            )
            self.histograms["lep_bjet_met_kin"].fill(
                region=region,
                electron_total_transverse_mass=normalize(
                    ele_total_transverse_mass[region_cut]
                ),
                muon_total_transverse_mass=normalize(
                    mu_total_transverse_mass[region_cut]
                ),
                weight=region_weight,
            )
            """
            # cutflow
            cutflow_selections = []
            for selection in regions[self.lepton_flavor][region]:
                cutflow_selections.append(selection)
                cutflow_cut = self.selections.all(*cutflow_selections)

                cutflow_weight = weights_container.weight()
                print(region, selection, np.sum(cutflow_weight[cutflow_cut]))

            """
            
        for region in regions[self.lepton_flavor]:
            fill(region)
            
        output["metadata"].update({"raw_initial_nevents": nevents})
        output["histograms"] = self.histograms

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator