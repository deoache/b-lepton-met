import json
import copy
import vector 
import pickle
import numpy as np
import awkward as ak
import importlib.resources
from coffea import processor
from coffea.analysis_tools import PackedSelection, Weights
from wprime_plus_b.processors.utils import histograms
from wprime_plus_b.corrections.jec import apply_jet_corrections
from wprime_plus_b.corrections.met import apply_met_phi_corrections
from wprime_plus_b.corrections.rochester import apply_rochester_corrections
from wprime_plus_b.corrections.tau_energy import apply_tau_energy_scale_corrections
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.l1prefiring import add_l1prefiring_weight
from wprime_plus_b.corrections.pujetid import add_pujetid_weight
from wprime_plus_b.corrections.btag import BTagCorrector
from wprime_plus_b.corrections.muon import MuonCorrector
from wprime_plus_b.corrections.muon_highpt import MuonHighPtCorrector
from wprime_plus_b.corrections.tau import TauCorrector
from wprime_plus_b.corrections.electron import ElectronCorrector
from wprime_plus_b.corrections.jetvetomaps import jetvetomaps_mask
from wprime_plus_b.selections.susy.electron_config import susy_electron_config
from wprime_plus_b.selections.susy.muon_config import susy_muon_config
from wprime_plus_b.selections.susy.muon_veto_config import susy_muon_veto_config
from wprime_plus_b.selections.susy.tau_config import susy_tau_config
from wprime_plus_b.selections.susy.bjet_config import susy_bjet_config
from wprime_plus_b.selections.susy.jet_config import susy_jet_config
from wprime_plus_b.selections.susy.electron_selection import select_good_electrons
from wprime_plus_b.selections.susy.muon_selection import select_good_muons
from wprime_plus_b.selections.susy.muon_veto_selection import select_good_veto_muons
from wprime_plus_b.selections.susy.tau_selection import select_good_taus
from wprime_plus_b.selections.susy.bjet_selection import select_good_bjets
from wprime_plus_b.selections.susy.jet_selection import select_good_jets
from wprime_plus_b.processors.utils.analysis_utils import (
    delta_r_mask,
    normalize,
    trigger_match,
)


class SusyAnalysis(processor.ProcessorABC):
    def __init__(
        self,
        year: str = "2017",
        output_type: str = "hist",
        overflow: str = "True",
    ):
        self.year = year
        self.output_type = output_type
        self.overflow = overflow
        # initialize dictionary of hists for control regions
        self.hist_dict = {
            "dimuon_kin": histograms.susy_dimuon_hist,
            "met_kin": histograms.susy_met_hist,
            "dijet_kin": histograms.susy_dijet_hist,
            "muon_kin": histograms.susy_muon_hist,
            "jet_kin": histograms.susy_jet_hist,
        }
    def process(self, events):
        # get dataset name
        dataset = events.metadata["dataset"]
        # get number of events before selection
        nevents = len(events)
        # check if sample is MC
        self.is_mc = hasattr(events, "genWeight")
        # create copies of histogram objects
        hist_dict = copy.deepcopy(self.hist_dict)
        # dictionary to store output data and metadata
        output = {}
        output["metadata"] = {}
        output["metadata"].update({"raw_initial_nevents": nevents})
        # define systematic variations shifts
        syst_variations = ["nominal"]
        for syst_var in syst_variations:
            # -------------------------------------------------------------
            # object corrections
            # -------------------------------------------------------------
            if self.is_mc:
                # apply JEC/JER corrections to jets (in data, the corrections are already applied)
                apply_jet_corrections(events, self.year)
                # apply energy corrections to taus (only to MC)
                apply_tau_energy_scale_corrections(
                    events=events, year=self.year, variation=syst_var
                )
            # apply rochester corretions to muons
            apply_rochester_corrections(
                events=events, is_mc=self.is_mc, year=self.year, variation=syst_var
            )
            # apply MET phi modulation corrections
            apply_met_phi_corrections(
                events=events,
                is_mc=self.is_mc,
                year=self.year,
            )
            # -------------------------------------------------------------
            # event SF/weights computation
            # -------------------------------------------------------------
            # get trigger mask
            with importlib.resources.path(
                "wprime_plus_b.data", "triggers.json"
            ) as path:
                with open(path, "r") as handle:
                    self._triggers = json.load(handle)[self.year]
            trigger_paths = self._triggers[susy_muon_config["muon_id_wp"]]
            trigger_mask = np.zeros(nevents, dtype="bool")
            for tp in trigger_paths:
                if tp in events.HLT.fields:
                    trigger_mask = trigger_mask | events.HLT[tp]
            # get DeltaR matched trigger objects mask
            trigger_match_mask = np.zeros(nevents, dtype="bool")
            for trigger_path in trigger_paths:
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
                add_l1prefiring_weight(events, weights_container, self.year, syst_var)
                # add pileup weigths
                add_pileup_weight(events, weights_container, self.year, syst_var)
                # add pujetid weigths
                add_pujetid_weight(
                    jets=events.Jet,
                    weights=weights_container,
                    year=self.year,
                    working_point=susy_bjet_config["jet_pileup_id"],
                    variation=syst_var,
                )
                # b-tagging corrector
                btag_corrector = BTagCorrector(
                    jets=events.Jet,
                    weights=weights_container,
                    sf_type="comb",
                    worging_point=susy_bjet_config["btag_working_point"],
                    tagger="deepJet",
                    year=self.year,
                    full_run=False,
                    variation=syst_var,
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
                electron_corrector.add_id_weight(
                    id_working_point=susy_electron_config["electron_id_wp"],
                )
                # add electron reco weights
                electron_corrector.add_reco_weight()

                # muon corrector
                if susy_muon_config["muon_id_wp"] == "highpt":
                    mu_corrector = MuonHighPtCorrector
                else:
                    mu_corrector = MuonCorrector
                muon_corrector = mu_corrector(
                    muons=events.Muon,
                    weights=weights_container,
                    year=self.year,
                    variation=syst_var,
                    id_wp=susy_muon_config["muon_id_wp"],
                    iso_wp=susy_muon_config["muon_iso_wp"],
                )
                # add muon RECO weights
                muon_corrector.add_reco_weight()
                # add muon ID weights
                muon_corrector.add_id_weight()
                # add muon iso weights
                muon_corrector.add_iso_weight()
                # add trigger weights
                muon_corrector.add_triggeriso_weight(
                    trigger_mask=trigger_mask,
                    trigger_match_mask=trigger_match_mask,
                )
                # add tau weights
                tau_corrector = TauCorrector(
                    taus=events.Tau,
                    weights=weights_container,
                    year=self.year,
                    tau_vs_jet=susy_tau_config["tau_vs_jet"],
                    tau_vs_ele=susy_tau_config["tau_vs_ele"],
                    tau_vs_mu=susy_tau_config["tau_vs_mu"],
                    variation=syst_var,
                )
                tau_corrector.add_id_weight_deeptauvse()
                tau_corrector.add_id_weight_deeptauvsmu()
                tau_corrector.add_id_weight_deeptauvsjet()
                
            if syst_var == "nominal":
                # save sum of weights before selections
                output["metadata"].update({"sumw": ak.sum(weights_container.weight())})
                # save weights statistics
                output["metadata"].update({"weight_statistics": {}})
                for weight, statistics in weights_container.weightStatistics.items():
                    output["metadata"]["weight_statistics"][weight] = statistics
            # -------------------------------------------------------------
            # object selection
            # -------------------------------------------------------------
            # select good electrons
            good_electrons = select_good_electrons(
                events=events,
                electron_pt_threshold=susy_electron_config["electron_pt_threshold"],
                electron_abs_eta=susy_electron_config["electron_abs_eta"],
                electron_id_wp=susy_electron_config["electron_id_wp"],
                electron_iso_wp=susy_electron_config["electron_iso_wp"],
            )
            electrons = events.Electron[good_electrons]

            # select good muons
            good_muons = select_good_muons(
                events=events,
                muon_pt_threshold=susy_muon_config["muon_pt_threshold"],
                abs_muon_eta=susy_muon_config["abs_muon_eta"],
                muon_id_wp=susy_muon_config["muon_id_wp"],
                muon_iso_wp=susy_muon_config["muon_iso_wp"],
            )
            muon_electron_dr = delta_r_mask(events.Muon, electrons, threshold=0.4)
            muons = events.Muon[good_muons & muon_electron_dr]
            
            
            good_veto_muons = select_good_veto_muons(
                events=events,
                muon_pt_threshold_min=susy_muon_veto_config["muon_pt_threshold_min"],
                muon_pt_threshold_max=susy_muon_veto_config["muon_pt_threshold_max"],
                abs_muon_eta=susy_muon_veto_config["abs_muon_eta"],
                muon_id_wp=susy_muon_veto_config["muon_id_wp"],
                muon_iso_wp=susy_muon_veto_config["muon_iso_wp"],
            )
            veto_muons = events.Muon[good_veto_muons & muon_electron_dr]
            
            # select good taus
            good_taus = select_good_taus(
                events=events,
                tau_pt_threshold=susy_tau_config["tau_pt_threshold"],
                tau_eta_threshold=susy_tau_config["tau_eta_threshold"],
                tau_dz_threshold=susy_tau_config["tau_dz_threshold"],
                tau_vs_jet=susy_tau_config["tau_vs_jet"],
                tau_vs_ele=susy_tau_config["tau_vs_ele"],
                tau_vs_mu=susy_tau_config["tau_vs_mu"],
                prong=susy_tau_config["prongs"],
            )
            good_taus = (
                (good_taus)
                & (delta_r_mask(events.Tau, electrons, threshold=0.3))
                & (delta_r_mask(events.Tau, muons, threshold=0.3))
            )
            taus = events.Tau[good_taus]

            # select good jets
            good_jets = select_good_jets(
                jets=events.Jet,
                year=self.year,
                jet_pt_threshold=susy_jet_config["jet_pt_threshold"],
                jet_abs_eta=susy_jet_config["jet_abs_eta"],
                jet_pileup_id=susy_jet_config["jet_pileup_id"],
            )
            good_jets = (
                good_jets
                & (delta_r_mask(events.Jet, electrons, threshold=0.4))
                & (delta_r_mask(events.Jet, muons, threshold=0.4))
                & (delta_r_mask(events.Jet, taus, threshold=0.4))
            )
            jets = events.Jet[good_jets]
            
            # select good bjets
            good_bjets = select_good_bjets(
                jets=events.Jet,
                year=self.year,
                btag_working_point=susy_bjet_config["btag_working_point"],
                jet_pt_threshold=susy_bjet_config["jet_pt_threshold"],
                jet_abs_eta=susy_bjet_config["jet_abs_eta"],
                jet_pileup_id=susy_bjet_config["jet_pileup_id"],
            )
            good_bjets = (
                good_bjets
                & (delta_r_mask(events.Jet, electrons, threshold=0.4))
                & (delta_r_mask(events.Jet, muons, threshold=0.4))
                & (delta_r_mask(events.Jet, taus, threshold=0.4))
            )
            # if self.year in ["2016APV", "2016", "2018"]:
            #    vetomask = jetvetomaps_mask(jets=events.Jet, year=self.year, mapname="jetvetomap")
            #    good_bjets = good_bjets & vetomask
            bjets = events.Jet[good_bjets]
            
            # -------------------------------------------------------------
            # composite object selection
            # -------------------------------------------------------------
            # add muons pT to MET to simulate a 0-lepton final state
            all_muons = ak.sum(muons, axis=1)
            muons2D = ak.zip(
                {
                    "pt": all_muons.pt,
                    "phi": all_muons.phi,
                },
                with_name="Momentum2D",
                behavior=vector.backends.awkward.behavior,
            )
            met2D = ak.zip(
                {
                    "pt": events.MET.pt,
                    "phi": events.MET.phi,
                },
                with_name="Momentum2D",
                behavior=vector.backends.awkward.behavior,
            )
            zl_state_met_pt = (met2D + muons2D).pt
            # create pair combinations with all muons
            dimuons = ak.combinations(muons, 2, fields=["mu1", "mu2"])
            # add dimuon 4-momentum field
            dimuons["p4"] = dimuons.mu1 + dimuons.mu2
            # impose some cuts on the dimuons
            dimuons = dimuons[
                ((dimuons.p4.mass > 60) & (dimuons.p4.mass < 120))
                & (dimuons.mu1.charge * dimuons.mu2.charge < 0)
            ]
            
            # create pair combinations with all jets (VBF selection)
            dijets = ak.combinations(jets, 2, fields=["j1", "j2"])
            # add dijet 4-momentum field
            dijets["p4"] = dijets.j1 + dijets.j2
            # impose some cuts on the dijets
            dijets = dijets[
                (np.abs(dijets.j1.eta - dijets.j2.eta) > 3.8)
                & (dijets.j1.eta * dijets.j2.eta < 0)
                & (dijets.p4.mass > 500)
            ]
            # get largest dijet mass
            largest_dijets_mass = ak.max(dijets.p4.mass, axis=1)
            
            # -------------------------------------------------------------
            # event selection
            # -------------------------------------------------------------
            # make a PackedSelection object to store selection masks
            self.selections = PackedSelection()
            # add luminosity calibration mask (only to data)
            with importlib.resources.path(
                "wprime_plus_b.data", "lumi_masks.pkl"
            ) as path:
                with open(path, "rb") as handle:
                    self._lumi_mask = pickle.load(handle)
            if not self.is_mc:
                lumi_mask = self._lumi_mask[self.year](
                    events.run, events.luminosityBlock
                )
            else:
                lumi_mask = np.ones(len(events), dtype="bool")
            self.selections.add("lumi", lumi_mask)
            # add lepton triggers masks
            self.selections.add("trigger", trigger_mask)
            # select events with at least one matched trigger object
            self.selections.add(
                "trigger_match", ak.sum(trigger_match_mask, axis=-1) > 0
            )
            # add MET filters mask
            with importlib.resources.path(
                "wprime_plus_b.data", "metfilters.json"
            ) as path:
                with open(path, "r") as handle:
                    self._metfilters = json.load(handle)[self.year]
            metfilters = np.ones(nevents, dtype="bool")
            metfilterkey = "mc" if self.is_mc else "data"
            for mf in self._metfilters[metfilterkey]:
                if mf in events.Flag.fields:
                    metfilters = metfilters & events.Flag[mf]
            self.selections.add("metfilters", metfilters)
            # select events with at least one good vertex
            self.selections.add("goodvertex", events.PV.npvsGood > 0)
            # cut on MET 
            self.selections.add("0lstate", zl_state_met_pt > 250)
            # stitching on DY inclusive samples
            dy_stitching = np.ones(nevents, dtype="bool")
            if dataset.startswith("DYJetsToLL_inclusive"):
                dy_stitching = events.LHE.HT < 70
            self.selections.add("dy_stitching", dy_stitching)
            # add number of leptons and jets
            self.selections.add("atleast_two_muons", ak.num(muons) > 1)
            self.selections.add("atleast_one_dimuon", ak.num(dimuons) > 0)
            self.selections.add("atleast_two_jets", ak.num(jets) > 1)
            self.selections.add("muon_veto", ak.num(veto_muons) == 0)
            self.selections.add("electron_veto", ak.num(electrons) == 0)
            self.selections.add("tau_veto", ak.num(taus) == 0)
            self.selections.add("bjet_veto", ak.num(bjets) == 0)

            if self.year == "2018":
                # hem-cleaning selection
                # https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
                # Due to the HEM issue in year 2018, we veto the events with jets and electrons in the
                # region -3 < eta <-1.3 and -1.57 < phi < -0.87 to remove fake MET

                hem_veto = ak.any(
                    (
                        (bjets.eta > -3.2)
                        & (bjets.eta < -1.3)
                        & (bjets.phi > -1.57)
                        & (bjets.phi < -0.87)
                    ),
                    -1,
                ) | ak.any(
                    (
                        (electrons.pt > 30)
                        & (electrons.eta > -3.2)
                        & (electrons.eta < -1.3)
                        & (electrons.phi > -1.57)
                        & (electrons.phi < -0.87)
                    ),
                    -1,
                )
                hem_cleaning = (
                    (
                        (events.run >= 319077) & (not self.is_mc)
                    )  # if data check if in Runs C or D
                    # else for MC randomly cut based on lumi fraction of C&D
                    | ((np.random.rand(len(events)) < 0.632) & self.is_mc)
                ) & (hem_veto)

                # self.selections.add("HEMCleaning", ~hem_cleaning)
                self.selections.add("HEMCleaning", np.ones(len(events), dtype="bool"))
            else:
                self.selections.add("HEMCleaning", np.ones(len(events), dtype="bool"))
            # define selection region

            region_selections = [
                "goodvertex",
                "lumi",
                "trigger",
                "trigger_match",
                "metfilters",
                "HEMCleaning",
                "dy_stitching",
                "0lstate",
                "atleast_two_muons",
                "atleast_one_dimuon",
                "atleast_two_jets",
                "muon_veto",
                "electron_veto",
                "tau_veto",
                "bjet_veto",
            ]
            self.selections.add("0lep", self.selections.all(*region_selections))
            region_selection = self.selections.all("0lep")

            # save cutflow
            if syst_var == "nominal":
                cut_names = region_selections
                output["metadata"].update({"cutflow": {}})
                selections = []
                for cut_name in cut_names:
                    selections.append(cut_name)
                    current_selection = self.selections.all(*selections)
                    output["metadata"]["cutflow"][cut_name] = ak.sum(
                        weights_container.weight()[current_selection]
                    )
            # -------------------------------------------------------------
            # event variables
            # -------------------------------------------------------------
            # check that there are events left after selection
            nevents_after = ak.sum(region_selection)
            if nevents_after > 0:
                # select region objects
                feature_map = {
                    "dimuon_mass": dimuons.p4.mass[region_selection],
                    "dimuon_pt": dimuons.p4.pt[region_selection],
                    "met": zl_state_met_pt[region_selection],
                    "dijet_mass": largest_dijets_mass[region_selection],
                    "mu1_pt": dimuons.mu1.pt[region_selection],
                    "mu2_pt": dimuons.mu1.pt[region_selection],
                    "jet_pt": jets.pt[region_selection],
                    "jet_eta": jets.eta[region_selection],
                    "lepton_pt": muons.pt[region_selection],
                    "lepton_eta": muons.eta[region_selection],
                }
                if syst_var == "nominal":
                    # save weighted events to metadata
                    output["metadata"].update(
                        {
                            "weighted_final_nevents": ak.sum(
                                weights_container.weight()[region_selection]
                            ),
                            "raw_final_nevents": nevents_after,
                        }
                    )
                # -------------------------------------------------------------
                # histogram filling
                # -------------------------------------------------------------
                if self.output_type == "hist":
                    # break up the histogram filling for event-wise variations and object-wise variations
                    # apply event-wise variations only for nominal
                    if self.is_mc and syst_var == "nominal":
                        # get event weight systematic variations for MC samples
                        variations = ["nominal"] + list(weights_container.variations)
                        for variation in variations:
                            if variation == "nominal":
                                region_weight = weights_container.weight()[
                                    region_selection
                                ]
                            else:
                                region_weight = weights_container.weight(
                                    modifier=variation
                                )[region_selection]
                            for kin in hist_dict:
                                hist_axes_names = [
                                    axis
                                    for axis in hist_dict[kin].axes.name
                                    if axis != "variation"
                                ]
                                hist_max_bin_edge = {
                                    axis: hist_dict[kin]
                                    .axes[axis]
                                    .edges[-1]
                                    - 0.1
                                    for axis in hist_axes_names
                                }
                                fill_args = {
                                    feature: (
                                        np.minimum(
                                            normalize(feature_map[feature]),
                                            hist_max_bin_edge[feature],
                                        )
                                        if self.overflow
                                        else normalize(feature_map[feature])
                                    )
                                    for feature in hist_dict[kin].axes.name
                                    if feature not in ["variation"]
                                }
                                hist_dict[kin].fill(
                                    **fill_args,
                                    variation=variation,
                                    weight=region_weight,
                                )
                    elif self.is_mc and syst_var != "nominal":
                        # object-wise variations
                        region_weight = weights_container.weight()[region_selection]
                        for kin in hist_dict:
                            # get filling arguments
                            hist_axes_names = [
                                axis
                                for axis in hist_dict[kin].axes.name
                                if axis != "variation"
                            ]
                            hist_max_bin_edge = {
                                axis: hist_dict[kin].axes[axis].edges[-1]
                                - 0.1
                                for axis in hist_axes_names
                            }
                            # get filling arguments
                            fill_args = {
                                feature: (
                                    np.minimum(
                                        normalize(feature_map[feature]),
                                        hist_max_bin_edge[feature],
                                    )
                                    if self.overflow
                                    else normalize(feature_map[feature])
                                )
                                for feature in hist_dict[kin].axes.name[
                                    :-1
                                ]
                                if feature not in ["variation"]
                            }
                            # fill histograms
                            hist_dict[kin].fill(
                                **fill_args,
                                variation=syst_var,
                                weight=region_weight,
                            )
                    elif not self.is_mc and syst_var == "nominal":
                        # object-wise variations
                        region_weight = weights_container.weight()[region_selection]
                        for kin in hist_dict:
                            # get filling arguments
                            hist_axes_names = [
                                axis
                                for axis in hist_dict[kin].axes.name
                                if axis != "variation"
                            ]
                            hist_max_bin_edge = {
                                axis: hist_dict[kin].axes[axis].edges[-1]
                                - 0.1
                                for axis in hist_axes_names
                            }
                            # get filling arguments
                            fill_args = {
                                feature: (
                                    np.minimum(
                                        normalize(feature_map[feature]),
                                        hist_max_bin_edge[feature],
                                    )
                                    if self.overflow
                                    else normalize(feature_map[feature])
                                )
                                for feature in hist_dict[kin].axes.name[
                                    :-1
                                ]
                                if feature not in ["variation"]
                            }
                            # fill histograms
                            hist_dict[kin].fill(
                                **fill_args,
                                variation=syst_var,
                                weight=region_weight,
                            )
        # define output dictionary accumulator
        output["histograms"] = hist_dict
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator