import json
import copy
import pickle
import numpy as np
import awkward as ak
import importlib.resources
from coffea import processor
from coffea.analysis_tools import PackedSelection, Weights
from wprime_plus_b.processors.utils import histograms
from wprime_plus_b.corrections.jec import jet_corrections
from wprime_plus_b.corrections.met import met_phi_corrections
from wprime_plus_b.corrections.btag import BTagCorrector
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.l1prefiring import add_l1prefiring_weight
from wprime_plus_b.corrections.pujetid import add_pujetid_weight
from wprime_plus_b.corrections.tau_energy import tau_energy_scale, met_corrected_tes
from wprime_plus_b.processors.utils.analysis_utils import delta_r_mask, normalize, trigger_match
from wprime_plus_b.selections.ttbar.jet_selection import select_good_bjets
from wprime_plus_b.corrections.rochester import apply_rochester_corrections
from wprime_plus_b.corrections.lepton import (
    ElectronCorrector,
    MuonCorrector,
    TauCorrector,
)
from wprime_plus_b.selections.ttbar.config import (
    ttbar_electron_selection,
    ttbar_muon_selection,
    ttbar_tau_selection,
    ttbar_jet_selection,
)
from wprime_plus_b.selections.ttbar.lepton_selection import (
    select_good_electrons,
    select_good_muons,
    select_good_taus,
)


class TtbarAnalysis(processor.ProcessorABC):
    """
    Ttbar Analysis processor

    Parameters:
    -----------
    channel:
        region channel {'2b1l', '1b1e1mu', '1b1l'}
    lepton_flavor:
        lepton flavor {'ele', 'mu'}
    year:
        year of the dataset {"2017"}
    year_mode:
        year modifier {""}
    syst:
        systematics to apply {"nominal", "jec", "jer", "jet", "met", "full"} 
    output_type:
        output object type {'hist', 'array'}
    """

    def __init__(
        self,
        channel: str = "2b1l",
        lepton_flavor: str = "ele",
        year: str = "2017",
        yearmod: str = "",
        syst: str = "nominal",
        output_type: str = "hist",
    ):
        self._year = year
        self._yearmod = yearmod
        self._lepton_flavor = lepton_flavor
        self._channel = channel
        self._syst = syst
        self._output_type = output_type

        # define region of the analysis
        self._region = f"{self._channel}_{self._lepton_flavor}"
        # initialize dictionary of hists for control regions
        self.hist_dict = {}
        self.hist_dict[self._region] = {
            "n_kin": histograms.ttbar_n_hist,
            "jet_kin": histograms.ttbar_jet_hist,
            "met_kin": histograms.ttbar_met_hist,
            "lepton_kin": histograms.ttbar_lepton_hist,
            "lepton_bjet_kin": histograms.ttbar_lepton_bjet_hist,
            "lepton_met_kin": histograms.ttbar_lepton_met_hist,
            "lepton_met_bjet_kin": histograms.ttbar_lepton_met_bjet_hist,
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
        # get number of events before selection
        nevents = len(events)
        # check if sample is MC
        self.is_mc = hasattr(events, "genWeight")
        # create copies of histogram objects
        hist_dict = copy.deepcopy(self.hist_dict)
        # create copy of array dictionary
        array_dict = copy.deepcopy(self.array_dict)
        # dictionary to store output data and metadata
        output = {}
        output["metadata"] = {}
        output["metadata"].update({"raw_initial_nevents": nevents})
                    
        # define systematic variations shifts
        syst_variations = ["nominal"]
        if self.is_mc:
            jet_jec_syst_variations = ["JESUp", "JESDown"]
            jet_jer_syst_variations = ["JERUp", "JERDown"]
            met_obj_syst_variations = ["UEUp", "UEDown"]

            if self._syst == "jec":
                syst_variations.extend(jet_jec_syst_variations)
            elif self._syst == "jer":
                syst_variations.extend(jet_jer_syst_variations)
            elif self._syst == "jet":
                syst_variations.extend(jet_jec_syst_variations)
                syst_variations.extend(jet_jer_syst_variations)
            elif self._syst == "met":
                syst_variations.extend(met_obj_syst_variations)
            elif self._syst == "full":
                syst_variations.extend(jet_jec_syst_variations)
                syst_variations.extend(jet_jer_syst_variations)
                syst_variations.extend(met_obj_syst_variations)
                
        for syst_var in syst_variations:
            
            # -------------------------------------------------------------  
            # object corrections
            # -------------------------------------------------------------  
            # apply JEC/JER corrections to jets (in data, the corrections are already applied)
            if self.is_mc:
                corrected_jets, met = jet_corrections(
                    events, self._year + self._yearmod
                )
                # jet JEC/JER shift
                if syst_var == "JESUp":
                    corrected_jets = corrected_jets.JES_Total.up
                elif syst_var == "JESDown":
                    corrected_jets = corrected_jets.JES_Total.down
                elif syst_var == "JERUp":
                    corrected_jets = corrected_jets.JER.up
                elif syst_var == "JERDown":
                    corrected_jets = corrected_jets.JER.down
                # MET UnclusteredEnergy shift
                elif syst_var == "UEUp":
                    met = met.MET_UnclusteredEnergy.up
                elif syst_var == "UEDown":
                    met = met.MET_UnclusteredEnergy.down
            else:
                corrected_jets, met = events.Jet, events.MET
                
            # apply MET phi corrections
            met_pt, met_phi = met_phi_corrections(
                met_pt=met.pt,
                met_phi=met.phi,
                npvs=events.PV.npvsGood,
                run=events.run,
                is_mc=self.is_mc,
                year=self._year,
                year_mod=self._yearmod,
            )
            met["pt"], met["phi"] = met_pt, met_phi
            
            # apply Tau energy corrections (only to MC)
            corrected_taus = events.Tau
            if self.is_mc:
                # Data does not have corrections
                corrected_taus["pt"], corrected_taus["mass"] = tau_energy_scale(
                    events, "2017", "", "DeepTau2017v2p1", "nom"
                )
                # Given the tau corrections. We need to recalculate the MET.
                # https://github.com/columnflow/columnflow/blob/16d35bb2f25f62f9110a8f1089e8dc5c62b29825/columnflow/calibration/util.py#L42
                # https://github.com/Katsch21/hh2bbtautau/blob/e268752454a0ce0089ff08cc6c373a353be77679/hbt/calibration/tau.py#L117
                met["pt"], met["phi"] = met_corrected_tes(
                    old_taus=events.Tau, new_taus=corrected_taus, met=met
                )
                
            # apply rochester corretions to muons
            corrected_muons = events.Muon 
            muon_pt = apply_rochester_corrections(
                corrected_muons, self.is_mc, self._year + self._yearmod
            )
            corrected_muons["pt"] = muon_pt
            met["pt"], met["phi"] = met_corrected_tes(
                events.Muon, corrected_muons, met
            )
            
            # -------------------------------------------------------------  
            # event SF/weights computation
            # -------------------------------------------------------------  
            # get triggers masks 
            with importlib.resources.path("wprime_plus_b.data", "triggers.json") as path:
                with open(path, "r") as handle:
                    self._triggers = json.load(handle)[self._year]
            trigger_mask = {}
            for ch in ["ele", "mu"]:
                trigger_mask[ch] = np.zeros(nevents, dtype="bool")
                for t in self._triggers[ch]:
                    if t in events.HLT.fields:
                        trigger_mask[ch] = trigger_mask[ch] | events.HLT[t]
        
            # get DeltaR matched trigger objects mask
            trigger_path = {
                "1b1l": {
                    "ele": self._triggers["ele"][0],
                    "mu": self._triggers["mu"][0],
                },
                "2b1l": {
                    "ele": self._triggers["ele"][0],
                    "mu": self._triggers["mu"][0],
                },
                "1b1e1mu": {
                    "ele": self._triggers["mu"][0],
                    "mu": self._triggers["ele"][0],
                }
            }
            trigger_leptons = {
                "1b1l": {
                    "ele": events.Electron,
                    "mu": corrected_muons,
                },
                "2b1l": {
                    "ele": events.Electron,
                    "mu": corrected_muons,
                },
                "1b1e1mu": {
                    "ele": corrected_muons,
                    "mu": events.Electron,
                }
            }
            trigger_match_mask = trigger_match(
                leptons=trigger_leptons[self._channel][self._lepton_flavor],
                trigobjs=events.TrigObj,
                trigger_path=trigger_path[self._channel][self._lepton_flavor],
            )
            # weights (for all channels): genweight, pileup, l1prefiring, pujetid, b-tagging
            # electron weights (for 2b1e, 1b1e or 1b1e1mu): electronId, electronReco
            # muon weights (for 2b1mu, 1b1mu, or 1b1e1mu): muonId, muonIso, muonTriggerIso
            weights_container = Weights(len(events), storeIndividual=True)
            if self.is_mc:
                # add gen weigths
                weights_container.add("genweight", events.genWeight)
                # add l1prefiring weigths
                add_l1prefiring_weight(events, weights_container, self._year, syst_var)
                # add pileup weigths
                add_pileup_weight(
                    events, weights_container, self._year, self._yearmod, syst_var
                )
                # add pujetid weigths
                add_pujetid_weight(
                    jets=corrected_jets,
                    genjets=events.GenJet,
                    weights=weights_container,
                    year=self._year,
                    year_mod=self._yearmod,
                    working_point=ttbar_jet_selection[self._channel][
                        self._lepton_flavor
                    ]["jet_pileup_id"],
                    variation=syst_var,
                )
                # b-tagging corrector
                btag_corrector = BTagCorrector(
                    jets=corrected_jets,
                    weights=weights_container,
                    sf_type="comb",
                    worging_point=ttbar_jet_selection[self._channel][
                        self._lepton_flavor
                    ]["btag_working_point"],
                    tagger="deepJet",
                    year=self._year,
                    year_mod=self._yearmod,
                    full_run=False,
                    variation=syst_var,
                )
                # add b-tagging weights
                btag_corrector.add_btag_weights(flavor="bc")
                # electron corrector
                electron_corrector = ElectronCorrector(
                    electrons=events.Electron,
                    weights=weights_container,
                    year=self._year,
                    year_mod=self._yearmod,
                    variation=syst_var,
                )
                # add electron ID weights
                electron_corrector.add_id_weight(
                    id_working_point=ttbar_electron_selection[self._channel][
                        self._lepton_flavor
                    ]["electron_id_wp"]
                )
                # add electron reco weights
                electron_corrector.add_reco_weight()
                # muon corrector
                muon_corrector = MuonCorrector(
                    muons=corrected_muons,
                    weights=weights_container,
                    year=self._year,
                    year_mod=self._yearmod,
                    variation=syst_var,
                    id_wp=ttbar_muon_selection[self._channel][self._lepton_flavor][
                        "muon_id_wp"
                    ],
                    iso_wp=ttbar_muon_selection[self._channel][self._lepton_flavor][
                        "muon_iso_wp"
                    ],
                )
                # add muon ID weights
                muon_corrector.add_id_weight()
                # add muon iso weights
                muon_corrector.add_iso_weight()
                
                # add trigger weights
                if self._channel == "1b1e1mu":
                    if self._lepton_flavor == "ele":
                        muon_corrector.add_triggeriso_weight(
                            trigger_mask=trigger_mask["mu"],
                            trigger_match_mask=trigger_match_mask
                        )
                    else:
                        pass
                        """
                        electron_corrector.add_trigger_weight(
                            trigger_mask=trigger_mask["ele"],
                            trigger_match_mask=trigger_match_mask
                        )
                        """
                else:
                    if self._lepton_flavor == "mu":
                        muon_corrector.add_triggeriso_weight(
                            trigger_mask=trigger_mask["mu"],
                            trigger_match_mask=trigger_match_mask
                        )
                    else:
                        pass
                        """
                        electron_corrector.add_trigger_weight(
                            trigger_mask=trigger_mask["ele"],
                            trigger_match_mask=trigger_match_mask
                        )
                        """
                # add tau weights
                tau_corrector = TauCorrector(
                    taus=corrected_taus,
                    weights=weights_container,
                    year=self._year,
                    year_mod=self._yearmod,
                    tau_vs_jet=ttbar_tau_selection[self._channel][self._lepton_flavor][
                        "tau_vs_jet"
                    ],
                    tau_vs_ele=ttbar_tau_selection[self._channel][self._lepton_flavor][
                        "tau_vs_ele"
                    ],
                    tau_vs_mu=ttbar_tau_selection[self._channel][self._lepton_flavor][
                        "tau_vs_mu"
                    ],
                    variation=syst_var,
                )
                tau_corrector.add_id_weight_DeepTau2017v2p1VSe()
                tau_corrector.add_id_weight_DeepTau2017v2p1VSmu()
                tau_corrector.add_id_weight_DeepTau2017v2p1VSjet()

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
                electron_pt_threshold=ttbar_electron_selection[self._channel][
                    self._lepton_flavor
                ]["electron_pt_threshold"],
                electron_id_wp=ttbar_electron_selection[self._channel][
                    self._lepton_flavor
                ]["electron_id_wp"],
                electron_iso_wp=ttbar_electron_selection[self._channel][
                    self._lepton_flavor
                ]["electron_iso_wp"],
            )
            electrons = events.Electron[good_electrons]
            
            # select good muons
            good_muons = select_good_muons(
                muons=corrected_muons,
                muon_pt_threshold=ttbar_muon_selection[self._channel][
                    self._lepton_flavor
                ]["muon_pt_threshold"],
                muon_id_wp=ttbar_muon_selection[self._channel][self._lepton_flavor][
                    "muon_id_wp"
                ],
                muon_iso_wp=ttbar_muon_selection[self._channel][self._lepton_flavor][
                    "muon_iso_wp"
                ],
            )
            good_muons = (good_muons) & (
                delta_r_mask(corrected_muons, electrons, threshold=0.4)
            )
            muons = corrected_muons[good_muons]
            
            # select good taus
            good_taus = select_good_taus(
                taus=corrected_taus,
                tau_pt_threshold=ttbar_tau_selection[self._channel][
                    self._lepton_flavor
                ]["tau_pt_threshold"],
                tau_eta_threshold=ttbar_tau_selection[self._channel][
                    self._lepton_flavor
                ]["tau_eta_threshold"],
                tau_dz_threshold=ttbar_tau_selection[self._channel][
                    self._lepton_flavor
                ]["tau_dz_threshold"],
                tau_vs_jet=ttbar_tau_selection[self._channel][self._lepton_flavor][
                    "tau_vs_jet"
                ],
                tau_vs_ele=ttbar_tau_selection[self._channel][self._lepton_flavor][
                    "tau_vs_ele"
                ],
                tau_vs_mu=ttbar_tau_selection[self._channel][self._lepton_flavor][
                    "tau_vs_mu"
                ],
                prong=ttbar_tau_selection[self._channel][self._lepton_flavor]["prongs"],
            )
            good_taus = (
                (good_taus)
                & (delta_r_mask(corrected_taus, electrons, threshold=0.4))
                & (delta_r_mask(corrected_taus, muons, threshold=0.4))
            )
            taus = corrected_taus[good_taus]
            
            # select good bjets
            good_bjets = select_good_bjets(
                jets=corrected_jets,
                year=self._year,
                btag_working_point=ttbar_jet_selection[self._channel][
                    self._lepton_flavor
                ]["btag_working_point"],
                jet_pt_threshold=ttbar_jet_selection[self._channel][
                    self._lepton_flavor
                ]["jet_pt_threshold"],
                jet_id=ttbar_jet_selection[self._channel][self._lepton_flavor][
                    "jet_id"
                ],
                jet_pileup_id=ttbar_jet_selection[self._channel][self._lepton_flavor][
                    "jet_pileup_id"
                ],
            )
            good_bjets = (
                good_bjets
                & (delta_r_mask(corrected_jets, electrons, threshold=0.4))
                & (delta_r_mask(corrected_jets, muons, threshold=0.4))
                & (delta_r_mask(corrected_jets, taus, threshold=0.4))
            )
            bjets = corrected_jets[good_bjets]
            
            
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
                lumi_mask = self._lumi_mask[self._year](
                    events.run, events.luminosityBlock
                )
            else:
                lumi_mask = np.ones(len(events), dtype="bool")
            self.selections.add("lumi", lumi_mask)

            # add lepton triggers masks
            self.selections.add("trigger_ele", trigger_mask["ele"])
            self.selections.add("trigger_mu", trigger_mask["mu"])

            # add MET filters mask
            with importlib.resources.path(
                "wprime_plus_b.data", "metfilters.json"
            ) as path:
                with open(path, "r") as handle:
                    self._metfilters = json.load(handle)[self._year]
            metfilters = np.ones(nevents, dtype="bool")
            metfilterkey = "mc" if self.is_mc else "data"
            for mf in self._metfilters[metfilterkey]:
                if mf in events.Flag.fields:
                    metfilters = metfilters & events.Flag[mf]
            self.selections.add("metfilters", metfilters)

            # check that there be a minimum MET greater than 50 GeV
            self.selections.add("met_pt", met.pt > 50)

            # add number of leptons and jets
            self.selections.add("one_electron", ak.num(electrons) == 1)
            self.selections.add("electron_veto", ak.num(electrons) == 0)
            self.selections.add("one_muon", ak.num(muons) == 1)
            self.selections.add("muon_veto", ak.num(muons) == 0)
            self.selections.add("tau_veto", ak.num(taus) == 0)
            self.selections.add("one_bjet", ak.num(bjets) == 1)
            self.selections.add("two_bjets", ak.num(bjets) == 2)

            # select events with at least one good vertex
            self.selections.add("goodvertex", events.PV.npvsGood > 0)
            
            # select events with at least one matched trigger object
            self.selections.add("trigger_match", ak.sum(trigger_match_mask, axis=-1) > 0)

            # define selection regions for each channel
            region_selection = {
                "2b1l": {
                    "ele": [
                        "goodvertex",
                        "lumi",
                        "trigger_ele",
                        "trigger_match",
                        "metfilters",
                        "met_pt",
                        "two_bjets",
                        "tau_veto",
                        "muon_veto",
                        "one_electron",
                    ],
                    "mu": [
                        "goodvertex",
                        "lumi",
                        "trigger_mu",
                        "trigger_match",
                        "metfilters",
                        "met_pt",
                        "two_bjets",
                        "tau_veto",
                        "electron_veto",
                        "one_muon",
                    ],
                },
                "1b1e1mu": {
                    "ele": [
                        "goodvertex",
                        "lumi",
                        "trigger_mu",
                        "trigger_match",
                        "metfilters",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "one_muon",
                        "one_electron",
                    ],
                    "mu": [
                        "goodvertex",
                        "lumi",
                        "trigger_ele",
                        "trigger_match",
                        "metfilters",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "one_electron",
                        "one_muon",
                    ],
                },
                "1b1l": {
                    "ele": [
                        "goodvertex",
                        "lumi",
                        "trigger_ele",
                        "trigger_match",
                        "metfilters",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "muon_veto",
                        "one_electron",
                    ],
                    "mu": [
                        "goodvertex",
                        "lumi",
                        "trigger_mu",
                        "trigger_match",
                        "metfilters",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "electron_veto",
                        "one_muon",
                    ],
                },
            }
            # save cutflow
            if syst_var == "nominal":
                cut_names = region_selection[self._channel][self._lepton_flavor]
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
            self.selections.add(
                self._region,
                self.selections.all(
                    *region_selection[self._channel][self._lepton_flavor]
                ),
            )
            region_selection = self.selections.all(self._region)
            # check that there are events left after selection
            nevents_after = ak.sum(region_selection)
            if nevents_after > 0:
                # select region objects
                region_bjets = bjets[region_selection]
                region_electrons = electrons[region_selection]
                region_muons = muons[region_selection]
                region_met = met[region_selection]
                
                # define region leptons
                region_leptons = (
                    region_electrons if self._lepton_flavor == "ele" else region_muons
                )
                # lepton relative isolation
                lepton_reliso = (
                    region_leptons.pfRelIso04_all
                    if hasattr(region_leptons, "pfRelIso04_all")
                    else region_leptons.pfRelIso03_all
                )
                # leading bjets
                leading_bjets = ak.firsts(region_bjets)
                # lepton-bjet deltaR and invariant mass
                lepton_bjet_dr = leading_bjets.delta_r(region_leptons)
                lepton_bjet_mass = (region_leptons + leading_bjets).mass
                # lepton-MET transverse mass and deltaPhi
                lepton_met_mass = np.sqrt(
                    2.0
                    * region_leptons.pt
                    * region_met.pt
                    * (
                        ak.ones_like(region_met.pt)
                        - np.cos(region_leptons.delta_phi(region_met))
                    )
                )
                lepton_met_delta_phi = np.abs(region_leptons.delta_phi(region_met))
                # lepton-bJet-MET total transverse mass
                lepton_met_bjet_mass = np.sqrt(
                    (region_leptons.pt + leading_bjets.pt + region_met.pt) ** 2
                    - (region_leptons + leading_bjets + region_met).pt ** 2
                )

                self.add_feature("lepton_pt", region_leptons.pt)
                self.add_feature("lepton_eta", region_leptons.eta)
                self.add_feature("lepton_phi", region_leptons.phi)
                self.add_feature("jet_pt", leading_bjets.pt)
                self.add_feature("jet_eta", leading_bjets.eta)
                self.add_feature("jet_phi", leading_bjets.phi)
                self.add_feature("met", region_met.pt)
                self.add_feature("met_phi", region_met.phi)
                self.add_feature("lepton_bjet_dr", lepton_bjet_dr)
                self.add_feature("lepton_bjet_mass", lepton_bjet_mass)
                self.add_feature("lepton_met_mass", lepton_met_mass)
                self.add_feature("lepton_met_delta_phi", lepton_met_delta_phi)
                self.add_feature("lepton_met_bjet_mass", lepton_met_bjet_mass)
                self.add_feature("njets", ak.num(corrected_jets)[region_selection])
                self.add_feature("npvs", events.PV.npvsGood[region_selection])

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
                if self._output_type == "hist":
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
                            for kin in hist_dict[self._region]:
                                fill_args = {
                                    feature: normalize(self.features[feature])
                                    for feature in hist_dict[self._region][
                                        kin
                                    ].axes.name
                                    if feature not in ["variation"]
                                }
                                hist_dict[self._region][kin].fill(
                                    **fill_args,
                                    variation=variation,
                                    weight=region_weight,
                                )
                    elif self.is_mc and syst_var != "nominal":
                        # object-wise variations
                        region_weight = weights_container.weight()[region_selection]
                        for kin in hist_dict[self._region]:
                            # get filling arguments

                            fill_args = {
                                feature: normalize(self.features[feature])
                                for feature in hist_dict[self._region][kin].axes.name[
                                    :-1
                                ]
                                if feature not in ["variation"]
                            }
                            # fill histograms
                            hist_dict[self._region][kin].fill(
                                **fill_args,
                                variation=syst_var,
                                weight=region_weight,
                            )
                    elif not self.is_mc and syst_var == "nominal":
                        # object-wise variations
                        region_weight = weights_container.weight()[region_selection]
                        for kin in hist_dict[self._region]:
                            # get filling arguments
                            fill_args = {
                                feature: normalize(self.features[feature])
                                for feature in hist_dict[self._region][kin].axes.name[
                                    :-1
                                ]
                                if feature not in ["variation"]
                            }
                            # fill histograms
                            hist_dict[self._region][kin].fill(
                                **fill_args,
                                variation=syst_var,
                                weight=region_weight,
                            )
                elif self._output_type == "array":
                    array_dict = {}
                    self.add_feature(
                        "weights", weights_container.weight()[region_selection]
                    )
                    # uncoment next two lines to save individual weights
                    # for weight in weights_container.weightStatistics:
                    #    self.add_feature(weight, weights_container.partial_weight(include=[weight]))
                    if syst_var == "nominal":
                        # select variables and put them in column accumulators
                        array_dict.update({
                            feature_name: processor.column_accumulator(
                                normalize(feature_array)
                            )
                            for feature_name, feature_array in self.features.items()
                        })
        # define output dictionary accumulator
        if self._output_type == "hist":
            output["histograms"] = hist_dict[f"{self._channel}_{self._lepton_flavor}"]
        elif self._output_type == "array":
            output["arrays"] = array_dict
            
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator