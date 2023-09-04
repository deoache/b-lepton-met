import json
import copy
import pickle
import numpy as np
import awkward as ak
import importlib.resources
from coffea import processor
from coffea.analysis_tools import Weights, PackedSelection
from wprime_plus_b.processors.utils import histograms
from wprime_plus_b.processors.utils.analysis_utils import delta_r_mask, normalize
from wprime_plus_b.corrections.btag import BTagCorrector
from wprime_plus_b.corrections.jec import jet_corrections
from wprime_plus_b.corrections.met import met_phi_corrections
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.lepton import ElectronCorrector, MuonCorrector
from wprime_plus_b.selections.ttbar.jet_selection import select_good_bjets
from wprime_plus_b.selections.ttbar.config import (
    electron_selection,
    muon_selection,
    jet_selection,
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
        year of the dataset {"2016", "2017", "2018"}
    year_mode:
        year modifier {"", "APV"}
    btag_wp:
        working point of the deepJet tagger
    syst:
        systematics to apply
    output_type:

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

        # define systematic variations
        syst_variations = ["nominal"]
        if self.is_mc:
            jet_jec_syst_variations = ["JESUp", "JESDown"]
            jet_jer_syst_variations = ["JERUp", "JERDown"]
            met_obj_syst_variations = ["UEUp", "UEDown"]

            if self._syst == "jec":
                syst_variations.extend(jet_jec_syst_variations)
            elif self._syst == "jer":
                syst_variations.extend(jet_jer_syst_variations)
            elif self._syst == "jec":
                syst_variations.extend(jet_jec_syst_variations)
                syst_variations.extend(jet_jer_syst_variations)
            elif self._syst == "met":
                syst_variations.extend(met_obj_syst_variations)
            elif self._syst == "full":
                syst_variations.extend(jet_jec_syst_variations)
                syst_variations.extend(jet_jer_syst_variations)
                syst_variations.extend(met_obj_syst_variations)

        for syst_var in syst_variations:
            # ------------------
            # event preselection
            # ------------------

            # ------------------
            # leptons
            # -------------------
            # select good electrons
            good_electrons = select_good_electrons(
                events=events,
                electron_pt_threshold=electron_selection[self._channel][
                    self._lepton_flavor
                ]["electron_pt_threshold"],
                electron_id_wp=electron_selection[self._channel][self._lepton_flavor][
                    "electron_id_wp"
                ],
                electron_iso_wp=electron_selection[self._channel][self._lepton_flavor][
                    "electron_iso_wp"
                ],
            )
            electrons = events.Electron[good_electrons]

            # select good muons
            good_muons = select_good_muons(
                events=events,
                muon_pt_threshold=muon_selection[self._channel][self._lepton_flavor][
                    "muon_pt_threshold"
                ],
                muon_id_wp=muon_selection[self._channel][self._lepton_flavor][
                    "muon_id_wp"
                ],
                muon_iso_wp=muon_selection[self._channel][self._lepton_flavor][
                    "muon_iso_wp"
                ],
            )
            good_muons = (good_muons) & (
                delta_r_mask(events.Muon, electrons, threshold=0.4)
            )
            muons = events.Muon[good_muons]

            # select good taus
            good_taus = (
                select_good_taus(events)
                & (delta_r_mask(events.Tau, electrons, threshold=0.4))
                & (delta_r_mask(events.Tau, muons, threshold=0.4))
            )
            taus = events.Tau[good_taus]

            # ------------------
            # jets
            # -------------------
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

            # select good bjets
            good_bjets = select_good_bjets(
                events=events,
                year=self._year,
                btag_working_point=jet_selection[self._channel][self._lepton_flavor][
                    "btag_working_point"
                ],
                jet_pt_threshold=jet_selection[self._channel][self._lepton_flavor][
                    "jet_pt_threshold"
                ],
                jet_id=jet_selection[self._channel][self._lepton_flavor]["jet_id"],
                jet_pileup_id=jet_selection[self._channel][self._lepton_flavor][
                    "jet_pileup_id"
                ],
            )
            good_bjets = (
                good_bjets
                & (delta_r_mask(corrected_jets, electrons, threshold=0.4))
                & (delta_r_mask(corrected_jets, muons, threshold=0.4))
            )
            bjets = corrected_jets[good_bjets]

            # apply MET phi corrections
            met_pt, met_phi = met_phi_corrections(
                met_pt=met.pt,
                met_phi=met.phi,
                npvs=events.PV.npvs,
                is_mc=self.is_mc,
                year=self._year,
                year_mod=self._yearmod,
            )
            met["pt"], met["phi"] = met_pt, met_phi

            # ---------------
            # event selection
            # ---------------
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
            with importlib.resources.path(
                "wprime_plus_b.data", "triggers.json"
            ) as path:
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
            # open and load met filters
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

            # define selection regions for each channel
            region_selection = {
                "2b1l": {
                    "ele": [
                        "lumi",
                        "metfilters",
                        "trigger_ele",
                        "met_pt",
                        "two_bjets",
                        "tau_veto",
                        "muon_veto",
                        "one_electron",
                    ],
                    "mu": [
                        "lumi",
                        "metfilters",
                        "trigger_mu",
                        "met_pt",
                        "two_bjets",
                        "tau_veto",
                        "electron_veto",
                        "one_muon",
                    ],
                },
                "1b1e1mu": {
                    "ele": [
                        "lumi",
                        "metfilters",
                        "trigger_mu",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "one_muon",
                        "one_electron",
                    ],
                    "mu": [
                        "lumi",
                        "metfilters",
                        "trigger_ele",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "one_electron",
                        "one_muon",
                    ],
                },
                "1b1l": {
                    "ele": [
                        "lumi",
                        "metfilters",
                        "trigger_ele",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "muon_veto",
                        "one_electron",
                    ],
                    "mu": [
                        "lumi",
                        "metfilters",
                        "trigger_mu",
                        "met_pt",
                        "one_bjet",
                        "tau_veto",
                        "electron_veto",
                        "one_muon",
                    ],
                },
            }
            # ---------------
            # event variables
            # ---------------
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

                # -------------
                # event weights
                # -------------
                weights_container = Weights(
                    len(events[region_selection]), storeIndividual=True
                )
                if self.is_mc:
                    # add gen weigths
                    gen_weight = events.genWeight[region_selection]
                    weights_container.add("genweight", gen_weight)

                    # add L1prefiring weights
                    if self._year in ("2016", "2017"):
                        if syst_var == "nominal":
                            weights_container.add(
                                "L1Prefiring",
                                weight=events.L1PreFiringWeight.Nom[region_selection],
                                weightUp=events.L1PreFiringWeight.Up[region_selection],
                                weightDown=events.L1PreFiringWeight.Dn[
                                    region_selection
                                ],
                            )
                        else:
                            weights_container.add(
                                "L1Prefiring",
                                weight=events.L1PreFiringWeight.Nom[region_selection],
                            )
                    # add pileup reweighting
                    add_pileup_weight(
                        n_true_interactions=ak.to_numpy(
                            events.Pileup.nPU[region_selection]
                        ),
                        weights=weights_container,
                        year=self._year,
                        year_mod=self._yearmod,
                        variation=syst_var,
                    )
                    # b-tagging corrector
                    btag_corrector = BTagCorrector(
                        jets=region_bjets,
                        njets=2 if self._channel == "2b1l" else 1,
                        weights=weights_container,
                        sf_type="comb",
                        worging_point=jet_selection[self._channel][self._lepton_flavor][
                            "btag_working_point"
                        ],
                        tagger="deepJet",
                        year=self._year,
                        year_mod=self._yearmod,
                        full_run=False,
                        variation=syst_var,
                    )
                    # add b-tagging weights
                    btag_corrector.add_btag_weights(flavor="bc")

                    # electron corrector
                    if (self._channel == "1b1e1mu") or (self._lepton_flavor == "ele"):
                        electron_corrector = ElectronCorrector(
                            electrons=ak.firsts(region_electrons),
                            weights=weights_container,
                            year=self._year,
                            year_mod=self._yearmod,
                            tag="leading_electron",
                            variation=syst_var,
                        )
                        # add electron ID weights
                        electron_corrector.add_id_weight(
                            id_working_point=electron_selection[self._channel][
                                self._lepton_flavor
                            ]["electron_id_wp"]
                        )
                        # add electron reco weights
                        electron_corrector.add_reco_weight()
                    # muon corrector
                    if (self._channel == "1b1e1mu") or (self._lepton_flavor == "mu"):
                        muon_corrector = MuonCorrector(
                            muons=ak.firsts(region_muons),
                            weights=weights_container,
                            year=self._year,
                            year_mod=self._yearmod,
                            tag="leading_muon",
                            variation=syst_var,
                            id_wp=muon_selection[self._channel][self._lepton_flavor][
                                "muon_id_wp"
                            ],
                            iso_wp=muon_selection[self._channel][self._lepton_flavor][
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
                            muon_corrector.add_triggeriso_weight()
                        else:
                            pass  # electron_corrector.add_trigger_weight()
                    else:
                        if self._lepton_flavor == "ele":
                            pass  # electron_corrector.add_trigger_weight()
                        else:
                            muon_corrector.add_triggeriso_weight()

                # ------------------
                # histogram filling
                # ------------------
                if self._output_type == "hist":
                    # break up the histogram filling for event-wise variations and object-wise variations
                    # apply event-wise variations only for nominal
                    if self.is_mc and syst_var == "nominal":
                        # get event weight systematic variations for MC samples
                        event_weights = [
                            weight
                            for weight in weights_container.weightStatistics
                            if "genweight" not in weight
                        ]
                        event_weight_syst_variations_up = [
                            f"{event_weight}Up" for event_weight in event_weights
                        ]
                        event_weight_syst_variations_down = [
                            f"{event_weight}Down" for event_weight in event_weights
                        ]
                        event_weight_syst = ["nominal"]
                        event_weight_syst.extend(event_weight_syst_variations_up)
                        event_weight_syst.extend(event_weight_syst_variations_down)

                        for variation in event_weight_syst:
                            # get weight
                            if variation == "nominal":
                                syst_weight = weights_container.weight()
                            else:
                                syst_weight = weights_container.weight(variation)

                            for kin in hist_dict[self._region]:
                                # get filling arguments
                                fill_args = {
                                    feature: normalize(self.features[feature])
                                    for feature in hist_dict[self._region][
                                        kin
                                    ].axes.name[:-1]
                                    if "dataset" not in feature
                                }
                                # fill histograms
                                hist_dict[self._region][kin].fill(
                                    **fill_args,
                                    dataset=dataset,
                                    variation=variation,
                                    weight=syst_weight,
                                )
                    # object-wise variations
                    syst_weight = weights_container.weight()
                    for kin in hist_dict[self._region]:
                        # get filling arguments
                        fill_args = {
                            feature: normalize(self.features[feature])
                            for feature in hist_dict[self._region][kin].axes.name[:-1]
                            if "dataset" not in feature
                        }
                        # fill histograms
                        hist_dict[self._region][kin].fill(
                            **fill_args,
                            dataset=dataset,
                            variation=syst_var,
                            weight=syst_weight,
                        )
                elif self._output_type == "array":
                    self.add_feature("weights", weights_container.weight())
                    self.add_feature("genweights", weights_container.partial_weight("genweight"))
                    # select variables and put them in column accumulators
                    array_dict = {
                        feature_name: processor.column_accumulator(
                            normalize(feature_array)
                        )
                        for feature_name, feature_array in self.features.items()
                    }
        # define output dictionary accumulator
        output = {}
        if self._output_type == "hist":
            output["histograms"] = hist_dict[f"{self._channel}_{self._lepton_flavor}"]
        elif self._output_type == "array":
            output["arrays"] = array_dict
        # save metadata
        output["metadata"] = {
            "events_before": nevents,
            "events_after": nevents_after,
        }
        # save sumw for MC samples
        if self.is_mc:
            output["metadata"].update({"sumw": ak.sum(events.genWeight)})
        return output

    def postprocess(self, accumulator):
        return accumulator