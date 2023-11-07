import awkward as ak
from coffea.analysis_tools import Weights
from wprime_plus_b.corrections.btag import BTagCorrector
from wprime_plus_b.corrections.pileup import add_pileup_weight
from wprime_plus_b.corrections.pujetid import add_pujetid_weight
from wprime_plus_b.corrections.lepton import ElectronCorrector, MuonCorrector


def add_l1prefiring_sf(events, weights_container, year="2017", variation="nominal"):
    # add L1prefiring weights
    if year in ("2016", "2017"):
        if variation == "nominal":
            weights_container.add(
                "L1Prefiring",
                weight=events.L1PreFiringWeight.Nom,
                weightUp=events.L1PreFiringWeight.Up,
                weightDown=events.L1PreFiringWeight.Dn,
            )
        else:
            weights_container.add(
                "L1Prefiring",
                weight=events.L1PreFiringWeight.Nom,
            )


def add_pileup_sf(
    events, weights_container, year="2017", yearmod="", variation="nominal"
):
    add_pileup_weight(
        n_true_interactions=ak.to_numpy(events.Pileup.nPU),
        weights=weights_container,
        year=year,
        year_mod=yearmod,
        variation=variation,
    )


def add_pujetid_sf(
    bjets,
    weights_container,
    bjets_selection,
    year="2017",
    yearmod="",
    variation="nominal",
    channel="2b1l",
    lepton_flavor="mu",
):
    add_pujetid_weight(
        jets=bjets,
        weights=weights_container,
        year=year,
        year_mod=yearmod,
        working_point=bjets_selection[channel][lepton_flavor]["btag_working_point"],
        variation=variation,
    )


def add_btag_sf(
    bjets,
    weights_container,
    bjets_selection,
    year="2017",
    yearmod="",
    variation="nominal",
    channel="2b1l",
    lepton_flavor="mu",
):
    # b-tagging corrector
    btag_corrector = BTagCorrector(
        jets=bjets,
        weights=weights_container,
        sf_type="comb",
        worging_point=bjets_selection[channel][lepton_flavor]["btag_working_point"],
        tagger="deepJet",
        year=year,
        year_mod=yearmod,
        full_run=False,
        variation=variation,
    )
    # add b-tagging weights
    btag_corrector.add_btag_weights(flavor="bc")


def add_lepton_sf(
    electrons,
    muons,
    weights_container,
    electron_selection,
    muon_selection,
    year="2017",
    yearmod="",
    variation="nominal",
    channel="2b1l",
    lepton_flavor="mu",
):
    if (channel == "1b1e1mu") or (lepton_flavor == "ele"):
        electron_corrector = ElectronCorrector(
            electrons=electrons,
            weights=weights_container,
            year=year,
            year_mod=yearmod,
            variation=variation,
        )
        # add electron ID weights
        electron_corrector.add_id_weight(
            id_working_point=electron_selection[channel][lepton_flavor][
                "electron_id_wp"
            ]
        )
        # add electron reco weights
        electron_corrector.add_reco_weight()

    # muon corrector
    if (channel == "1b1e1mu") or (lepton_flavor == "mu"):
        muon_corrector = MuonCorrector(
            muons=muons,
            weights=weights_container,
            year=year,
            year_mod=yearmod,
            variation=variation,
            id_wp=muon_selection[channel][lepton_flavor]["muon_id_wp"],
            iso_wp=muon_selection[channel][lepton_flavor]["muon_iso_wp"],
        )
        # add muon ID weights
        muon_corrector.add_id_weight()

        # add muon iso weights
        muon_corrector.add_iso_weight()

    # add trigger weights
    if channel == "1b1e1mu":
        if lepton_flavor == "ele":
            muon_corrector.add_triggeriso_weight()
        else:
            pass  # electron_corrector.add_trigger_weight()
    else:
        if lepton_flavor == "ele":
            pass  # electron_corrector.add_trigger_weight()
        else:
            muon_corrector.add_triggeriso_weight()


def event_weights(
    events,
    electrons,
    muons,
    bjets,
    met,
    electron_selection,
    muon_selection,
    bjets_selection,
    is_mc=True,
    channel="2b1l",
    lepton_flavor="mu",
    year="2017",
    yearmod="",
    variation="nominal",
):
    weights_container = Weights(len(events), storeIndividual=True)
    if is_mc:
        # add gen weigths
        weights_container.add("genweight", events.genWeight)

        # add l1prefiring weigths
        add_l1prefiring_sf(events, weights_container, year, variation)

        # add pileup weigths
        add_pileup_sf(events, weights_container, year, yearmod, variation)

        # add pujetid weigths
        add_pujetid_sf(
            bjets,
            weights_container,
            bjets_selection,
            year,
            yearmod,
            variation,
            channel,
            lepton_flavor,
        )

        # add b-tagging weigths
        add_btag_sf(
            bjets,
            weights_container,
            bjets_selection,
            year,
            yearmod,
            variation,
            channel,
            lepton_flavor,
        )

        # add leptons weigths
        add_lepton_sf(
            electrons,
            muons,
            weights_container,
            electron_selection,
            muon_selection,
            year,
            yearmod,
            variation,
            channel,
            lepton_flavor,
        )

    return weights_container