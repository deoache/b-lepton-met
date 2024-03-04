import copy
import correctionlib
import numpy as np
import awkward as ak
from wprime_plus_b.corrections.utils import get_pog_json

# ----------------------------------------------------------------------------------- #
# -- The tau energy scale (TES) corrections for taus are provided  ------------------ #
# --  to be applied to reconstructed tau_h Lorentz vector ----------------------------#
# --  It should be applied to a genuine tau -> genmatch = 5 --------------------------#
# -----------  (pT, mass and energy) in simulated data -------------------------------#
# tau_E  *= tes
# tau_pt *= tes
# tau_m  *= tes
# https://github.com/cms-tau-pog/TauIDSFs/tree/master
# ----------------------------------------------------------------------------------- #


def mask_energy_corrections(tau):
    # https://github.com/cms-tau-pog/TauFW/blob/4056e9dec257b9f68d1a729c00aecc8e3e6bf97d/PicoProducer/python/analysis/ETauFakeRate/ModuleETau.py#L320
    # https://gitlab.cern.ch/cms-tau-pog/jsonpog-integration/-/blob/TauPOG_v2/POG/TAU/scripts/tau_tes.py
    tau_mask_gm = (
        (tau.genPartFlav == 5)  # Genuine tau
        | (tau.genPartFlav == 1)  # e -> fake
        | (tau.genPartFlav == 2)  # mu -> fake
        | (tau.genPartFlav == 6)  # unmached
    )
    tau_mask_dm = (
        (tau.decayMode == 0)
        | (tau.decayMode == 1)  # 1 prong
        | (tau.decayMode == 2)  # 1 prong
        | (tau.decayMode == 10)  # 1 prong
        | (tau.decayMode == 11)  # 3 prongs  # 3 prongs
    )
    # I have change eta -> np.abs(eta)
    tau_eta_mask = (tau.eta >= 0) & (tau.eta < 2.5)
    tau_mask = tau_mask_gm & tau_mask_dm  # & tau_eta_mask
    return tau_mask


def tau_energy_scale(
    events: ak.Array,
    year: str = "2017",
    year_mod: str = "",
    id: str = "DeepTau2017v2p1",
    sys: str = "nom",
):

    # Corrections works with flatten values
    ntaus = ak.num(copy.deepcopy(events.Tau))
    taus_flatten = ak.flatten(copy.deepcopy(events.Tau))
    # It is defined the taus will be corrected with the energy scale factor: Only a subset of the initial taus.
    mask = mask_energy_corrections(taus_flatten)
    taus_filter = taus_flatten.mask[mask]
    # We fill the None values with valid entries
    pt = ak.fill_none(taus_filter.pt, 0)
    eta = ak.fill_none(taus_filter.eta, 0)
    dm = ak.fill_none(taus_filter.decayMode, 0)
    genmatch = ak.fill_none(taus_filter.genPartFlav, 2)
    # Define correction set_id
    cset = correctionlib.CorrectionSet.from_file(
        get_pog_json(json_name="tau", year=year + year_mod)
    )
    SF = cset["tau_energy_scale"].evaluate(pt, eta, dm, genmatch, id, sys)
    taus_new_pt = taus_filter.pt * SF
    taus_new_mass = taus_filter.mass * SF
    new_tau_pt = ak.where(mask, taus_new_pt, taus_flatten.pt)
    new_tau_mass = ak.where(mask, taus_new_mass, taus_flatten.mass)
    # We have to unflatten the taus
    tau_pt = ak.unflatten(new_tau_pt, ntaus)
    tau_mass = ak.unflatten(new_tau_mass, ntaus)
    return tau_pt, tau_mass


def met_corrected_tes(old_taus: ak.Array, new_taus: ak.Array, met: ak.Array):
    # Helper function to compute new MET based on tau pts corrections.

    # build px and py sums before and after: we sum the time at x and the time at y of each event
    # Initial
    old_px_tau = old_taus.pt * np.cos(old_taus.phi)
    old_py_tau = old_taus.pt * np.sin(old_taus.phi)
    old_px = ak.sum(old_px_tau, axis=1)
    old_py = ak.sum(old_py_tau, axis=1)
    # Final
    new_px_tau = new_taus.pt * np.cos(new_taus.phi)
    new_py_tau = new_taus.pt * np.sin(new_taus.phi)
    new_px = ak.sum(new_px_tau, axis=1)
    new_py = ak.sum(new_py_tau, axis=1)
    # Increase
    Delta_x = new_px - old_px
    Delta_y = new_py - old_py
    # propagate to met
    met_px = met.pt * np.cos(met.phi) - Delta_x
    met_py = met.pt * np.sin(met.phi) - Delta_y
    # compute new components
    met_pt = np.sqrt((met_px**2.0 + met_py**2.0))
    met_phi = np.arctan2(met_py, met_px)
    return met_pt, met_phi
