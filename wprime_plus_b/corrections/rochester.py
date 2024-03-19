import numpy as np
import awkward as ak
from coffea.lookup_tools import txt_converters, rochester_lookup


def apply_rochester_corrections(muons, is_mc, year):
    rochester_data = txt_converters.convert_rochester_file(
        f"wprime_plus_b/data/RoccoR{year}UL.txt", loaduncs=True
    )
    rochester = rochester_lookup.rochester_lookup(rochester_data)

    if is_mc:
        hasgen = ~np.isnan(ak.fill_none(muons.matched_gen.pt, np.nan))
        mc_rand = np.random.rand(*ak.to_numpy(ak.flatten(muons.pt)).shape)
        mc_rand = ak.unflatten(mc_rand, ak.num(muons.pt, axis=1))
        corrections = np.array(ak.flatten(ak.ones_like(muons.pt)))
        mc_kspread = rochester.kSpreadMC(
            muons.charge[hasgen],
            muons.pt[hasgen],
            muons.eta[hasgen],
            muons.phi[hasgen],
            muons.matched_gen.pt[hasgen],
        )
        mc_ksmear = rochester.kSmearMC(
            muons.charge[~hasgen],
            muons.pt[~hasgen],
            muons.eta[~hasgen],
            muons.phi[~hasgen],
            muons.nTrackerLayers[~hasgen],
            mc_rand[~hasgen],
        )
        hasgen_flat = np.array(ak.flatten(hasgen))
        corrections[hasgen_flat] = np.array(ak.flatten(mc_kspread))
        corrections[~hasgen_flat] = np.array(ak.flatten(mc_ksmear))
        corrections = ak.unflatten(corrections, ak.num(muons.pt, axis=1))
    else:
        corrections = rochester.kScaleDT(muons.charge, muons.pt, muons.eta, muons.phi)
    return muons.pt * corrections
