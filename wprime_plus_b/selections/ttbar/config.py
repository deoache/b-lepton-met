ttbar_electron_selection = {
    "1b1e1mu": {
        "ele": {
            "electron_pt_threshold": 30,
            "electron_id_wp": "wp80iso",
            "electron_iso_wp": None,
        },
        "mu": {
            "electron_pt_threshold": 55,
            "electron_id_wp": "wp90iso",
            "electron_iso_wp": None,
        },
    },
    "2b1l": {
        "ele": {
            "electron_pt_threshold": 55,
            "electron_id_wp": "wp80iso",
            "electron_iso_wp": None,
        },
        "mu": {
            "electron_pt_threshold": 30,
            "electron_id_wp": "wp90iso",
            "electron_iso_wp": None,
        },
    },
    "1b1l": {
        "ele": {
            "electron_pt_threshold": 55,
            "electron_id_wp": "wp80iso",
            "electron_iso_wp": None,
        },
        "mu": {
            "electron_pt_threshold": 50,
            "electron_id_wp": "wp90iso",
            "electron_iso_wp": None,
        },
    },
}

ttbar_muon_selection = {
    "1b1e1mu": {
        "ele": {
            "muon_pt_threshold": 35,
            "muon_id_wp": "tight",
            "muon_iso_wp": "tight",
        },
        "mu": {
            "muon_pt_threshold": 35,
            "muon_id_wp": "tight",
            "muon_iso_wp": "tight",
        },
    },
    "2b1l": {
        "ele": {
            "muon_pt_threshold": 35,
            "muon_id_wp": "tight",
            "muon_iso_wp": "tight",
        },
        "mu": {
            "muon_pt_threshold": 35,
            "muon_id_wp": "tight",
            "muon_iso_wp": "tight",
        },
    },
    "1b1l": {
        "ele": {
            "muon_pt_threshold": 35,
            "muon_id_wp": "tight",
            "muon_iso_wp": "tight",
        },
        "mu": {
            "muon_pt_threshold": 35,
            "muon_id_wp": "tight",
            "muon_iso_wp": "tight",
        },
    },
}

ttbar_tau_selection = {
    "2b1l": {
        "mu": {
            "tau_pt_threshold": 20,  
            "tau_eta_threshold": 2.1, 
            "tau_dz_threshold": 0.2, 
            "tau_vs_jet": "Loose",
            "tau_vs_ele": "VVLoose", 
            "tau_vs_mu": "Loose",
            "prongs": 13,            
        },        
        "ele": {
            "tau_pt_threshold": 20,  
            "tau_eta_threshold": 2.1, 
            "tau_dz_threshold": 0.2, 
            "tau_vs_jet": "Loose",
            "tau_vs_ele": "VVLoose", 
            "tau_vs_mu": "Loose",
            "prongs": 13,            
        },        
    },
    "1b1e1mu": {
        "mu": {
            "tau_pt_threshold": 20,  
            "tau_eta_threshold": 2.1, 
            "tau_dz_threshold": 0.2, 
            "tau_vs_jet": "Loose",
            "tau_vs_ele": "VVLoose", 
            "tau_vs_mu": "Loose",
            "prongs": 13,            
        },        
        "ele": {
            "tau_pt_threshold": 20,  
            "tau_eta_threshold": 2.1, 
            "tau_dz_threshold": 0.2, 
            "tau_vs_jet": "Loose",
            "tau_vs_ele": "VVLoose", 
            "tau_vs_mu": "Loose",
            "prongs": 13,            
        },        
    },
    "1b1l": {
        "mu": {
            "tau_pt_threshold": 20,  
            "tau_eta_threshold": 2.1, 
            "tau_dz_threshold": 0.2, 
            "tau_vs_jet": "Loose",
            "tau_vs_ele": "VVLoose", 
            "tau_vs_mu": "Loose",
            "prongs": 13,            
        },        
        "ele": {
            "tau_pt_threshold": 20,  
            "tau_eta_threshold": 2.1, 
            "tau_dz_threshold": 0.2, 
            "tau_vs_jet": "Loose",
            "tau_vs_ele": "VVLoose", 
            "tau_vs_mu": "Loose",
            "prongs": 13,            
        },        
    },
}

ttbar_jet_selection = {
    "1b1e1mu": {
        "ele": {
            "jet_pt_threshold": 20,
            "btag_working_point": "M",
            "jet_id": 6,
            "jet_pileup_id": "T"
        },
        "mu": {
            "jet_pt_threshold": 20,
            "btag_working_point": "M",
            "jet_id": 6,
            "jet_pileup_id": "T"
        },
    },
    "2b1l": {
        "ele": {
            "jet_pt_threshold": 20,
            "btag_working_point": "M",
            "jet_id": 6,
            "jet_pileup_id": "T"
        },
        "mu": {
            "jet_pt_threshold": 20,
            "btag_working_point": "M",
            "jet_id": 6,
            "jet_pileup_id": "T"
        },
    },
    "1b1l": {
        "ele": {
            "jet_pt_threshold": 20,
            "btag_working_point": "M",
            "jet_id": 6,
            "jet_pileup_id": "T"
        },
        "mu": {
            "jet_pt_threshold": 20,
            "btag_working_point": "M",
            "jet_id": 6,
            "jet_pileup_id": "T"
        },
    },
}