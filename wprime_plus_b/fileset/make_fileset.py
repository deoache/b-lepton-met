import os
import re
import json
import subprocess

datasets = {
    "2016APV_UL": {},
    "2016_UL": {},
    "2017_UL": {
        # Drell-Yan
        "DYJetsToLL_M-50_HT-70to100": "DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-100to200": "DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-1200to2500": "DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-200to400": "DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-2500toInf": "DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-400to600": "DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-600to800": "DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "DYJetsToLL_M-50_HT-800to1200": "DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        # Single Top
        "ST_s-channel_4f_leptonDecays": "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "ST_t-channel_antitop_5f_InclusiveDecays": "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",
        "ST_t-channel_top_5f_InclusiveDecays": "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "ST_tW_antitop_5f_inclusiveDecays": "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",
        "ST_tW_top_5f_inclusiveDecays": "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",
        # W+Jets
        "WJetsToLNu_HT-100To200": "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "WJetsToLNu_HT-200To400": "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "WJetsToLNu_HT-400To600": "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "WJetsToLNu_HT-600To800": "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "WJetsToLNu_HT-800To1200": "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v3/NANOAODSIM",
        "WJetsToLNu_HT-1200To2500": "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "WJetsToLNu_HT-2500ToInf": "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",
        # TT
        # no xsec
        "TTTo2L2Nu": "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "TTToHadronic": "TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "TTToSemiLeptonic": "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        # Diboson
        "WW": "WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "WZ": "WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "ZZ": "ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        # Higgs
        "VBFHToWWTo2L2Nu": "VBFHToWWTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV714_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",
        "VBFHToWWToLNuQQ": "VBFHToWWToLNuQQ_M-125_TuneCP5_withDipoleRecoil_13TeV_powheg_jhugen751_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        "GluGluHToWWToLNuQQ": "GluGluHToWWToLNuQQ_M-125_TuneCP5_13TeV_powheg_jhugen751_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",
        # Data
        "SingleMuon": [
            "SingleMuon/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleMuon/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleMuon/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleMuon/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleMuon/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
        ],
        "SingleElectron": [
            "SingleElectron/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleElectron/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleElectron/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleElectron/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
            "SingleElectron/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD",
        ],
    },
    "2018_UL": {},
}

if __name__ == "__main__":
    reco = "UL"
    filesets = {}
    years = ["2016APV", "2016", "2017", "2018"]
    for year in years:
        yreco = f"{year}_{reco}"
        filesets[yreco] = {}
        for dataset_key, dataset in datasets[yreco].items():
            querys = []
            if isinstance(dataset, list):
                for x in dataset:
                    querys.append(f"file dataset=/{x}")
            else:
                querys.append(f"file dataset=/{dataset}")
            filearray = []
            for newquery in querys:
                farray = subprocess.run(
                    ["dasgoclient", f"-query={newquery}"],
                    stdout=subprocess.PIPE,
                    universal_newlines=True,
                )
                stdout = farray.stdout
                stdout_array = stdout.split("\n")
                stdout_array = stdout_array[:-1]
                stdout_array[-1] = stdout_array[-1].replace(",", "")
                filearray.extend(stdout_array)
            filesets[yreco][dataset_key] = filearray
        with open(f"fileset_{yreco}_NANO.json", "w") as json_file:
            json.dump(filesets[yreco], json_file, indent=4, sort_keys=True)