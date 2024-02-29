import json
import glob
import subprocess
from pathlib import Path
from collections import OrderedDict


def get_filesets(sample: str, year: str, facility: str):
    """return a dictionary with sample names as keys and .json files as values"""
    main_dir = Path.cwd()
    fileset_path = Path(f"{main_dir}/wprime_plus_b/fileset/{facility}")
    file_list = glob.glob(f"{fileset_path}/*.json")
    filesets = {}
    for file in file_list:
        file_name = file.split("/")[-1].replace(".json", "")
        if file_name.startswith(sample):
            filesets[file_name] = file
    if len(filesets) != 1:
        # sort the dictionary keys based on the number after the "_" in ascending order
        sorted_keys = sorted(filesets.keys(), key=lambda x: int(x.split("_")[1]))
        # create an ordered dictionary using the sorted keys
        ordered_filesets = OrderedDict((key, filesets[key]) for key in sorted_keys)
        return ordered_filesets
    return filesets


def run_checker(args: dict) -> None:
    # check that datasets exists
    fileset_path = Path(f"{Path.cwd()}/wprime_plus_b/fileset/{args['facility']}")
    assert (
        fileset_path.exists()
    ), f"Filesets have not been created. Type 'python3 process_filesets.py --facility {args['facility']}'"
    
    # check processor
    available_processors = ["ttbar", "ztoll", "qcd", "btag_eff", "trigger_eff"]
    assert (
        args["processor"] in available_processors
    ), f"Incorrect processor. Available processors are: {available_processors}"

    if args["processor"] == "ttbar":
        # check channel
        available_channels = ["2b1l", "1b1e1mu", "1b1l"]
        assert (
            args["channel"] in available_channels
        ), f"Incorrect channel. Available channels are: {available_channels}"
        # check lepton flavor
        available_lepton_flavors = ["ele", "mu"]
        assert (
            args["lepton_flavor"] in available_lepton_flavors
        ), f"Incorrect lepton flavor. Available lepton flavors are: {available_lepton_flavors}"

        # check Data sample
        if args["channel"] == "1b1e1mu":
            if args["lepton_flavor"] == "mu":
                assert (
                    args["sample"] != "SingleMuon"
                ), "1b1e1mu muon channel should be run with SingleElectron dataset"
            else:
                assert (
                    args["sample"] != "SingleElectron"
                ), "1b1e1mu electron channel should be run with SingleMuon dataset"
        else:
            if args["lepton_flavor"] == "mu":
                assert (
                    args["sample"] != "SingleElectron"
                ), "2b1l muon channel should be run with SingleMuon dataset"
            else:
                assert (
                    args["sample"] != "SingleMuon"
                ), "2b1l electron channel should be run with SingleElectron dataset"
    if args["processor"] == "qcd":
        # check channel
        available_channels = ["A", "B", "C", "D"]
        assert (
            args["channel"] in available_channels
        ), f"Incorrect channel. Available channels are: {available_channels}"
        assert (
            args["lepton_flavor"] == "mu" and args["output_type"] == "hist"
        ), "Only muon channel and histograms are available"


def manage_args(args: dict) -> dict:
    processor_args_mapping = {
        "ztoll": ["channel", "syst"],
        "qcd": ["syst"],
        "btag_eff": ["lepton_flavor", "channel", "syst"],
        "trigger_eff": ["channel", "syst", "output_type"],
    }
    processor = args.get("processor")
    if processor in processor_args_mapping:
        for arg in processor_args_mapping[processor]:
            args[arg] = None
    return args