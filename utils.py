import os
import json
import glob
import subprocess
from pathlib import Path
from collections import OrderedDict
from wprime_plus_b.utils.load_config import load_dataset_config


def get_command(args):
    """return command to submit jobs at coffea-casa or lxplus"""
    cmd = f"python submit.py"
    for arg in args:
        if args[arg]:
            cmd += f" --{arg} {args[arg]}"
    cmd += f" --username {os.environ['USER']}"
    return cmd

def divide_list(lst: list, n: int):
    """
    Divide a list into n sublists
    """
    size = len(lst) // n
    remainder = len(lst) % n
    result = []
    start = 0
    for i in range(n):
        if i < remainder:
            end = start + size + 1
        else:
            end = start + size
        result.append(lst[start:end])
        start = end
    return result

def build_filesets(facility: str):
    main_dir = Path.cwd()
    fileset_path = Path(f"{main_dir}/wprime_plus_b/fileset")
    
    # make output filesets directory
    output_directory = Path(f"{fileset_path}/{facility}")
    if output_directory.exists():
        for file in output_directory.glob("*"):
            if file.is_file():
                file.unlink()
    else:
        output_directory.mkdir(parents=True)
    
    # load filesets
    with open(f"{fileset_path}/das_datasets.json", "r") as f:
        datasets = json.load(f)
    for yreco in datasets:
        year = yreco.replace("_UL", "")
        for sample in datasets[yreco]:
            if facility == "lxplus":
                json_file = f"{fileset_path}/fileset_{year}_UL_NANO_lxplus.json"
            else:
                json_file = f"{fileset_path}/fileset_{year}_UL_NANO.json"
            with open(json_file, "r") as handle:
                data = json.load(handle)

            # split fileset and save filesets
            filesets = {}
            # load dataset config
            dataset_config = load_dataset_config(config_name=sample)
            if dataset_config.nsplit == 1:
                filesets[sample] = f"{output_directory}/{sample}.json"
                sample_data = {sample: data[sample]}
                with open(f"{output_directory}/{sample}.json", "w") as json_file:
                    json.dump(sample_data, json_file, indent=4, sort_keys=True)
            else:
                root_files_list = divide_list(data[sample], dataset_config.nsplit)
                keys = ".".join(f"{sample}_{i}" for i in range(1, dataset_config.nsplit + 1)).split(".")
                for key, value in zip(keys, root_files_list):
                    sample_data = {}
                    sample_data[key] = list(value)

                    filesets[key] = f"{output_directory}/{key}.json"
                    with open(f"{output_directory}/{key}.json", "w") as json_file:
                        json.dump(sample_data, json_file, indent=4, sort_keys=True)

                        
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
        sorted_keys = sorted(filesets.keys(), key=lambda x: int(x.split("_")[-1]))
        # create an ordered dictionary using the sorted keys
        ordered_filesets = OrderedDict((key, filesets[key]) for key in sorted_keys)
        return ordered_filesets
    return filesets


def run_checker(args: dict) -> None:
    # check args
    requiered_args = ["processor", "sample", "executor", "year", "facility", "output_type", "nfiles"]
    for arg in requiered_args:
        assert arg in args, f"You must provide the {arg} argument!"
            
    # check facility
    available_facilitys = ["coffea-casa", "lxplus"]
    assert (
        args["facility"] in available_facilitys
    ), f"Incorrect facility. Available facilities are: {available_facilitys}"
    
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
        
    # check executor
    available_executors = ["iterative", "futures"]
    assert (
        args["executor"] in available_executors
    ), f"Incorrect executor. Available executors are: {available_executors}"
    
    # check years
    available_years = ["2016", "2017", "2018"]
    assert (
        args["year"] in available_years
    ), f"Incorrect year. Available years are: {available_years}"
    
    available_yearmods = ["APV", ""]
    assert (
        args["yearmod"] in available_yearmods
    ), f"Incorrect yearmod. Available yearmods are: {available_yearmods}"
    
    # check sample
    fileset_path = Path(f"{Path.cwd()}/wprime_plus_b/fileset")
    with open(f"{fileset_path}/das_datasets.json", "r") as f:
        datasets = json.load(f)[args["year"] + args["yearmod"] + "_UL"]
    available_samples = list(datasets.keys())
    assert (
        args["sample"] in available_samples
    ), f"Incorrect sample. Available samples are: {available_samples}"
    
    # check nsample
    dataset_config = load_dataset_config(config_name=args["sample"])
    available_nsamples = [""] + [str(i) for i in range(1, dataset_config.nsplit + 1)]
    nsamples = args["nsample"].split(",")
    for nsample in nsamples:
        assert (
            nsample in available_nsamples
        ), f"Incorrect nsample. Available nsamples are: {available_nsamples}"
    
    # check output type
    available_output_types = ["hist", "array"]
    assert (
        args["output_type"] in available_output_types
    ), f"Incorrect output_type. Available output_types are: {available_output_types}"
    
    # check systematics
    available_systs = ['nominal', 'jet', 'met', 'full']
    assert (
        args["syst"] in available_systs
    ), f"Incorrect syst. Available systs are: {available_systs}"
    

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