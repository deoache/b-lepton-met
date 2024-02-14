import json
from pathlib import Path

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


def divide_fileset(sample: str, year: str, nsplit: int):
    # load fileset  
    main_dir = Path.cwd()
    fileset_path = Path(f"{main_dir}/wprime_plus_b/fileset")
    with open(f"{fileset_path}/fileset_{year}_UL_NANO.json", "r") as handle:
        data = json.load(handle)
        
    # make output filesets directory
    output_directory = Path(f"{fileset_path}/filesets")
    if output_directory.exists():
        for file in output_directory.glob("*"):
            if file.is_file():
                file.unlink()
    else:
        output_directory.mkdir(parents=True)
        
    # split fileset and save filesets
    filesets = {}
    if nsplit == 1:
        filesets[sample] = f"{output_directory}/{sample}.json"
        sample_data = {sample: data[sample]}
        with open(f"{output_directory}/{sample}.json", "w") as json_file:
            json.dump(sample_data, json_file, indent=4, sort_keys=True)
    else:
        root_files_list = divide_list(data[sample], nsplit) 
        keys = ".".join(f"{sample}_{i}" for i in range(1, nsplit + 1)).split(".")
        for key, value in zip(keys, root_files_list):
            sample_data = {}
            sample_data[key] = list(value)

            filesets[key] = f"{output_directory}/{key}.json"
            with open(f"{output_directory}/{key}.json", "w") as json_file:
                json.dump(sample_data, json_file, indent=4, sort_keys=True)

    return filesets

def get_filesets(sample: str, year: str, nsplit: int):
    return divide_fileset(sample, year, nsplit)

def run_checker(args: dict) -> None:
    if args["processor"] == "ttbar":
        assert args[
            "channel"
        ], "No 'channel' argument specified! You must to provide a channel"
        assert args["channel"] in [
            "2b1l",
            "1b1e1mu",
            "1b1l",
        ], "Incorrect 'channel' argument. Options are '2b1l', '1b1e1mu', '1b1l'"
        assert args["lepton_flavor"] in [
            "ele",
            "mu",
        ], "Incorrect 'lepton_flavor' argument. Options are 'ele' or 'mu'"
        
    if args["processor"] == "qcd":
        assert (
            args["lepton_flavor"] == "mu" and args["output_type"] == "hist"
        ), "Only muon channel and histograms are available"


def manage_args(args: dict) -> dict:
    processor_args_mapping = {
        "ztoll": ["channel", "syst"],
        "qcd": ["channel", "syst"],
        "btag_eff": ["lepton_flavor", "channel", "syst"],
        "trigger_eff": ["channel", "syst", "output_type"],
    }
    processor = args.get("processor")
    if processor in processor_args_mapping:
        for arg in processor_args_mapping[processor]:
            args[arg] = None
    return args