import json
import argparse
from pathlib import Path
from wprime_plus_b.utils.load_config import load_dataset_config


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

def main(args):
    main_dir = Path.cwd()
    fileset_path = Path(f"{main_dir}/wprime_plus_b/fileset")
    
    # make output filesets directory
    output_directory = Path(f"{fileset_path}/{args.facility}")
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
            if args.facility == "lxplus":
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
            
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--facility",
        dest="facility",
        type=str,
        default="coffea-casa",
        help="facility where jobs will be launched {coffea-casa, lxplus}",
    )
    args = parser.parse_args()
    main(args)