import os
import json
from coffea.dataset_tools.dataset_query import DataDiscoveryCLI


ERAS = {
    "2016APV": ["B", "C", "D", "F"],
    "2016": ["F", "G", "H"],
    "2017": ["B", "C", "D", "E", "F"],
    "2018": ["A", "B", "C", "D"],
}

SITES = {
    "2016APV": [
        "T3_US_FNALLPC",
        "T1_US_FNAL_Disk",
        "T2_US_Vanderbilt",
        "T2_US_Purdue",
        "T2_US_Nebraska",
        "T2_DE_DESY",
        "T2_BE_IIHE",
        "T2_CH_CERN",
        "T1_DE_KIT_Disk",
        "T2_DE_RWTH",
    ],
    "2016": [
        "T3_US_FNALLPC",
        "T1_US_FNAL_Disk",
        "T2_US_Vanderbilt",
        "T2_US_Purdue",
        "T2_US_Nebraska",
        "T2_DE_DESY",
        "T2_BE_IIHE",
        "T2_CH_CERN",
        "T1_DE_KIT_Disk",
        "T2_DE_RWTH",
    ],
    "2017": [
        "T3_US_FNALLPC",
        "T1_US_FNAL_Disk",
        "T2_US_Vanderbilt",
        "T2_US_Purdue",
        "T2_US_Nebraska",
        "T2_DE_DESY",
        "T2_BE_IIHE",
        "T2_CH_CERN",
    ],
    "2018": [
        "T3_US_FNALLPC",
        "T1_US_FNAL_Disk",
        "T2_US_Vanderbilt",
        "T2_US_Purdue",
        "T2_US_Nebraska",
        "T2_DE_DESY",
        "T2_BE_IIHE",
        "T2_CH_CERN",
        "T1_DE_KIT_Disk",
        "T2_DE_RWTH",
    ],
}


def main():
    with open("das_datasets.json", "r") as f:
        datasets = json.load(f)
    for year in ERAS.keys():
        # create a dataset_definition dict for each year
        yreco = f"{year}_UL"
        if not datasets[yreco]:
            continue
        dataset_definition = {}
        for dataset_key, dataset in datasets[yreco].items():
            if isinstance(dataset, list):
                for _dataset, era in zip(dataset, ERAS[year]):
                    dataset_definition[f"/{_dataset}"] = {
                        "short_name": f"{dataset_key}_{era}",
                        "metadata": {"isMC": True},
                    }
            else:
                dataset_definition[f"/{dataset}"] = {
                    "short_name": dataset_key,
                    "metadata": {"isMC": False},
                }
        # the dataset definition is passed to a DataDiscoveryCLI
        ddc = DataDiscoveryCLI()
        # set the allow sites to look for replicas
        ddc.do_allowlist_sites(SITES[year])
        # query rucio and get replicas
        ddc.load_dataset_definition(
            dataset_definition,
            query_results_strategy="all",
            replicas_strategy="round-robin",
        )
        ddc.do_save(f"dataset_discovery_{yreco}.json")
        
        # load and reformat generated fileset
        with open(f"dataset_discovery_{yreco}.json", "r") as f:
            dataset_discovery = json.load(f)
        new_dataset = {key: [] for key in datasets[yreco]}
        for dataset in dataset_discovery:
            root_files = list(dataset_discovery[dataset]["files"].keys())
            dataset_key = dataset_discovery[dataset]["metadata"]["short_name"]
            if dataset_key.startswith("Single"):
                new_dataset[dataset_key.split("_")[0]] += root_files
            else:
                new_dataset[dataset_key] = root_files
        # save new fileset and drop 'dataset_discovery' fileset
        os.remove(f"dataset_discovery_{yreco}.json")
        with open(f"fileset_{yreco}_NANO_lxplus.json", "w") as json_file:
            json.dump(new_dataset, json_file, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()