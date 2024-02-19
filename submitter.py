import json
import time
import dask
import pickle
import argparse
import datetime
import numpy as np
import wprime_plus_b.utils
from pathlib import Path
from coffea import processor
from dask.distributed import Client
from humanfriendly import format_timespan
from utils import run_checker, manage_args, get_filesets
from distributed.diagnostics.plugin import UploadDirectory
from wprime_plus_b.utils import paths
from wprime_plus_b.utils.load_config import load_processor_config, load_dataset_config
from wprime_plus_b.processors.trigger_efficiency_processor import (
    TriggerEfficiencyProcessor,
)
from wprime_plus_b.processors.btag_efficiency_processor import BTagEfficiencyProcessor
from wprime_plus_b.processors.ttbar_analysis import TtbarAnalysis
from wprime_plus_b.processors.ztoll_processor import ZToLLProcessor
from wprime_plus_b.processors.qcd_analysis import QcdAnalysis
from wprime_plus_b.selections.ttbar.config import (
    ttbar_electron_selection,
    ttbar_muon_selection,
    ttbar_jet_selection,
)
from wprime_plus_b.selections.ztoll.config import (
    ztoll_electron_selection,
    ztoll_muon_selection,
    ztoll_jet_selection,
)
from wprime_plus_b.selections.qcd.config import (
    qcd_electron_selection,
    qcd_muon_selection,
    qcd_jet_selection,
)


def main(args):
    # check and manage args
    run_checker(vars(args))
    args = manage_args(vars(args))
    
    # load dataset config
    dataset_config = load_dataset_config(config_name=args["sample"])
    
    # divide dataset into 'nsplit' json files
    filesets = get_filesets(
        sample=args["sample"],
        year=args["year"] + args["yearmod"],
        nsplit=dataset_config.nsplit,
    )
    # run over each sample partition
    for sample, fileset_path in filesets.items():
        if len(args['nsample']) != 0:
            if sample.split("_")[-1] not in args['nsample']:
                continue
        print(f"Processing {sample}")
        fileset = {}
        with open(fileset_path, "r") as handle:
            data = json.load(handle)
        for root_file in data.values():
            if args["nfiles"] != -1:
                root_file = root_file[: args["nfiles"]]
        fileset[sample] = [f"root://{args['redirector']}/" + file for file in root_file]

        # define processors and its kwargs
        processors = {
            "ttbar": TtbarAnalysis,
            "ztoll": ZToLLProcessor,
            "qcd": QcdAnalysis,
            "btag_eff": BTagEfficiencyProcessor,
            "trigger_eff": TriggerEfficiencyProcessor,
        }
        processor_args = [
            "year",
            "yearmod",
            "channel",
            "lepton_flavor",
            "output_type",
            "syst",
        ]
        processor_kwargs = {k: args[k] for k in processor_args if args[k]}

        # define executors
        executors = {
            "iterative": processor.iterative_executor,
            "futures": processor.futures_executor,
            "dask": processor.dask_executor,
        }
        executor_args = {
            "schema": processor.NanoAODSchema,
        }
        if args["executor"] == "futures":
            executor_args.update({"workers": args["workers"]})
        if args["executor"] == "dask":
            client = Client("tls://localhost:8786")
            executor_args.update({"client": client})
            # upload local directory to dask workers
            try:
                client.register_worker_plugin(
                    UploadDirectory(f"{Path.cwd()}", restart=True, update_path=True),
                    nanny=True,
                )
                print(f"Uploaded {Path.cwd()} succesfully")
            except OSError:
                print("Failed to upload the directory")
                
        # get processor config
        processor_config_name = "_".join(
            [
                i
                for i in [args["processor"], args["channel"], args["lepton_flavor"]]
                if i
            ]
        )
        processor_config = load_processor_config(config_name=processor_config_name)
        processor_output_path = paths.processor_path(
            processor_name=processor_config.name,
            processor_lepton_flavour=processor_config.lepton_flavor,
            processor_channel=processor_config.channel,
            dataset_year=args["year"] + args["yearmod"],
            mkdir=True,
        )
        
        # run processor
        t0 = time.monotonic()
        out = processor.run_uproot_job(
            fileset,
            treename="Events",
            processor_instance=processors[args["processor"]](**processor_kwargs),
            executor=executors[args["executor"]],
            executor_args=executor_args,
        )
        exec_time = format_timespan(time.monotonic() - t0)

        # get metadata
        metadata = {"walltime": exec_time}
        metadata.update({"fileset": fileset[sample]})
        if "metadata" in out[sample]:
            output_metadata = out[sample]["metadata"]
            # save number of raw initial events
            metadata.update(
                {"events_before": float(output_metadata["events_before"])}
            )
            # save number of weighted initial events
            metadata.update({"sumw": float(output_metadata["sumw"])})

            if args["processor"] in ["qcd"]:
                metadata.update({"nevents": {}})
                for region in ["A", "B", "C", "D"]:
                    metadata["nevents"].update({region: {}})
                    metadata["nevents"][region]["events_after"] = str(
                        output_metadata[region]["events_after"]
                    )
                    metadata["nevents"][region]["events_after_weighted"] = str(
                        output_metadata[region]["events_after_weighted"]
                    )

            if args["processor"] in ["ttbar", "ztoll"]:
                metadata.update(
                    {"events_after": float(output_metadata["events_after"])}
                )
                # save cutflow to metadata
                for cut_selection, nevents in output_metadata["cutflow"].items():
                    output_metadata["cutflow"][cut_selection] = str(nevents)
                metadata.update({"cutflow": output_metadata["cutflow"]})

                for weight, statistics in output_metadata[
                    "weight_statistics"
                ].items():
                    output_metadata["weight_statistics"][weight] = str(statistics)
                metadata.update(
                    {"weight_statistics": output_metadata["weight_statistics"]}
                )

            if args["processor"] in ["ttbar", "ztoll", "qcd"]:
                # save selectios to metadata
                selections = {
                    "ttbar": {
                        "electron_selection": ttbar_electron_selection[
                            args["channel"]
                        ][args["lepton_flavor"]],
                        "muon_selection": ttbar_muon_selection[args["channel"]][
                            args["lepton_flavor"]
                        ],
                        "jet_selection": ttbar_jet_selection[args["channel"]][
                            args["lepton_flavor"]
                        ],
                    },
                    "ztoll": {
                        "electron_selection": ztoll_electron_selection,
                        "muon_selection": ztoll_muon_selection,
                        "jet_selection": ztoll_jet_selection,
                    },
                    "qcd": {
                        "electron_selection": qcd_electron_selection,
                        "muon_selection": qcd_muon_selection,
                        "jet_selection": qcd_jet_selection,
                    },
                }
                metadata.update(
                    {
                        "electron_selection": selections[args["processor"]][
                            "electron_selection"
                        ]
                    }
                )
                metadata.update(
                    {
                        "muon_selection": selections[args["processor"]][
                            "muon_selection"
                        ]
                    }
                )
                metadata.update(
                    {
                        "jet_selection": selections[args["processor"]][
                            "jet_selection"
                        ]
                    }
                )
            # save args to metadata
            args_dict = args.copy()
            metadata.update(args_dict)

            # save metadata
            metadata_path = Path(f"{processor_output_path}/metadata")
            if not metadata_path.exists():
                metadata_path.mkdir(parents=True)
            with open(f"{metadata_path}/{sample}_metadata.json", "w") as f:
                f.write(json.dumps(metadata))
                
        # save output
        with open(f"{processor_output_path}/{sample}.pkl", "wb") as handle:
            pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--redirector",
        dest="redirector",
        type=str,
        default="xcache",
        help="redirector to find CMS datasets {use 'xcache' at coffea-casa. use 'cmsxrootd.fnal.gov', 'xrootd-cms.infn.it' or 'cms-xrd-global.cern.ch' at lxplus} (default xcache)",
    )
    parser.add_argument(
        "--processor",
        dest="processor",
        type=str,
        default="ttbar",
        help="processor to be used {trigger, ttbar, candle, btag_eff} (default ttbar)",
    )
    parser.add_argument(
        "--executor",
        dest="executor",
        type=str,
        default="iterative",
        help="executor to be used {iterative, futures, dask} (default iterative)",
    )
    parser.add_argument(
        "--workers",
        dest="workers",
        type=int,
        default=4,
        help="number of workers to use with futures executor (default 4)",
    )
    parser.add_argument(
        "--year",
        dest="year",
        type=str,
        default="2017",
        help="year of the data {2016, 2017, 2018} (default 2017)",
    )
    parser.add_argument(
        "--yearmod",
        dest="yearmod",
        type=str,
        default="",
        help="year modifier {'', 'APV'} (default '')",
    )
    parser.add_argument(
        "--channel",
        dest="channel",
        type=str,
        default="",
        help="channel to be processed {'2b1l', '1b1e1mu'}",
    )
    parser.add_argument(
        "--lepton_flavor",
        dest="lepton_flavor",
        type=str,
        default="mu",
        help="lepton flavor to be processed {'mu', 'ele'}",
    )
    parser.add_argument(
        "--nfiles",
        dest="nfiles",
        type=int,
        default=1,
        help="number of .root files to be processed by sample. To run all files use -1 (default 1)",
    )
    parser.add_argument(
        "--nsample",
        dest="nsample",
        type=list,
        default=[],
        help="nsample",
    )
    parser.add_argument(
        "--chunksize",
        dest="chunksize",
        type=int,
        default=50000,
        help="number of chunks to process",
    )
    parser.add_argument(
        "--output_type",
        dest="output_type",
        type=str,
        default="hist",
        help="type of output {hist, array}",
    )
    parser.add_argument(
        "--syst",
        dest="syst",
        type=str,
        default="nominal",
        help="systematic to apply {'nominal', 'jet', 'met', 'full'}",
    )
    parser.add_argument(
        "--sample",
        dest="sample",
        type=str,
        default="all",
        help="sample key to be processed {'all', 'mc' or <sample_name>} (default all)",
    )
    args = parser.parse_args()
    main(args)