import json
import glob
import argparse
import subprocess
from pathlib import Path


def main(args):
    """Helper function to resubmit condor jobs"""

    main_dir = Path.cwd()
    fileset_path = Path(f"{main_dir}/wprime_plus_b/fileset")
    with open(f"{fileset_path}/das_datasets.json", "r") as f:
        datasets = json.load(f)[f"{args.year}_UL"]

    # get jobs done
    jobs_done = []
    out_path = f"{args.output_path}/{args.processor}"
    if args.label:
        out_path += f"/{args.label}"
    if args.processor not in ["btag_eff", "trigger_eff", "susy"]:
        out_path += f"/{args.channel}/{args.lepton_flavor}"
    if args.processor == "trigger_eff":
        out_path += f"/{args.lepton_flavor}"
    out_path += f"/{args.year}"
    print(f"{out_path=}")
    
    for sample in datasets:
        output_list = glob.glob(f"{out_path}/{sample}*.pkl")
        for output in output_list:
            jobs_done.append(output.split("/")[-1].replace(".pkl", ""))
    n_jobs_done = len(jobs_done)
    
    # get jobs to be run
    condor_path = f"{main_dir}/condor/{args.processor}"
    if args.processor not in ["btag_eff", "trigger_eff", "susy"]:
        condor_path += f"/{args.channel}/{args.lepton_flavor}"
    if args.processor == "trigger_eff":
        condor_path += f"/{args.lepton_flavor}"
    condor_path += f"/{args.year}"
    print(f"{condor_path=}")
    condor_files = glob.glob(f"{condor_path}/*/*.sub", recursive=True)
    n_jobs = len(condor_files)
    
    # get jobs names
    to_replace = f"{args.processor}_"
    if args.processor not in ["btag_eff", "trigger_eff", "susy"]:
        to_replace += f"{args.channel}_{args.lepton_flavor}_"
    if args.processor == "trigger_eff":
        to_replace += f"{args.lepton_flavor}_"
    condor_files_keys = [
        f.split("/")[-1]
        .replace(to_replace, "")
        .replace(".sub", "")
        for f in condor_files
    ]

    for job, sub_file in zip(condor_files_keys, condor_files):
        if job not in jobs_done:
            print(job)
            if args.resubmit == "True":
                # resubmit missing jobs
                subprocess.run(["condor_submit", sub_file])

    print("")
    print(f"{n_jobs=}")
    print(f"{n_jobs_done=}")
    print(f"missing files: {n_jobs - n_jobs_done}", "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_path",
        dest="output_path",
        type=str,
        #default="/eos/user/../../susy_vbf/outs"
        help="path to the 'outs' folder",
    )
    parser.add_argument(
        "--label",
        dest="label",
        type=str,
        default="",
        help="label of the run",
    )
    parser.add_argument(
        "--processor",
        dest="processor",
        type=str,
        default="susy",
        help="processor to be used",
    )
    parser.add_argument(
        "--channel",
        dest="channel",
        type=str,
        default="",
        help="channel to be processed {'2b1l', '1b1e1mu', '1b1l'}",
    )
    parser.add_argument(
        "--lepton_flavor",
        dest="lepton_flavor",
        type=str,
        default="",
        help="lepton flavor to be processed {'mu', 'ele'}",
    )
    parser.add_argument(
        "--year",
        dest="year",
        type=str,
        default="2017",
        help="year of the data {2016APV, 2016, 2017, 2018} (default 2017)",
    )
    parser.add_argument(
        "--resubmit",
        dest="resubmit",
        type=str,
        default="False",
        help="if True resubmit the jobs. if False only print the missing jobs",
    )
    args = parser.parse_args()
    main(args)