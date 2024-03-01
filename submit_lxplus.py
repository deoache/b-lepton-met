import os
import argparse
import subprocess
from pathlib import Path
from utils import get_command, run_checker, build_filesets
from wprime_plus_b.utils.load_config import load_dataset_config


def move_X509():
    """move x509 proxy file from /tmp to /afs/private. Returns the afs path"""
    try:
        x509_localpath = (
            [
                line
                for line in os.popen("voms-proxy-info").read().split("\n")
                if line.startswith("path")
            ][0]
            .split(":")[-1]
            .strip()
        )
    except Exception as err:
        raise RuntimeError(
            "x509 proxy could not be parsed, try creating it with 'voms-proxy-init --voms cms'"
        ) from err
    x509_path = f"{Path.home()}/private/{x509_localpath.split('/')[-1]}"
    subprocess.run(["cp", x509_localpath, x509_path])
    return x509_path


def submit_condor(jobname:str , cmd: str, flavor: str):
    """build condor and executable files, and submit condor job"""
    main_dir = Path.cwd()
    condor_dir = Path(f"{main_dir}/condor")

    # create logs and output directories
    log_dir = Path(str(condor_dir) + "/logs")
    if not log_dir.exists():
        log_dir.mkdir()
    # make condor file
    condor_template_file = open(f"{condor_dir}/submit.sub")
    local_condor = f"{condor_dir}/{jobname}.sub"
    condor_file = open(local_condor, "w")
    for line in condor_template_file:
        line = line.replace("DIRECTORY", str(condor_dir))
        line = line.replace("JOBNAME", jobname)
        line = line.replace("JOBFLAVOR", f'"{flavor}"')
        condor_file.write(line)
    condor_file.close()
    condor_template_file.close()

    # make executable file
    x509_path = move_X509()
    sh_template_file = open(f"{condor_dir}/submit.sh")
    local_sh = f"{condor_dir}/{jobname}.sh"
    sh_file = open(local_sh, "w")
    for line in sh_template_file:
        line = line.replace("MAINDIRECTORY", str(main_dir))
        line = line.replace("COMMAND", cmd)
        line = line.replace("X509PATH", x509_path)
        sh_file.write(line)
    sh_file.close()
    sh_template_file.close()

    # submit jobs
    print(f"submitting {jobname}")
    subprocess.run(["condor_submit", local_condor])


def main(args):
    build_filesets(facility="lxplus")
    args = vars(args)
    run_checker(args)
    jobname = f'{args["processor"]}_{args["channel"]}_{args["lepton_flavor"]}_{args["sample"]}'
    dataset_config = load_dataset_config(config_name=args["sample"])
    if dataset_config.nsplit == 1:
        cmd = get_command(args)
        submit_condor(jobname, cmd, flavor="microcentury")
    else:
        for nsplit in range(1, dataset_config.nsplit + 1):
            args["nsample"] = nsplit
            cmd = get_command(args)
            submit_condor(jobname + f"_{nsplit}", cmd, flavor="longlunch")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--processor",
        dest="processor",
        type=str,
        default="ttbar",
        help="processor to be used {ttbar, ztoll, qcd, trigger_eff, btag_eff} (default ttbar)",
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
        default="mu",
        help="lepton flavor to be processed {'mu', 'ele'}",
    )
    parser.add_argument(
        "--sample",
        dest="sample",
        type=str,
        default="all",
        help="sample key to be processed",
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
        "--nfiles",
        dest="nfiles",
        type=int,
        default=1,
        help="number of .root files to be processed by sample. To run all files use -1 (default 1)",
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
        "--facility",
        dest="facility",
        type=str,
        default="lxplus",
        help="facility to launch jobs {coffea-casa, lxplus}",
    )
    parser.add_argument(
        "--nsample",
        dest="nsample",
        type=str,
        default="",
        help="partitions to run (--nsample 1,2,3 will only run partitions 1,2 and 3)",
    )
    args = parser.parse_args()
    main(args)