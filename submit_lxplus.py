import os
import argparse
import subprocess
from pathlib import Path
from wprime_plus_b.utils.load_config import load_dataset_config


def get_command(args, nsample=None):
    cmd = f"python submit.py"
    cmd += f" --processor {args.processor}"
    cmd += f" --channel {args.channel}"
    cmd += f" --lepton_flavor {args.lepton_flavor}"
    cmd += f" --year {args.year}"
    cmd += f" --sample {args.sample}"
    cmd += f" --executor {args.executor}"
    cmd += f" --workers {args.workers}"
    cmd += f" --nfiles {args.nfiles}"
    cmd += f" --output_type {args.output_type}"
    cmd += f" --syst {args.syst}"
    cmd += f" --facility lxplus"
    if nsample:
        cmd += f" --nsample {nsample}"
    return cmd


def submit_condor(jobname, cmd):
    main_dir = Path.cwd()
    condor_dir = Path(f"{main_dir}/condor")

    # create logs and output directories
    log_dir = Path(str(condor_dir) + "/logs")
    if not log_dir.exists():
        log_dir.mkdir()
        
    # define EOS area to move output files
    username = os.environ["USER"]
    eos_dir = Path(f"/eos/user/{username[0]}/{username}")
    
    # make condor file
    condor_template_file = open(f"{condor_dir}/submit.sub")
    local_condor = f"{condor_dir}/{jobname}.sub"
    condor_file = open(local_condor, "w")
    for line in condor_template_file:
        line = line.replace("DIRECTORY", str(condor_dir))
        line = line.replace("JOBNAME", jobname)
        condor_file.write(line)
    condor_file.close()
    condor_template_file.close()

    # make executable file
    sh_template_file = open(f"{condor_dir}/submit.sh")
    local_sh = f"{condor_dir}/{jobname}.sh"
    sh_file = open(local_sh, "w")
    for line in sh_template_file:
        line = line.replace("MAINDIRECTORY", str(main_dir))
        line = line.replace("COMMAND", cmd)
        line = line.replace("EOSDIR", str(eos_dir))
        sh_file.write(line)
    sh_file.close()
    sh_template_file.close()

    # submit jobs
    #subprocess.run(["condor_submit", local_condor])

    
def main(args):
    # load dataset config
    dataset_config = load_dataset_config(config_name=args.sample)
    nsplits = dataset_config.nsplit
    if nsplits == 1:
        jobname = f"{args.processor}_{args.channel}_{args.lepton_flavor}_{args.sample}"
        cmd = get_command(args)
        submit_condor(jobname, cmd)
    else:
        for nsplit in range(1, dataset_config.nsplit + 1):
            jobname = f"{args.processor}_{args.channel}_{args.lepton_flavor}_{args.sample}_{nsplit}"
            cmd = get_command(args, nsplit)
            submit_condor(jobname, cmd)
    
    
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
    args = parser.parse_args()
    main(args)