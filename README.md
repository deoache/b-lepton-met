# W' + b

[![Codestyle](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<p align="left">
  <img width="300" src="https://i.imgur.com/OWhX13O.jpg" />
</p>

Python package for analyzing W' + b in the electron and muon channels. The analysis uses a columnar framework to process input tree-based [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) files using the [coffea](https://coffeateam.github.io/coffea/) and [scikit-hep](https://scikit-hep.org) Python libraries.


- [Data/MC filesets](Data/MC-filesets)
    * [Making the input filesets for Coffea-Casa](#Making-the-input-filesets-for-Coffea-Casa)
    * [Making the input filesets for Lxplus](#Making-the-input-filesets-for-Lxplus)
- [Submitting jobs](#Submitting-jobs)
    * [Submitting jobs at Coffea-Casa](#Submitting-jobs-at-Coffea-Casa)
    * [Submitting jobs at lxplus](#Submitting-jobs-at-lxplus)
- [Processors](#Processors)
    * [Trigger efficiency processor](#Trigger-efficiency-processor)
    * [tt processor](#tt-processor)
- [Corrections and scale factors](#Corrections-and-scale-factors)
- [Luminosity](#luminosity)


## Data/MC filesets

We use the recomended Run-2 UltraLegacy [datasets](https://github.com/deoache/wprime_plus_b/blob/main/wprime_plus_b/fileset/das_datasets.json). See https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis. 

#### Making the input filesets for Coffea-Casa

To build the input data/MC filesets to be used in Coffea-Casa use the [make_fileset.py](https://github.com/deoache/wprime_plus_b/blob/main/wprime_plus_b/fileset/make_fileset.py) script:
```
# connect to lxplus 
ssh <your_username>@lxplus.cern.ch

# then activate your proxy
voms-proxy-init --voms cms

# clone the repository 
git clone https://github.com/deoache/wprime_plus_b.git

# move to the fileset directory
cd wprime_plus_b/wprime_plus_b/fileset/

# run the 'make_fileset' script
python3 make_fileset.py
```

#### Making the input filesets for Lxplus

It has been observed that, in lxplus, opening files through a concrete xrootd endpoint rather than a redirector is far more robust. Use the [make_fileset_lxplus.py](https://github.com/deoache/wprime_plus_b/blob/main/wprime_plus_b/fileset/make_fileset_lxplus.py) script to build the input filesets with xrootd endpoints:
```
# connect to lxplus 
ssh <your_username>@lxplus.cern.ch

# then activate your proxy
voms-proxy-init --voms cms

# clone the repository 
git clone https://github.com/deoache/wprime_plus_b.git

# move to the fileset directory
cd wprime_plus_b/wprime_plus_b/fileset/

# get the singularity shell 
singularity shell -B /afs -B /eos -B /cvmfs /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask:latest-py3.10

# run the 'make_fileset_lxplus' script
python make_fileset_lxplus.py
```
We use the [dataset discovery tools](https://coffeateam.github.io/coffea/notebooks/dataset_discovery.html) from Coffea 2024, that's why we need to use a singularity shell in which we can use these tools.


The json files containing the datasets will be saved in the `wprime_plus_b/fileset/coffea-casa` or `wprime_plus_b/fileset/lxplus` directories, depending on which script is executed.


## Submitting jobs

[Coffea-Casa](https://coffea-casa.readthedocs.io/en/latest/cc_user.html) is easier to use and more convenient for beginners, however still somewhat experimental, so for large inputs and/or processors which may require heavier cpu/memory using HTCondor at lxplus is recommended. 

Before submitting jobs, make sure you have the facility (`coffea-casa` or `lxplus`) datasets by typing:

```bash
python3 process_filesets.py --<facility>
```

### Submitting jobs at Coffea-Casa

The `submit.py` file executes a desired processor with user-selected options. To see a list of arguments needed to run this script please enter the following in the terminal:

```bash
python3 submit.py --help
```
The output should look something like this:

```
usage: submit.py [-h] [--processor PROCESSOR] [--channel CHANNEL] [--lepton_flavor LEPTON_FLAVOR] [--sample SAMPLE] [--year YEAR] [--yearmod YEARMOD]
                 [--executor EXECUTOR] [--workers WORKERS] [--nfiles NFILES] [--nsample NSAMPLE] [--chunksize CHUNKSIZE] [--output_type OUTPUT_TYPE]
                 [--syst SYST] [--facility FACILITY] [--tag TAG]

optional arguments:
  -h, --help            show this help message and exit
  --processor PROCESSOR
                        processor to be used {ttbar, ztoll, qcd, trigger_eff, btag_eff} (default ttbar)
  --channel CHANNEL     channel to be processed {'2b1l', '1b1e1mu', '1b1l'}
  --lepton_flavor LEPTON_FLAVOR
                        lepton flavor to be processed {'mu', 'ele'}
  --sample SAMPLE       sample key to be processed
  --year YEAR           year of the data {2016, 2017, 2018} (default 2017)
  --yearmod YEARMOD     year modifier {'', 'APV'} (default '')
  --executor EXECUTOR   executor to be used {iterative, futures, dask} (default iterative)
  --workers WORKERS     number of workers to use with futures executor (default 4)
  --nfiles NFILES       number of .root files to be processed by sample. To run all files use -1 (default 1)
  --nsample NSAMPLE     partitions to run (--nsample 1,2,3 will only run partitions 1,2 and 3)
  --chunksize CHUNKSIZE
                        number of chunks to process
  --output_type OUTPUT_TYPE
                        type of output {hist, array}
  --syst SYST           systematic to apply {'nominal', 'jet', 'met', 'full'}
  --facility FACILITY   facility to launch jobs {coffea-casa, lxplus}
  --tag TAG             tag to reference output files directory
```

* The processor to be run is selected using the `--processor` flag. 
* According to the processor, you can choose channel and lepton flavor by means of the `--channel` and `--lepton_flavor` flags
* You can select a particular sample with `--sample <sample_name>` (see samples names [here]((https://github.com/deoache/wprime_plus_b/blob/main/wprime_plus_b/fileset/das_datasets.json)))
* The year can be selected using the `--year` flag, and the `--yearmod` flag is used to specify whether the dataset uses APV or not.
* You can select the executor to run the processor using the `--executor` flag. Three executors are available: `iterative`, `futures`, and `dask`. The `iterative` executor uses a single worker, while the `futures` executor uses the number of workers specified by the `--workers` flag. The `dask` executor uses Dask functionalities to scale up the analysis (only available at coffea-casa).
* To lighten the workload of jobs, the fileset is divided into sub-filesets. The number of partitions per dataset can be defined [here](https://github.com/deoache/wprime_plus_b/blob/main/wprime_plus_b/configs/dataset/datasets_configs.yaml). Set `--nfiles -1` to use all `.root` files.
* You can set `--nsample <n>` to run only the `n` partition of the selected dataset.
* The output type of the processor (histograms or arrays) is defined with the `output_type` flag.
* If you choose histograms as output, you can add some systematics to the output. With `--syst nominal`, variations of the scale factors will be added. With `jet` or `met`, JEC/JER or MET variations will be added, respectively. Use `full` to add all variations. 
* The selected processor is executed at some facility, defined by the `--facility` flag.  

Let's assume we want to execute the `ttbar` processor, in the `2b1l` electron control region, using the `TTTo2L2Nu` sample from 2017. To test locally first, can do e.g.:

```bash
python3 submit.py --processor ttbar --channel 2b1l --lepton_flavor ele --executor iterative --sample TTTo2L2Nu --nfiles 1
```
Then, if everything is ok, you can run the full dataset with:

```bash
python submit.py --processor ttbar --channel 2b1l --lepton_flavor ele --executor futures --sample TTTo2L2Nu --nfiles -1
```
The results will be stored in the `wprime_plus_b/outs` folder

### Submitting condor jobs at lxplus 

To submit jobs at lxplus using HTCondor, you need to have a valid grid proxy in the CMS VO. (see [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideLcgAccess) for details on how to register in the CMS VO). The needed grid proxy is obtained via the usual command
```bash
voms-proxy-init --voms cms
```
To execute a processor using some sample of a particular year type:
```bash
python3 submit_lxplus.py --processor ttbar --channel 2b1l --lepton_flavor ele --sample TTTo2L2Nu --year 2017 --nfiles -1
```
The script will create the condor and executable files (using the `submit.sub` and `submit.sh` templates) needed to submit jobs, as well as the folders containing the logs and outputs within the `/condor` folder (click [here](https://batchdocs.web.cern.ch/local/quick.html) for more info). After submitting the jobs, you can watch their status typing
```bash
watch condor_q
```
The output will be save to your EOS area. 

#### Notes: 
* Currently, the processors are only functional for the year 2017. 

## Processors

### [Trigger efficiency](processors/trigger_efficiency_processor.py) 

Processor use to compute trigger efficiencies. 

The processor applies the following pre-selection cuts




| $$\textbf{Object}$$    | $$\textbf{Variable}$$          | $$\textbf{Cut}$$                                                    | 
| ---------------------  | ------------------------------ | ------------------------------------------------------------------- |
| $$\textbf{Electrons}$$ |                                |                                                                     |
|                        | $p_T$                        | $\geq 30$ GeV                                           |
|                        | $\eta$                       | $\| \eta \| < 1.44$ and $1.57 < \| \eta \| < 2.5$       |
|                        | pfRelIso04_all                 | $\lt 0.25$                                                        |
|                        | mvaFall17V2Iso_WP80 (ele) mvaFall17V2Iso_WP90 (mu) | $\text{True}$|
| $$\textbf{Muons}$$     |                                |                                                                     |
|                        | $p_T$                        | $\geq 30$ GeV                                          |
|                        | $\eta$                       | $\lt 2.4$                                                 |
|                        | pfRelIso04_all       | $\lt 0.25$                                                        |
|                        | mediumId (ele) tightId (mu) | $\text{True}$                   |
| $$\textbf{Jets}$$      |                                |                                                                     |
|                        | $p_T$                        |  $\geq 20$ GeV                                       |
|                        | $\eta$                       | $\lt 2.4$                                                 |
|                        | JetId                        | $6$                                                               |
|                        | puId                          | $7$                                                               |
|                        | btagDeepFlavB                | $\gt$ Medium WP                                          |



The trigger efficiency is computed as:


$$\epsilon = \frac{\text{selection cuts and reference trigger and main trigger}}{\text{selection cuts and reference trigger}}$$


so we define two regions for each channel: ```numerator``` and ```denominator```. We use the following triggers:


$\text{Analysis triggers}$


| Channel        | 2016           |    | 2017           |   | 2018           |
|----------------|----------------|--- |----------------|---|----------------|
| Muon           | IsoMu24        |    | IsoMu27        |   | IsoMu24        |
| Electron       | Ele27\_WPTight\_Gsf |   | Ele35\_WPTight\_Gsf |   | Ele32\_WPTight\_Gsf |


The reference and main triggers, alongside the selection criteria applied to establish each region, are presented in the following tables:


#### Electron channel

| Trigger        | 2016           |   | 2017           |   | 2018           |
|----------------|----------------|---|----------------|---|----------------|
| Reference trigger   | IsoMu24        |   | IsoMu27        |   | IsoMu24        |
| Main trigger         | Ele27\_WPTight\_Gsf |   | Ele35\_WPTight\_Gsf |   | Ele32\_WPTight\_Gsf |


| Selection cuts      | 
| ---------------------------------|
| Luminosity calibration            |
| MET filters                       |
| $N(bjet) \geq 1$                   |
| $N(\mu) = 1$                      |
| $N(e) = 1$                       |

#### Muon channel

| Trigger        | 2016           |   | 2017           |   | 2018           |
|----------------|----------------|---|----------------|---|----------------|
| Reference trigger         | Ele27\_WPTight\_Gsf |   | Ele35\_WPTight\_Gsf |   | Ele32\_WPTight\_Gsf |
| Main trigger   | IsoMu24        |   | IsoMu27        |   | IsoMu24        |



| Selection cuts      | 
| ------------------------------------------------------------------- |
| Luminosity calibration                                    |
| MET filters                        |
| $\Delta R (\mu, bjet) \gt 0.4$                           |
| $N(bjet) \geq 1$                           |
| $N(\mu) = 1$                        |
| $N(e) = 1$                       |


### [tt processor](processors/ttbar_analysis.py) 

Processor use to estimate backgrounds in two $t\bar{t}$ control regions (`2b1l` and `1b1e1mu`) and signal region (`1b1l`), in both $e$ and $\mu$ lepton channels.

**`2b1l` region**: The processor applies the following pre-selection cuts for the electron (ele) and muon (mu) channels:

| $$\textbf{Object}$$    | $$\textbf{Variable}$$          | $$\textbf{Cut}$$                                                    | 
| ---------------------  | ------------------------------ | ------------------------------------------------------------------- |
| $$\textbf{Electrons}$$ |                                |                                       |
|                        | $p_T$                        | $\geq 55$ GeV (ele)  $\geq 30$ GeV (mu)                                        |
|                        | $\eta$                       | $\| \eta \| < 1.44$ and $1.57 < \| \eta \| < 2.5$       |
|                        | pfRelIso04_all                 | $\lt 0.25$                                                        |
|                        | mvaFall17V2Iso_WP80 (ele) mvaFall17V2Iso_WP90 (mu) | $\text{True}$|
| $$\textbf{Muons}$$     |                                |                                                                     |
|                        | $p_T$                        | $\geq 35$ GeV                                         |
|                        | $\eta$                       | $\lt 2.4$                                                 |
|                        | pfRelIso04_all               | $\lt 0.25$                                                        |
|                        | tightId                      | $\text{True}$                   |
| $$\textbf{Taus}$$      |                                |                                                                     |
|                        | $p_T$                        | $\geq 20$ GeV                                               |
|                        | $\eta$                       | $\lt 2.3$                                                 |
|                        | $dz$                         | $\lt 0.2$                                                        | 
|                        | idDeepTau2017v2p1VSjet       | $\gt 8$                                                           |
|                        | idDeepTau2017v2p1VSe         | $\gt 8$                                                           |
|                        | idDeepTau2017v2p1VSmu        | $\gt 1$                                                           |
| $$\textbf{Jets}$$      |                                |                                                                     |
|                        | $p_T$                        |  $\geq 20$ GeV                                            |
|                        | $\eta$                       | $\lt 2.4$                                                 |
|                        | JetId                        | $6$                                                               |
|                        | puId                          | $7$                                                               |
|                        | btagDeepFlavB                | $\gt$ Medium WP                                          |




and additional selection cuts for each channel:

#### Electron channel

| Selection cuts      | 
| ---------------------------------|
| Electron Trigger      |
| Luminosity calibration                  |
| MET filters           |
| $p_T^{miss}\gt 50$ GeV |
| $N(bjet) = 2$                  |
| $N(\tau) = 0$                  |
| $N(\mu) = 0$                   |
| $N(e) = 1$                    |
| $\Delta R (e, bjet_0) \gt 0.4$ |

expected to be run with the `SingleElectron` dataset.

#### Muon channel


| Selection cuts      | 
| ---------------------------------|
| Muon Trigger          |
| Luminosity calibration                  |
| MET filters           |
| $p_T^{miss}\gt 50$ GeV |
| $N(bjet) = 2$                  |
| $N(\tau) = 0$                  |
| $N(e) = 0$                     |
| $N(\mu) = 1$                   |
| $\Delta R (\mu, bjet_0) \gt 0.4$ |

expected to be run with the `SingleMuon` dataset.
 

**`1b1e1mu` region:** Processor use to estimate backgrounds in a $t\bar{t}$ control region. 

The processor applies the following pre-selection cuts for the electron (ele) and muon (mu) channels:

| $$\textbf{Object}$$    | $$\textbf{Variable}$$          | $$\textbf{Cut}$$                                                    | 
| ---------------------  | ------------------------------ | ------------------------------------------------------------------- |
| $$\textbf{Electrons}$$ |                                |                                       |
|                        | $p_T$                        | $\geq 55$ GeV (mu)  $\geq 30$ GeV (ele)                                        |
|                        | $\eta$                       | $\| \eta \| < 1.44$ and $1.57 < \| \eta \| < 2.5$       |
|                        | pfRelIso04_all                 | $\lt 0.25$                                                        |
|                        | mvaFall17V2Iso_WP80 (ele) mvaFall17V2Iso_WP90 (mu) | $\text{True}$|
| $$\textbf{Muons}$$     |                                |                                                                     |
|                        | $p_T$                        | $\geq 35$ GeV                                         |
|                        | $\eta$                       | $\lt 2.4$                                                 |
|                        | pfRelIso04_all               | $\lt 0.25$                                                        |
|                        | tightId                      | $\text{True}$                   |
| $$\textbf{Taus}$$      |                                |                                                                     |
|                        | $p_T$                        | $\geq 20$ GeV                                               |
|                        | $\eta$                       | $\lt 2.3$                                                 |
|                        | $dz$                         | $\lt 0.2$                                                        | 
|                        | idDeepTau2017v2p1VSjet       | $\gt 8$                                                           |
|                        | idDeepTau2017v2p1VSe         | $\gt 8$                                                           |
|                        | idDeepTau2017v2p1VSmu        | $\gt 1$                                                           |
| $$\textbf{Jets}$$      |                                |                                                                     |
|                        | $p_T$                        |  $\geq 20$ GeV                                            |
|                        | $\eta$                       | $\lt 2.4$                                                 |
|                        | JetId                        | $6$                                                               |
|                        | puId                          | $7$                                                               |
|                        | btagDeepFlavB                | $\gt$ Medium WP                                          |



and additional selection cuts for each channel:

#### Electron channel

| Selection cuts      | 
| ---------------------------------|
| Muon Trigger      |
| Luminosity calibration                  |
| MET filters           |
| $p_T^{miss}\gt 50$ GeV |
| $N(bjet) = 1$                  |
| $N(\tau) = 0$                  |
| $N(\mu) = 1$                   |
| $N(e) = 1$                    |
| $\Delta R (e, bjet_0) \gt 0.4$ |
| $\Delta R (\mu, bjet_0) \gt 0.4$ |

expected to be run with the `SingleMuon` dataset.

#### Muon channel


| Selection cuts      | 
| ---------------------------------|
| Electron Trigger          |
| Luminosity calibration                  |
| MET filters           |
| $p_T^{miss}\gt 50$ GeV |
| $N(bjet) = 1$                  |
| $N(\tau) = 0$                  |
| $N(e) = 1$                     |
| $N(\mu) = 1$                   |
| $\Delta R (\mu, bjet_0) \gt 0.4$ |
| $\Delta R (e, bjet_0) \gt 0.4$ |

expected to be run with the `SingleElectron` dataset.

## Corrections and scale factors

We implemented particle-level corrections and event-level scale factors

### Particle-level corrections 

**JEC/JER corrections**: The basic idea behind the JEC corrections at CMS is the following: *"The detector response to particles is not linear and therefore it is not straightforward to translate the measured jet energy to the true particle or parton energy. The jet corrections are a set of tools that allows the proper mapping of the measured jet energy deposition to the particle-level jet energy"* (see https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC).

We follow the recomendations by the Jet Energy Resolution and Corrections (JERC) group (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC#Recommended_for_MC). In order to apply these corrections to the MC (in data, the corrections are already applied) we use the `jetmet_tools` from Coffea (https://coffeateam.github.io/coffea/modules/coffea.jetmet_tools.html). With these tools, we construct the [Jet and MET factories](wprime_plus_b/data/scripts/build_jec.py) which contain the JEC/JER corrections that are eventually loaded in the function [`jet_corrections`](wprime_plus_b/corrections/jec.py), which is the function we use in the processors to apply the corrections to the jet and MET objects.

**Note**: Since we modify the kinematic properties of jets, we must recalculate the MET. That's the work of the MET factory: it takes the corrected jets as an argument, and use them to recalculate the MET.

**Note:** These corrections must be applied before performing any kind of selection.

**MET phi modulation:** The distribution of true MET is independent of $\phi$ because of the rotational symmetry of the collisions around the beam axis. However, we observe that the reconstructed MET does depend on $\phi$. The MET $\phi$ distribution has roughly a sinusoidal curve with the period of $2\pi$. The possible causes of the modulation include anisotropic detector responses, inactive calorimeter cells, the detector misalignment, the displacement of the beam spot. The amplitude of the modulation increases roughly linearly with the number of the pile-up interactions.

We implement this correction [here](wprime_plus_b/corrections/met.py). This correction reduces the MET $\phi$ modulation. It is also a mitigation for the pile-up effects. 

(taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#7_7_6_MET_Corrections)

### Event-level scale factors (SF)

We use the common json format for scale factors (SF), hence the requirement to install [correctionlib](https://github.com/cms-nanoAOD/correctionlib). The SF themselves can be found in the central [POG repository](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration), synced once a day with CVMFS: `/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration`. A summary of their content can be found [here](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/). The SF implemented are:

* [Pileup SF](wprime_plus_b/corrections/pileup.py)
* [Electron ID, Reconstruction and Trigger SF](wprime_plus_b/corrections/lepton.py) (see the `ElectronCorrector` class)
* [Muon ID, Iso and TriggerIso ](wprime_plus_b/corrections/lepton.py) (see the `MuonCorrector` class)
* [PileupJetId SF](wprime_plus_b/corrections/pujetid.py)
* L1PreFiring SF: These are read from the NanoAOD events as `events.L1PreFiringWeight.Nom/Up/Dn`.

*We derive our own set of trigger scale factors. 

* B-tagging: b-tagging weights are computed as (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods):

  $$w = \prod_{i=\text{tagged}} \frac{SF_{i} \cdot \varepsilon_i}{\varepsilon_i} \prod_{j=\text{not tagged}} \frac{1 - SF_{j} \cdot \varepsilon_j}{1-\varepsilon_j} $$
  
  where $\varepsilon_i$ is the MC b-tagging efficiency and $\text{SF}$ are the b-tagging scale factors. $\text{SF}_i$ and $\varepsilon_i$ are functions of the jet flavor, jet $p_T$, and jet $\eta$. It's important to notice that the two products are 1. over jets tagged at the respective working point, and 2. over jets not tagged at the respective working point. **This is not to be confused with the flavor of the jets**.
  
  We can see, then, that the calculation of these weights require the knowledge of the MC b-tagging efficiencies, which depend on the event kinematics. It's important to emphasize that **the BTV POG only provides the scale factors and it is the analyst responsibility to compute the MC b-tagging efficiencies for each jet flavor in their signal and background MC samples before applying the scale factors**. The calculation of the MC b-tagging efficiencies is describe [here](https://github.com/deoache/wprime_plus_b/blob/refactor/corrections/binder/btag_eff.ipynb).

  The computation of the b-tagging weights can be found [here](wprime_plus_b/corrections/btag.py)



## Luminosity

See luminosity recomendations for Run2 at https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2. To obtain the integrated luminosity type (on lxplus):

```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH
pip uninstall brilws
pip install --install-option="--prefix=$HOME/.local" brilws
```

* SingleMuon: type

```
brilcalc lumi -b "STABLE BEAMS" --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --hltpath HLT_IsoMu27_v*
```

output:
```
#Summary: 
+-----------------+-------+------+--------+-------------------+------------------+
| hltpath         | nfill | nrun | ncms   | totdelivered(/fb) | totrecorded(/fb) |
+-----------------+-------+------+--------+-------------------+------------------+
| HLT_IsoMu27_v10 | 13    | 36   | 8349   | 2.007255669       | 1.870333304      |
| HLT_IsoMu27_v11 | 9     | 21   | 5908   | 1.383159994       | 1.254273727      |
| HLT_IsoMu27_v12 | 47    | 122  | 46079  | 8.954672794       | 8.298296788      |
| HLT_IsoMu27_v13 | 91    | 218  | 124447 | 27.543983745      | 26.259684708     |
| HLT_IsoMu27_v14 | 2     | 13   | 4469   | 0.901025085       | 0.862255849      |
| HLT_IsoMu27_v8  | 2     | 3    | 1775   | 0.246872270       | 0.238466292      |
| HLT_IsoMu27_v9  | 11    | 44   | 14260  | 2.803797063       | 2.694566730      |
+-----------------+-------+------+--------+-------------------+------------------+
#Sum delivered : 43.840766620
#Sum recorded : 41.477877399
```

* SingleElectron: type

```
brilcalc lumi -b "STABLE BEAMS" --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --hltpath HLT_Ele35_WPTight_Gsf_v*
```
output:

```
#Summary: 
+--------------------------+-------+------+--------+-------------------+------------------+
| hltpath                  | nfill | nrun | ncms   | totdelivered(/fb) | totrecorded(/fb) |
+--------------------------+-------+------+--------+-------------------+------------------+
| HLT_Ele35_WPTight_Gsf_v1 | 2     | 3    | 1775   | 0.246872270       | 0.238466292      |
| HLT_Ele35_WPTight_Gsf_v2 | 11    | 44   | 14260  | 2.803797063       | 2.694566730      |
| HLT_Ele35_WPTight_Gsf_v3 | 13    | 36   | 8349   | 2.007255669       | 1.870333304      |
| HLT_Ele35_WPTight_Gsf_v4 | 9     | 21   | 5908   | 1.383159994       | 1.254273727      |
| HLT_Ele35_WPTight_Gsf_v5 | 20    | 66   | 22775  | 5.399580877       | 4.879405647      |
| HLT_Ele35_WPTight_Gsf_v6 | 27    | 56   | 23304  | 3.555091917       | 3.418891141      |
| HLT_Ele35_WPTight_Gsf_v7 | 93    | 231  | 128916 | 28.445008830      | 27.121940558     |
+--------------------------+-------+------+--------+-------------------+------------------+
#Sum delivered : 43.840766620
#Sum recorded : 41.477877399
```