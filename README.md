**Overview**

**CLUES2Companion** is a comprehensive three-phase pipeline for inferring and dating selection coefficients using **Relate** (Speidel et al. 2019) or **SINGER** (Stern et al. 2024) as alternative ARG (ancestral recombination graph) inference engines.  

Depending on the chosen pathway, the pipeline automatically manages input preparation, genealogical inference, derived allele frequency estimation, CLUES2 execution, and visualization of results.  

---

- **Phase 1 (Relate pathway)**:  
  - Converts `.vcf` to `.haps` & `.sample`  
  - Prepares input (masking, polarization with ancestral FASTA, recombination maps)  
  - Runs **Relate** (producing `.anc`, `.mut`, `.coal`)  
  - Extracts SNP list, per-SNP derived allele encodings, and frequency table  

- **Phase 1 (SINGER pathway)**:  
  - Runs **SINGER** directly on phased `.vcf` files over a user-specified genomic interval  
  - Requires minimal input (mutation rate, optional Ne, recombination/mutation maps)  
  - Produces ARGs in `.trees` format  
  - Extracts SNP list and computes allele frequency table (ALT-based by default, polarization optional)  

---

- **Phase 2 (Relate pathway)**:  
  - Applies importance sampling of branch lengths with Relate (using `.coal`)  
  - Generates per-SNP `.newick` genealogies  
  - Runs **RelateToCLUES.py** to create per-SNP `*_times.txt`  
  - Runs **CLUES2** (`inference.py`) for selection inference (`*_inference.txt` + optional `*_CI.txt`)  
  - Merges SNP results into a `.tsv` summary file  
  - Generates integrative plots with allele frequencies, logLR, −log₁₀(p), and confidence intervals  

- **Phase 2 (SINGER pathway)**:  
  - Uses **SingerToCLUES.py** to extract per-SNP genealogical times from ARGs  
  - Runs **CLUES2** (`inference.py`) with user-defined Ne (`--N`) instead of Relate `.coal`  
  - Supports the same options as Relate (df, tCutoff, ancient DNA, noTraj, CI, timeBins, dominance h)  
  - Produces SNP-level `.inference.txt` and merged `.tsv` summary  
  - Generates plots identically to Relate pathway for consistency  

---

- **Phase 3** (common to both Relate and SINGER):  
  - Scans across time windows to date the onset of selection at a target SNP  
  - Produces onset estimates in generations and years  
  - Optionally performs bootstraps with user-defined settings to compute confidence intervals  
  - Saves results in JSON (`*_Dating.json`, `*_Bootstrap_onset.json`) and log files  

All outputs for Phase 1, Phase 2, and Phase 3 are automatically saved to organized folders created by CLUES2Companion, e.g.:  


---

**[Important]**: Users must keep the following folders exactly where they are:

- `Relate/` or `Relate-Linux/` (binaries + helper scripts)  
- `Singer/` or `Singer-Linux/` (binaries)  
- `CLUES2/` (inference.py, RelateToCLUES.py, SingerToCLUES.py, etc.)  
- `required_files/` (ancestral FASTA, recombination maps, masks, etc.)  
- auto-generated `phase1/`, `phase2/`, `phase3/` directories  

Moving or renaming **any** of these folders – or the main script itself – will interfere with the functioning of Phase 1, Phase 2, and Phase 3.

---

## Table of Contents

- [Installation](#Installation)  
- [Dependencies](#dependencies)  
- [Input files & folder layout](#inputs)  
- [Quick start](#quick-start)  
- [Phase 1 (details)](#Phase1)  
- [Phase 2 (details)](#Phase2)  
- [Phase 3 (details)](#Phase3)  
- [Tips & troubleshooting](#tips)  
- [License](#license)  
- [Contact](#contact)  

---

<a name="Installation"></a>
## Installation

There are multiple ways to install and use **CLUES2Companion**, depending on your preferences and computing environment.  

Before running any part of the pipeline, please read the file `README_before_use.txt`.  

---

### 1. GitHub (recommended)

Clone the repository from GitHub and move into the directory:

```bash
git clone https://github.com/<YourRepo>/CLUES2Companion.git
cd CLUES2Companion
```

### 2. Release package (precompiled)

Download the latest release (.tar.gz or .zip) from the Releases page.
Each release includes:
- CLUES2Companion shell and Python scripts
- Required folders (CLUES2/, required_files/)
- Example input data
- A Clues2Companion.yml Conda environment file

Unpack the archive and set permissions:
```
tar -xvzf CLUES2Companion_vX.Y.tar.gz
cd CLUES2Companion
chmod +x CLUES2Companion-Linux.sh
```

### 3. Dropbox mirror (backup)
As a convenience, we also maintain a Dropbox mirror of the package:

https://www.dropbox.com/scl/fo/m5y6aek0twd1jz9grg4p3/ALxMgIljUJRIZZNQXaGU-OE?rlkey=mbbh36ondftnqg0x07eao57eg&st=lryscqdm&dl=0


<a name="dependencies"></a>

## Dependencies

We provide a Conda YAML file (Clues2Companion.yml) that installs all Python dependencies for the pipeline:

```
conda env create -f Clues2Companion.yml
conda activate Clues2Companion
```
Note: Relate and SINGER are not Conda packages.
You can either:
Use the precompiled versions included in the GitHub release (Relate-Linux/, Singer-Linux/), or download them from their official repositories and rename the folders as Relate/ and Singer/ (or Relate-Linux/, Singer-Linux/) to match our default scripts.


Required folder layout

```
Ensure your working directory contains the following subfolders:
Relate/ or Relate-Linux/   contains Relate binaries and helper scripts  
Singer/ or Singer-Linux/   contains SINGER binary  
CLUES2/                    contains CLUES2 scripts (inference.py, RelateToCLUES.py, SingerToCLUES.py, …)  
required_files/            contains ancestral FASTA, recombination maps, masks, etc.  
phase1/, phase2/, phase3/  generated automatically by the pipeline
Do not move or rename these folders, otherwise Phase 1–3 will fail.
If you need to reorganize, update the paths inside the scripts accordingly.
```
If using Conda, these packages will be installed automatically from Clues2Companion.yml.
If installing manually, you can install the Python dependencies with:
```
pip3 install numpy pandas matplotlib adjustText biopython cyvcf2 numba tskit
```
Complete list of dependencies to install if manually installation was choosen:

Python ≥ 3.8 \
Relate ≥ v1.2.2 \
SINGER (latest master) \
CLUES2 (master branch) \
GNU parallel (optional, for bootstraps and batch runs) \

Python packages:
numpy \
pandas \
matplotlib \
adjustText \
biopython \
cyvcf2 \
numba \
tskit \
```

<a name="inputs"></a>

Input file & folder layout tree

```
CLUES2 Companion/
├── CLUES2Companion.sh         # main driver
├── Relate/                    # Relate bin/, scripts/
├── CLUES2/                    # CLUES v2 scripts
├── required_files/            # reference data
│   ├ ancestor/…
│   ├ mask/…
│   └ map/…
└── example/
    ├ Finnish_chr2.vcf.gz      # phased and indexed vcf
    └ Finnish.poplabels        # sampleID, population, group, SEX
```
Users must use fully phased vcf and indexed files (e.g., *.vcf.gz and *.tbi files). If your input data are GRCh38/hg38, please ensure chromosomes are encoded with prefix 'chr' (e.g. chr20). Alternatively, if your input data are GRCh37/hg19 please ensure chromosomes are encoded without prefix (e.g. 20)

Furthermore, the vcf file must contain only one population and one chromosome (e.g., example/Finnish_chr2.vcf.gz)

In addition, the *.poplabels file must contain four columns (namely, sampleID, population, group, SEX)

Example of *.poplabels file for diploid organisms:
```
sample population group sex
UNR1 FIN EUR NA
UNR2 FIN EUR NA
UNR3 FIN EUR NA
UNR4 FIN EUR NA
```
For more information please refer to: `https://myersgroup.github.io/relate/input_data.html#Prepare`

---

<a name="quick-start"></a>

## Quick start:
```
./CLUES2Companion.sh
```

Menu prompts: choose Phase 1, Phase 2 or Phase 3.


<a name="Phase1"></a>

## Phase 1 – Run Relate to create input files for CLUES2 

1 - Convert vcf to `*.haps`, `*.sample` \
2 - Apply PrepareInputFiles.sh in Relate (to mask, flip, and filter SNPs) \
3 - Run Relate mode All (with a user-specified Ne to generate the `*.anc` and `*.mut` files \
4 - Apply EstimatePopulationSize.sh to generate the `*.coal` file \
5 - Re-estimate branch lengths using the `*.coal` to generated updated `*.anc` and `*.mut` files \
6 - Use the cyvcf2 package to extract, SNPs and corresponding positions within a user-specified target region \
7 - Polarize derived alleles and compute derived allele frequency \

## Example of usage of Phase 1

```
./CLUES2Companion.sh
******  CLUES2 Companion – please cite CLUES2 and CLUES2Companion  ******
Choose phase to run
  1) Phase 1  : Apply Relate (*.mut, *.anc, *.coal files along with derived allele frequency) 
  2) Phase 2  : Apply Relate (BranchLengths) and CLUES2 (RelateToCLUES.py and inference.py) 
  3) Phase 3  : Date onset of selective sweeps of target SNP(s)
Enter option (1/2/3):
```
Here, users have to make their choice. That is to say, they will have to enter either 1, 2, or 3.

`then`

```
Enter option (1/2/3): 1   


          ╔═════════════════════════════════════════════════════════╗
          ║        🚀  PHASE 1: RELATE & SNP EXTRACTION  🚀        ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose the chromosome to analyze (e.g. 2, 17, X): 2
Prefix of phased vcf/bcf file (e.g. example/Finnish without _chrN.vcf.gz): example/Finnish
Path to population‑labels file (*.poplabels)(e.g. example/Finnish.poplabels): example/Finnish.poplabels
Start bp of target region: 135839626
End bp of target region: 135876443
Prefix of output name (e.g. Finnish): Finnish
```
### Explanation of Phase 1: required inputs and on-screen feedback

`1 - Select the input data:`

| Prompt                                                          | What to type                   |
| --------------------------------------------------------------- | ------------------------------ |
| `Chromosome to analyze (e.g. 2, 17, X):`                        | chromosome number (e.g., 2)    |
| `Prefix of phased vcf/bcf (Finnish without _chrN.vcf.gz):`      | e.g. `example/Finnish`         |


`2 - Define the target region"`

| Prompt                       | Description                                       |
| ---------------------------- | ------------------------------------------------- |
| `Start bp of target region:` | genomic start coordinate of the window to analyse |
| `End bp of target region:`   | genomic end coordinate of the window to analyse   |


All SNPs between these positions are extracted from the fully phased vcf.

`3 - Choose an output prefix:`

| Prompt                                  | Description                                                                                     |
| --------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `Prefix of output name:`                | e.g. `example/Finnish` **Important: You must use the exact same prefix in Phase 2 and Phase 3** |


### What you will see on the screen:

Every major step prints an INFO line that indicates where the resulting files are saved (e.g., output_C2Companion/phase1/)

While a step is running, a progress bar will move across the terminal. If something fails, you get a [WARN] or [ERROR] message; otherwise the step ends with the word `done!`.


<a name="Phase2"></a>

## Phase-2 – Selection coefficient inference

1 - Apply SampleBranchLengths.sh (Relate) to sample ancestral recombination graphs (ARGs) to generate Newick trees \
2 - Apply RelateToCLUES.py (CLUES2) to generate *_times.txt \
3 - Apply inference.py (CLUES2) to estimate selection coefficients and confidence intervals \
4 - Merge summary statistics, rsID, genomic coordinates, and derived allele frequency into *.tsv file \
5 - Generate tabular and graphical outputs \

---

### Example of usage for Phase 2

```
******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******
Choose phase to run
  1) Phase 1  : Apply Relate (*.mut, *.anc, *.coal files along with derived allele frequency) 
  2) Phase 2  : Apply Relate (BranchLengths) and CLUES2 (RelateToCLUES.py and inference.py) 
  3) Phase 3  : Date onset of selective sweeps of target SNP(s)
Enter option (1/2/3): 2

```

`then`

```
          ╔═════════════════════════════════════════════════════════╗
          ║        🧬  PHASE-2 – (requires Phase‑1 outputs)  🧬     ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose the chromosome to analyze (e.g. 2, 17, X): 2
Enter population prefix used in Phase‑1 (e.g., Finnish_MCM6): Finnish
Phase 1 auto-detect directory [ENTER = /~/CLUES2Companion2/output_C2Companion/phase1/Finnish_MCM6_chr2] or provide a different folder with phase1 outputs :

→ Using frequency file: ~/CLUES2Companion2/output_C2Companion/phase1/FIN_MCM6_chr2/Finnish_MCM6_Frequency_chr2_135839626_135876443.txt
→ Using SNPs file: ~/CLUES2Companion2/output_C2Companion/phase1/FIN_MCM6_chr2/Finnish_MCM6_SNPs_chr2_135839626_135876443.txt

tCutoff (e.g. 1000): 500
df (e.g. 600): 600
importance sampling of branch lengths: 200
AncientSamps file (optional, ENTER to skip): 
AncientHaps  file (optional, ENTER to skip): 
Disable allele trajectory? (y/N): y
Confidence interval (e.g. 0.95) [ENTER to skip]: 0.95
TimeBins (a list of epoch breakpoints; optional, ENTER to skip):
```

## Explanation of above command lines

### Script prompts:
Select the chromosome to analyze:
**Chromosome to analyze (e.g., 2, 17, X):**
If you run Phase 1 on several chromosomes, the pipeline automatically searches for the Phase 1 output that matches the chromosome that you enter above.

**Enter the population prefix:**
Provide the same prefix you used in Phase 1 (e.g. Finnish). The Phase 2 script now knows the exact folder structure (e.g., output_C2Companion/phase1/Finnish_chr2/).

**Confirm detected Phase 1 folders:**
The script shows the path to the location of the following file: *.anc, *.mut, *.coal, SNP list (`*.txt`), and frequency table (`*.txt`).
If the path is correct, press Enter to accept.
If users have manually moved Phase 1 folder (which is not recommended), type the path to the location of the new directory and press Enter.

### Mandatory CLUES2 parameters (please see CLUES2 and Relate manuals):

| Prompt                                    | Description                                          | Reference              |
| ----------------------------------------- | ---------------------------------------------------- | ---------------------- |
| `tCutoff (e.g. 1000):`                    | The maximum time (in generations)                    | CLUES2 `--tCutoff`    |
| `df (e.g. 600):`                          | number of allele frequency bins                      | CLUES2 `--df`         |
| `Importance sampling of branch lengths:`  | number of trees to sample from *Relate*              | Relate `--num_samples` |

None of the above prompts can be left blank.

### Optional files (please see CLUES2 and Relate manuals):

`AncientSamps — table of sample ages`
`AncientHaps — ancient haplotypes`
`Press Enter to skip either file`

Please see the CLUES2 manual for examples of formatting

### Optional parameters to include (please see CLUES2 and Relate manuals):

| prompt                                  | effect                                                              |
| --------------------------------------- | ------------------------------------------------------------------- |
| `Disable allele trajectory? (y/N):`     | answer **y** to *skip* posterior trajectories speeds up inference   |
| `Confidence intervals (e.g. 0.95):`     | compute CI for the selection coefficient (MLE); **Enter** = no CI   |
| `TimeBins (e.g., 200 300):`             | split 0-*tCutoff* into custom intervals; *Enter* = single epoch     |


The example TimeBins 200 300 for a tCutoff 500 yields three epochs: 0-200, 200-300, 300-500.

### Internal housekeeping:

NOTE: All output files will be saved to `~/Output_C2Companion/phase2/{prefix}_chr{N}/`. In this example output files are saved in `~/Output_C2Companion/phase2/FIN_MCM6_original_chr2/)`

### The resulting plot will contain:

SNPs analyzed in Phase 2 \
Confidence intervals \ 
Color intensity bar `[−log₁₀(p)]` \ 
Astericks indicating significance above each SNP (* (0.05), ** (0.01), *** (0.001)) 


### A *.tsv file containing SNP information, including related statistics calculated by CLUES2. 

**Example of `*.tsv file`**
```
rsID	POS	der_freq	logLR	-log10(p)	Epoch1_start	Epoch1_end	SelectionMLE1	95%_lower	95%_upper
rs55809728	135842606	0.0909	3.5508	2.11	0	536	0.09981	0.05336	0.14626
rs3754686	135845706	0.6515	5.3420	2.97	0	536	0.02312	0.00405	0.04218
rs4988243	135850133	0.1263	2.0500	1.37	0	536	0.09984	-0.01286	0.21254
rs4954490	135850661	0.6515	5.5070	3.04	0	536	0.02308	0.00369	0.04247
rs4988235	135851076	0.5909	17.2338	8.36	0	536	0.09986	0.08678	0.11294
rs4954493	135852405	0.6515	5.7506	3.16	0	536	0.02420	0.00069	0.04772
rs4988226	135853028	0.6515	5.3420	2.97	0	536	0.02312	0.00405	0.04218
rs309178	135854054	0.6515	5.5244	3.05	0	536	0.02288	0.00368	0.04208
rs309179	135856210	0.6515	5.4909	3.04	0	536	0.02248	0.00139	0.04356
```
[INFO] In addition to the final outputs, users will find three additional folders (`{prefix_inference_chr{N}/`, `prefix_times_chr{N}/`, `prefix_trees_chr{N}/` in which the output from inference.py, RelateToClues.py, and Relate are stored, respectively)

**FuAdditional reading**

CLUES2 manual: https://github.com/avaughn271/CLUES2#command-line-arguments \

Relate SampleBranchLengths: https://myersgroup.github.io/relate/modules.html#SampleBranchLengths \

---

<a name="Phase3"></a>


## Phase 3 – Dating selective sweeps

1 - Apply inference.py (CLUES2) to estimate selection coefficients and confidence intervals in indipendt multi-epochs bin \
2 - Generate Initial onset and writing result into a *.json file \
3 - (optional) Generate bootstraps to estimate confidence interval around first initial onset \
4 - Writing results into a *.json file with summary statistics \

The dating algorithm is based on inferring a selection coefficient in multiple epochs starting from the present and going back into the past. A selection coefficient is computed independently in each window. The assumption for dating is the observation of a selection coefficient that starts deviating from zero (s > 0) and maintains this trend for at least for four younger window. If the first criterion is not satisfied, the script searches for three consecutive windows in which s remains > 0. If there are not even 3 consecutive windows the script searches for 2 consecutive windows but will generate a warning message.


---

### Esample of usage for Phase 3

```
******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******
Choose phase to run
  1) Phase 1  : Apply Relate (*.mut, *.anc, *.coal files along with derived allele frequency) 
  2) Phase 2  : Apply Relate (BranchLengths) and CLUES2 (RelateToCLUES.py and inference.py) 
  3) Phase 3  : Date onset of selective sweeps of target SNP(s)
Enter option (1/2/3): 3

```

`then`

```
          ╔═════════════════════════════════════════════════════════╗
          ║        ⏱️  PHASE 3 – Dating a selective sweep  ⏱️         ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose chromosome to analyze (e.g. 2, 17, X): 2
Same prefix population name as used in Phase-2 (e.g. Finnish): Finnish 
rsID of SNP to date   (e.g. rs123): rs49988235
df score for CLUES2 (e.g. 600): 600
Initial epoch to scan for initial onset (e.g. 0 or 50): 0
Final epoch to scan for initial onset (e.g. 500 or 1000): 1000
Non overlapping windows size  [default is 50]: 50
```

## Explanation of above command lines

### Script prompts:
Select the chromosome to analyze:
**Chromosome to analyze (e.g., 2, 17, X):**
If you run Phase 1 and Phase 2 on several chromosomes, the pipeline automatically searches for the Phase 1 and Phase 2 outputs that matches the chromosome that you enter above.

**Enter the population prefix:**
Provide the same prefix you used in Phase 1 and Phase 2 (e.g., Finnish). Users must enter the same prefix as the one used in Phase 1 and Phase 2. The Phase 3 script automatically knows the location of input files (e.g., ~/output_C2Companion/phase1/Finnish_chr2/).

**Enter the rs identifier of SNP to date:**
Provide the `rs` identifier of SNP to date (e.g., rs4988235). Users can only analyze one variant per run. The rs identifier must be present in the phased vcf used in Phase 1.

**Enter the df score (e.g., 600):**
Provide the `df` (--df in CLUES2) score (which is the number of allele frequency bins). We recommend that users apply the same value as used in Phase 2 for comparability.

**Enter the starting time in generations ago to scan for the onset of selection:**
Provide the `start` time point from which to scan (e.g., 0 or 50)

**Enter the end time in generations ago to scan for the onset of selection:**
Provide the `end` time point from which to scan (e.g, 500 or 1000). Users should make sure that the end time in generstions ago ≥ putative sweep age. When unknown, 0–1000 generations ago is a safe default. Otherwise a final end time of 2000 generations ago will cover most possible scenarios

**Enter the non overlapping time bin size:**
Provide a value representing the initial breakpoint of two consecutive epochs (e.g, 25 or 50). If, as in the example, a starting epoch of 0 and a final (--tcutoff) of 500 and a windows size of 50 are provided, the script will automatically create: 0-50, 50-100, 100-150, 150-200 ..... 450-500 time bins in which selection coefficient will be estimated independentey. Smaller time bins mean longer computational time, but this approach will result in less noisy estimates.

```
e.g, :
Break-points : 50 (0-50)
50 (50-100)
100 (100-150)
150 (150-200)
200 (200-250)
250 (250-300)
300 (300-350)
350 (350-400)
400 (400-450)
450 (450-500)
500 (500-550)
550 (550-600)
600 (600-650)
650 (650-700)
700 (700-750)
750 (750-800)
800 (800-850)
850 (850-900)
900 (900-950)
950 (950-1000)
1000 (1000-1050)
tCutoff      : 1050
Window size  : 50
```

When the first dating is complete, a message with the summary results will be printed to the terminal along with a *.json file containing the first dating results:

```
For example, users will see the following on the terminal:

Initial onset ≈ 300 generations (≈ 8400 years)
epoch : 300 – 350
s(MLE) : 0.02187
```
```
e.g, on *.json:
{
  "rsID": "rs4988235",
  "population": "FIN_MCM63",
  "chromosome": "2",
  "onset_gen": 300,
  "onset_years": 8400,
  "epoch_start": 300,
  "epoch_end": 350,
  "epoch_start_years": 8400,
  "epoch_end_years": 9800,
  "s_MLE": 0.02187
}
```

`After the point estimate, the script will ask users if they would like to calculate the confidence intervals using a bootstrap approach:`

```
Proceed with bootstrap dating? [y/N]: y
```
`If the users type "y", they will be prompted to provide additioanl information (e.g., time bin size and number of replicates):`

```
▶  Bootstrap settings
Start point to scan in generations ago before initial onset (e.g. 200): 50
End point to scan in generations ago after initial onset (e.g. 500): 500
Non-overlapping time bin size  [default is 25]: 25
Number of bootstrap replicates  [default is 100]: 100
df score for CLUES2  [default is 450]: 600
Importance sampling of branch lengths: 200

Break-points : 25 (50-75)
75 (75-100)
100 (100-125)
125 (125-150)
150 (150-175)
175 (175-200)
200 (200-225)
225 (225-250)
250 (250-275)
275 (275-300)
300 (300-325)
325 (325-350)
350 (350-375)
375 (375-400)
400 (400-425)
425 (425-450)
450 (450-475)
475 (475-500)
500 (500-525)
tCutoff      : 525
Replicates   : 100
```
**[INFO]**: Runtime note: 100 bootstraps with default parameters typically require 16–32 CPU‑hours per SNP on a Desktop Intel 24-core computer with 32 Gb RAM.

To calculate confidence intervals, users will be prompted to provide a start and end time around the estimated onset of selection. For example, if the onset of selection was inferred to be 300 gen. ago (8,400 years ago), users could choose lower bound of 0 and an upper bound of 2,000 generations ago as starting point to calculate the confidence intervals around the onset of selection (in this case, 300 generations ago). Alternatively, to decrease the computational time, users could choose a smaller range, for example a lower bound of 0 and upper bound of 500 gen. ago to search for the upper and lower bounds of selection onset. Regardless of the choice, within the given range, our method will calculate s for each segmented non-overlapping window of size m (m=25 gen. is the default window size) and scan for consecutive windows with s > 0 and P < 0.05 to infer selection onset using the same trend-based procedure described above. The entire process of running Relate and CLUES2, including the calculation of s and P–values in segmented windows, is repeated n times (n the number of replicates defined by the user with a random seed, generating an empirical distribution of age of selection onset. CLUES2 Companion then will calculate the 2.5th and 97.5th percentiles of this distribution, which represent the lower and upper bounds of the confidence interval, respectively. 

```
For example, if users selected 100 replicates,

▶  Generating 100 bootstrap trees & CLUES runs
→ [bootstrap 1]  SampleBranchLengths (seed=973675772)
→   RelateToCLUES
→     CLUES inference • window [50 - 75]
→     CLUES inference • window [75 - 100]
→     CLUES inference • window [100 - 125]
→     CLUES inference • window [125 - 150]
→     CLUES inference • window [150 - 175]
→     CLUES inference • window [175 - 200]
→     CLUES inference • window [200 - 225]
→     CLUES inference • window [225 - 250]
→     CLUES inference • window [250 - 275]
→     CLUES inference • window [275 - 300]
→     CLUES inference • window [300 - 325]
→     CLUES inference • window [325 - 350]
→     CLUES inference • window [350 - 375]
→     CLUES inference • window [375 - 400]
→     CLUES inference • window [400 - 425]
→     CLUES inference • window [425 - 450]
→     CLUES inference • window [450 - 475]
→     CLUES inference • window [475 - 500]
```

Users must also decide the number of trees to sample (importance sampling using the SampleBrenchLenghts.sh script in Relate). For the `df` score and the `importance sampling`, we recommend that users apply the same parameters that were used in Phase 2 and Phase 3.

At the end of the bootstrap phase, the onset inferred by CLUES2 for each bootraps replication is printed to the terminal. These values ​​are used to calculate the summary statistics along with the *.json.

```
e.g, on *.json:
{
  "rsID": "rs4988235",
  "population": "Finnish_MCM6",
  "chromosome": "2",
  "CI95_low_gen": 155,
  "CI95_high_gen": 344,
  "CI95_low_year": 4340,
  "CI95_high_year": 9632,
  "method": "4 consecutive windows with s > 0"
}
```

<a name="license"></a>

License

MIT License

<a name="contact"></a>

Contact

Alessandro Lisi alisi@usc.edu (alisi1989)
Michael C. Campbell mc44680@usc.edu (mc44680)
