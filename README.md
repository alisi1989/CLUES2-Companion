**Overview**

**CLUES2 Companion** is a comprehensive set of pipelines for inferring selection coefficients using **Relate** (Speidel et al. 2019) or **SINGER** (Stern et al. 2024). In addition, our pipeline can infer the date of onset of selection on derived alleles using the **Relate** coalescent-based approach. 

Regardless of the software of choice (i.e., **Relate** or **SINGER**), CLUES2 Companion will manage input preparation, the genealogical inference, derived allele frequency estimation, CLUES2 execution, and visualization of results. To infer the age of selection onset, this method can only be performed using **Relate**.

---
**Phase 1**
- ### Relate-based Approach:  
  - Converts `.vcf` to `.haps` & `.sample`  
  - Prepares input (masking, polarization with ancestral FASTA, recombination maps)  
  - Runs **Relate** (producing `.anc`, `.mut`, `.coal`)  
  - Generates a SNP list, per-SNP derived allele polarization, and derived allele frequency table  

- ### SINGER-based Approach:  
  - Runs **SINGER** directly on phased `.vcf` files with a user-specified genomic interval  
  - Requires minimal input (mutation rate); optional parameters include Ne, recombination/mutation maps among others
  - Infers ARGs in `.trees` format file  
  - Generates a SNP list and computes allele frequency table (ALT-based by default)  

---

- **Phase 2**:
- ### Relate-based Approach:  
  - Applies importance sampling of branch lengths, which is given in the `.coal` file  
  - Generates per-SNP genealogies in `.newick` file  
  - Runs **RelateToCLUES.py** to create per-SNP genealogical times `*_times.txt` file
  - Runs **CLUES2** (`inference.py`) for selection inference (resulting files are `*_inference.txt` and the `*_CI.txt` which is optional)  
  - Merges SNP results into a `.tsv` summary file  
  - Generates integrative plots showing selection coefficients (s), associated confidence intervals, and −log₁₀(p)

- ### SINGER-based Approach:  
  - Uses **SingerToCLUES.py** to extract per-SNP genealogical times in `*_times.txt` file 
  - Runs **CLUES2** (`inference.py`) with user-defined Ne (`--N`) instead of the Relate `.coal` file to infer selection coefficient
  - Merges SNP results into a `.tsv` summary file  
  - Generates integrative plots showing selection coefficients (s), associated confidence intervals, and −log₁₀(p)  

---

- **Phase 3** (with **Relate**):  
  - Scans across time windows to date the onset of selection on a target SNP  
  - Produces onset estimates in generations and years  
  - Optionally performs bootstraps with user-defined settings to compute confidence intervals  
  - Saves results in JSON format (`*_Dating.json`, `*_Bootstrap_onset.json`) and associated log files  

---

**[Important]**: Users must keep the following directories exactly where they are:

- `Relate/` or `Relate-Linux/` (binaries + helper scripts)  
- `Singer/` or `Singer-Linux/` (binaries)  
- `CLUES2/` (e.g., inference.py, RelateToCLUES.py, SingerToCLUES.py)  
- `required_files/` (e.g., ancestral FASTA, recombination maps, masks)  
- auto-generated `phase1/`, `phase2/`, `phase3/` directories  

Moving or renaming **any** of these folders – or the main script itself – will interfere with the functioning of the Phase 1, Phase 2, and Phase 3 pipelines.

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
- [CLI version by command line](#cliversion)
- [Contact](#contact)  

---

<a name="Installation"></a>
## Installation

There are multiple ways to install and use **CLUES2Companion**, depending on your preferences and computing environment.  

Before running any part of this package, please read the `README_before_use.txt` file.  

---

### 1. GitHub (recommended)

Clone the repository from GitHub and move into the directory:

```bash
git clone https://github.com/<YourRepo>/CLUES2Companion.git
cd CLUES2Companion
```

### 2. Release package (precompiled)

Download the latest release (*.tar.gz or *.zip) from the "Releases" page.
Each release includes:
- CLUES2Companion shell (for non-expert users) and Python scripts (for expert users)
- Required folders (CLUES2/, required_files/)
- A Clues2Companion.yml Conda environment file (which includes required dependencies)

To use the CLUES2Companion shell script, users will need to apply the following steps: 
```
tar -xvzf CLUES2Companion_vX.Y.tar.gz
cd CLUES2Companion
chmod +x CLUES2Companion.sh
```

### 3. Dropbox mirror (backup)
Alternatively, we also maintain a Dropbox mirror of the package as a back-up:

https://www.dropbox.com/scl/fo/m5y6aek0twd1jz9grg4p3/ALxMgIljUJRIZZNQXaGU-OE?rlkey=mbbh36ondftnqg0x07eao57eg&st=lryscqdm&dl=0

Please note that the example files are present in this Dropbox repository

## Required folder layout

After downloading and decompressing the CLUES2 Companion package, users must ensure that the following directories and files are present:

```
Ensure your working directory contains the following subfolders:
Relate/ or Relate-Linux/   contains Relate binaries and helper scripts  
Singer/ or Singer-Linux/   contains SINGER binary  
CLUES2/                    contains CLUES2 scripts (inference.py, RelateToCLUES.py, SingerToCLUES.py, …)  
required_files/            contains ancestral FASTA, recombination maps, masks, etc.  

```
Do not move or rename these folders, otherwise Phases 1 through 3 will fail.

<a name="dependencies"></a>

## Python dependencies

Prior to installing the dependencies, users must ensure that their version of Python is 3.8 or higher. 

To manually install the dependencies, users can apply the following command:

```
pip3 install numpy pandas matplotlib adjustText biopython cyvcf2 numba tskit
for, was choosenComplete list of dependencies to install:

```
The complete list of dependencies is as follows:

numpy \
pandas \
matplotlib \
adjustText \
biopython \
cyvcf2 \
numba \
tskit \

Alternatively, we also provide a Conda YAML file (Clues2Companion.yml) that installs all Python dependencies as follows:

```
conda env create -f Clues2Companion.yml
conda activate Clues2Companion

```

<a name="inputs"></a>

---

## Input file 

Users must provide fully phased VCF and corresponding indexed files (e.g., `*.vcf.gz` and `*.tbi` index).  

- If your input VCF is build **GRCh38/hg38**, please ensure chromosome numbers have the prefix `chr` (e.g., `chr20`).  
- If your input VCF is build **GRCh37/hg19**, please ensure chromosome numbers do not have the prefix (e.g., `20`).  

The VCF file must consist of only one population and one chromosome (e.g., `example/Finnish_chr2.vcf.gz`).  

---

### Relate-based requirements

In addition to the phased VCF, the Relate-based approach requires auxiliary resources:  

- Ancestral FASTA file (e.g., `homo_sapiens_ancestor_chrN.fa`)  
- Pilot mask file (e.g., `PilotMask_chrN.fasta`)  
- Recombination map in Relate format (e.g., `genetic_map_chrN.txt`)  
- A population labels file (`*.poplabels`) that contains four columns: sampleID, population, group, sex  

The following is an example of a `*.poplabels` file for diploid organisms:  

```
sample population group sex
UNR1 FIN EUR NA
UNR2 FIN EUR NA
UNR3 FIN EUR NA
UNR4 FIN EUR NA
```

For more information, please refer to the **Relate** manual:  
`https://myersgroup.github.io/relate/input_data.html#Prepare`


### SINGER-based requirements

For the SINGER-based approach, the required files are:  

- Phased VCF and index files (e.g., `*_chr2.vcf.gz` + `*.tbi`)  
- The genomic interval of interest (start position and end position in bp) must be specified by the user  

Unlike Relate, **SINGER** does not require ancestral FASTA sequences, pilot mask files, or recombination maps.   

---

<a name="quick-start"></a>

## Quick start:
```
./CLUES2Companion.sh
```

Menu prompts: choose Phase 1, Phase 2, or Phase 3.


<a name="Phase1"></a>

## Phase 1 – Run Relate or SINGER to create input files for CLUES2

---

### Relate-based Approach

1. Convert VCF to `*.haps` and `*.sample`  
2. Apply `PrepareInputFiles.sh` in Relate (to mask, flip, and filter SNPs)  
3. Run Relate mode `All` (with a user-specified Ne) to generate the `*.anc` and `*.mut` files  
4. Apply `EstimatePopulationSize.sh` to generate the `*.coal` file  
5. Re-estimate branch lengths using the `*.coal` to generate updated `*.anc` and `*.mut` files  
6. Use the **cyvcf2** package to extract SNPs and corresponding positions within a user-specified target region  
7. Polarize derived alleles and compute derived allele frequencies  

---

### SINGER-based Approach

1. Provide a phased VCF (`*.vcf.gz` + index `*.tbi`) and specify the genomic interval of interest (start position and end position in bp)  
2. Run the `singer_master` binary on the chosen interval, providing the mutation rate (`-m`) and optionally the following parameters:  
   - effective population size (`--Ne`)  
   - recombination-to-mutation ratio (`--ratio`)  
   - recombination and/or mutation maps in SINGER format  
3. SINGER directly infers the Ancestral Recombination Graph (ARG) and outputs the results in `*.branches` file  
4. Convert the results in `*.branches` files to `*.trees` files using `convert_to_tskit.py` script (these files can be found in the phase1 output folder)  
5. Extract SNPs and positions from the VCF witin the user-specified range using **cyvcf2**
6. Calculate allele frequencies for the ALT allele, which are saved in the `<OUT_PREFIX>_Frequency_chr<CHR>_<START>_<END>.txt` file
7. All output files are saved in `output_CLUES2Companion-Singer/phase1/<OUT_PREFIX>_chr<CHR>/` and will serve as input files for Phase 2  


## Example of Phase 1 usage

```
./CLUES2Companion.sh
******  CLUES2Companion – please cite CLUES2, Relate, SINGER, and CLUES2Companion  ******
Choose the Phase to run
  1) Phase 1  : Apply Relate/SINGER (*.mut, *.anc, *.coal, and derived allele frequency file for Relate or *.trees files for SINGER)
  2) Phase 2  : Apply Relate or SINGER, and then run CLUES2 to infer selection coefficient (s) along with table and figure
  3) Phase 3  : Date onset of selective sweeps on target SNP(s) using Relate
Enter option (1/2/3):
```
Here, users have to make their choice. That is to say, they will have to enter either 1, 2, or 3.

```
Enter option (1/2/3): 1  
```

then

``` 
Which method do you want to use for genealogical inference?
  1) Relate
  2) SINGER
```
if 1 is selected, then users will proceed with Relate:

```
Enter option (1/2): 1
```
The following prompts will appear:
```   
          ╔═════════════════════════════════════════════════════════╗
          ║         🚀  PHASE 1: RELATE & SNP EXTRACTION  🚀        ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose the chromosome to analyze (e.g. 2, 17, X): 2
Prefix of phased vcf/bcf file (e.g. example/Finnish without _chrN.vcf.gz): example/Finnish
Path to population-labels file (*.poplabels)(e.g. example/Finnish.poplabels): example/Finnish.poplabels
Start bp of target region: 135839626
End bp of target region: 135876443
Prefix of output name (e.g. Finnish): Finnish_MCM6
```

Alternatively, if 2 is selected, then users will proceed with Singer:

The following prompts will appear:
```
Enter option (1/2): 2

          ╔═════════════════════════════════════════════════════════╗
          ║     🧬  PHASE 1: SINGER – ARGs INFERENCE FROM VCF  🧬     ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose the chromosome to analyze (e.g. 2, 17, X): 2
Prefix of phased VCF file (without _chrN.vcf): example/Finnish
Output prefix name (e.g. Bedouin_MCM6): Finnish_MCM6
Start bp of target region: 0
End bp of target region: 135876443
Mutation rate (-m, e.g. 1.25e-8): 1.25e-8
Effective population size (-Ne) [Optional]: 
Recombination/mutation ratio (-ratio) [Default: 1]: 
Recombination map file (-recomb_map)[Optional]: 
Mutation rate map file (-mut_map) [Optional]: 
Number of MCMC samples (-n) [Default: 100]: 
Thinning interval (-thin) [Default: 20]: 
Site flip probability (-polar) [Default: 0.5]: 
Random seed (-seed) [Default: 42]:
```
[IMPORTANT]: To work with SINGER, you must have a *.vcf file, choose an output name (prefix), the mutation rate, as well as the the start and end positions of genomic region of interest. \
[Note: The start position in SINGER must be 0. The end position will be end of the region of interest which is inputted by the user.]

All other parameters are optional and can be left at default or modified by the user according to their specific needs. For more information, please refer to the SINGER manual:
https://github.com/popgenmethods/SINGER

### Detailed explanation of Phase 1 (required inputs and on-screen feedback)

When running **Phase 1**, CLUES2 Companion will guide the user through a series of interactive prompts. Depending on which software **Relate** or **SINGER** is chosen, the required and optional inputs will differ slightly.

---

#### Relate-based Approach (required inputs)

| Prompt                                                          | Description                                                                 |
| --------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `Chromosome to analyze (e.g. 2, 17, X):`                        | Chromosome identifier                                                       |
| `Prefix of phased vcf/bcf file (without _chrN.vcf.gz):`         | Prefix of input VCF file (e.g. `example/Finnish`)                           |
| `Path to population-labels file (*.poplabels):`                 | File containing 4 columns: sampleID, population, group, sex (see example)    |
| `Start bp of target region:`                                    | Genomic start coordinate of the region to analyze                            |
| `End bp of target region:`                                      | Genomic end coordinate of the region to analyze                              |
| `Prefix of output name:`                                        | Output prefix (must be the same in Phases 1 through 3)                       |

*Note:* Relate requires additional resources in the `required_files/` folder, such as the ancestral FASTA, pilot mask, and genetic maps.

---

#### SINGER-based Approach (required inputs + optional parameters)

| Prompt                                          | Description                                                                 |
| ----------------------------------------------- | --------------------------------------------------------------------------- |
| `Chromosome to analyze (e.g. 2, 17, X):`        | Chromosome identifier                                                       |
| `Prefix of phased VCF file (without _chrN.vcf):`| Prefix of input VCF file (e.g. `example/Finnish`)                           |
| `Output prefix name:`                           | Output prefix (must be the same in Phases 1 through 3)                       |
| `Start bp of target region:`                    | Genomic start coordinate of the region to analyze [must be 0]                 |
| `End bp of target region:`                      | Genomic end coordinate of the region to analyze                              |
| `Mutation rate (-m):`                           | Mutation rate per bp per generation (required, e.g. `1.25e-8`)               |
| `Effective population size (-Ne):`              | Optional effective population size                                           |
| `Recombination/mutation ratio (-ratio):`        | Optional recombination-to-mutation ratio (default: 1)                        |
| `Recombination map file (-recomb_map):`         | Optional recombination map in SINGER format                                  |
| `Mutation rate map file (-mut_map):`            | Optional mutation rate map in SINGER format                                  |
| `Number of MCMC samples (-n):`                  | Number of samples for inference (default: 100)                               |
| `Thinning interval (-thin):`                    | Interval for thinning the MCMC chain (default: 20)                           |
| `Site flip probability (-polar):`               | Probability of flipping alleles (default: 0.5)                               |
| `Random seed (-seed):`                          | Random seed for reproducibility (default: 42)                                |


[IMPORTANT]: At minimum, users must provide a phased VCF, the region boundaries (start/end bp), an output prefix, and the mutation rate (`-m`). All other parameters are optional and default values will automatically be applied if left blank.

---

During the execution of CLUES2 Companion, the Phase 1 pipeline provides **on-screen feedback** with formatted banners (e.g., 🚀 Relate, 🧬 SINGER), status messages for each step, and automatic creation of the phase1 output folder (`output_CLUES2Companion-Relate/phase1/...` or `output_CLUES2Companion-Singer/phase1/...`).  
All necessary intermediate and final files (VCF-derived SNP list, frequency table, ARGs, and derived allele frequency) are conveniently saved for downstream Phase 2 analysis.



<a name="Phase2"></a>

## Phase 2 – Selection coefficient inference

### Relate-based Approach
1 - Apply `SampleBranchLengths.sh` (Relate) to sample ancestral recombination graphs (ARGs) and generate per-SNP Newick trees in`*.newick` files. 
2 - Apply `RelateToCLUES.py` (CLUES2) to convert each Newick tree into a `<rsID>_times.txt` file. 
3 - Apply `inference.py` (CLUES2) using the `<rsID>_times.txt`file and the derived allele frequency to estimate selection coefficients (s), confidence intervals, and −log10(p).
4 - Merge summary statistics (rsID, genomic coordinates, derived allele frequency, logLR, −log10(p), s estimates, CIs) into a single `*.tsv` file.  
5 - Generate tabular and graphical outputs (multi-SNP selection plots, error bars, −log10(p), etc.).  

### SINGER-based Approach
1 - Apply `SingerToCLUES.py` to extract per-SNP genealogical times directly from the ARGs (in `*.trees` files) generated in Phase 1  
2 - For each SNP in the user-defined region, produce `<rsID>_times.txt` files.  
3 - Apply `inference.py` (CLUES2) using the `<rsID>_times.txt` file and the allele frequency table from Phase 1; here, the effective population size parameter (`--N <Ne>`) is specified by the user instead of a Relate `.coal` file.
4 - Merge summary statistics (rsID, genomic coordinates, ALT/derived frequency, logLR, −log10(p), s estimates, confidence intervals) into a single `*.tsv` file.
5 - tabular and graphical outputs (e.g., multi-SNP selection plots, error bars, −log10(p)).  

### Key differences between the Relate- and SINGER- based approaches
- To infer branch lengths with **Relate**, the `*.coal` and `*.newick` files are required prior to running CLUES2. 
- To infer branch lengths with **SINGER**, the `*.newick` file is not required but instead SINGER uses the `.trees` file directly to generate genealogical times (via `SingerToCLUES.py` script). Furthermore, SINGER uses Ne `--N` to infer the branch length instead of the `.coal` file. 
- In the end, both Relate and SINGER generate `*.tsv` file that summarizes statistics for each SNP along with the graphical outputs. 

---

### Example of Phase 2 usage 

```
./CLUES2Companion.sh
******  CLUES2Companion – please cite CLUES2, Relate, SINGER, and CLUES2Companion  ******
Choose the Phase to run
  1) Phase 1  : Apply Relate/SINGER (*.mut, *.anc, *.coal, and derived allele frequency file for Relate or *.trees files for SINGER)
  2) Phase 2  : Apply Relate or SINGER, and then run CLUES2 to infer selection coefficient (s) along with table and figure
  3) Phase 3  : Date onset of selective sweeps on target SNP(s) using Relate
Enter option (1/2/3):
```
Here, users have to make their choice. That is to say, they will have to enter either 1, 2, or 3.

```
Enter option (1/2/3): 2 
```

then

``` 
Which method do you want to continue with?
  1) Relate
  2) SINGER
```
if 1 is selected, then users will proceed with Relate:

```
Enter option (1/2): 1
```
The following prompts will appear:
```   
Choose the chromosome to analyze (e.g. 2, 17, X): 2
Enter population prefix used in Phase-1 (e.g. Finnish): Finnish_MCM6
Phase-1 auto-detect directory [ENTER = /Users/CLUES2Companion/output_CLUES2Companion-Relate/phase1/Finnish_MCM6_chr2] or provide a different folder with phase1 outputs: 

→ Using frequency file: /Users/CLUES2Companion/output_CLUES2Companion-Relate/phase1/Finnish_MCM6_chr2/Finnish_MCM6_Frequency_chr2_135839626_135876443.txt
→ Using SNPs file: /Users/CLUES2Companion_Review/CLUES2Companion/output_CLUES2Companion-Relate/phase1/Finnish_MCM6_chr2/Finnish_MCM6_SNPs_chr2_135839626_135876443.txt

tCutoff (e.g. 1000): 500 
df (e.g. 600): 600
AncientSamps file (optional, press ENTER to skip): 
AncientHaps  file (optional, press ENTER to skip): 
Disable allele trajectory? (y/N): y
Confidence interval (e.g. 0.95) [press ENTER to skip]: 0.95
TimeBins (a list of epoch breakpoints; optional, press ENTER to skip): 
Dominance coefficient (default: 0.5, which is the additive model): 
Importance sampling of branch lengths: 200
```

if 2 is selected, then users will proceed with Singer:

```
Which method do you want to continue with?
  1) Relate
  2) SINGER
Enter option (1/2): 2
Choose the chromosome to analyze (e.g. 2, 17, X): 2
Enter population prefix used in Phase-1 (e.g. Finnish): Finnish_MCM6
Phase-1 auto-detect directory [ENTER = /Users/CLUES2Companion/output_CLUES2Companion-Singer/phase1/Finnish_chr22] or provide a different folder with phase1 outputs: 

→ Using frequency file: /Users/CLUES2Companion/output_CLUES2Companion-Singer/phase1/Finnish_chr2/Finnish_Frequency_chr2_35604500_35643818.txt
→ Using SNPs file: /Users/CLUES2Companion_Review/CLUES2Companion/output_CLUES2Companion-Singer/phase1/Finnish_chr2/Bedouin_SNPs_chr2_35604500_35643818.txt

tCutoff (e.g. 1000): 500
df (e.g. 600): 600
AncientSamps file (optional, press ENTER to skip): 
AncientHaps  file (optional, press ENTER to skip): 
Disable allele trajectory? (y/N): y
Confidence interval (e.g. 0.95) [press ENTER to skip]: 0.95
TimeBins (a list of epoch breakpoints; optional, press ENTER to skip): 
Dominance coefficient (default: 0.5, which is the additive model): 
Enter START bp of region: 135839626
Enter END bp of region: 135876443
Effective population size (Ne) for inference.py: 20000


```
## Explanation of above parameters

### Script prompts (shared between Relate and SINGER):
**Chromosome to analyze (e.g., 2, 17, X):**  
If you ran Phase 1 on several chromosomes, the pipeline automatically searches for the Phase 1 output that matches the chromosome in Phase 2.

**Population prefix (e.g. Finnish):**  
Provide the same prefix used in Phase 1. The Phase 2 script auto-detects the folder structure (e.g., `output_CLUES2Companion-Relate/phase1/Finnish_chr2/` or `output_CLUES2Companion-Singer/phase1/Finnish_chr2/`).

**Confirm detected Phase 1 folders:**  
The script shows the path to the SNP list (`*_SNPs.txt`) and frequency file (`*_Frequency.txt`). If the path is correct, press Enter to accept. Please note that if users manually move Phase 1 directory (not recommended), they must provide the new path to the directory.

---

### Mandatory CLUES2 parameters (shared):
| Prompt                                  | Description                                    | Reference           |
| --------------------------------------- | ---------------------------------------------- | ------------------- |
| `tCutoff (e.g. 1000):`                  | Maximum generations to consider                | CLUES2 `--tCutoff` |
| `df (e.g. 600):`                        | Number of allele frequency bins                | CLUES2 `--df`      |

These parameters cannot be left blank.

---

### Relate-specific prompts:
| Prompt                                    | Description                                              | Reference               |
| ----------------------------------------- | -------------------------------------------------------- | ----------------------- |
| `Importance sampling of branch lengths:`  | Number of trees to sample from Relate                    | Relate `--num_samples` |

Phase 2 will then:  
1. Run **SampleBranchLengths.sh** using the `.coal` file from Phase 1.  
2. Generate per-SNP Newick trees.  
3. Call **RelateToCLUES.py** to produce `*_times.txt`.  
4. Pass these outputs to **inference.py** for selection coefficient inference.

---

### SINGER-specific prompts:
| Prompt                                    | Description                                              | Reference               |
| ----------------------------------------- | -------------------------------------------------------- | ----------------------- |
| `Enter START bp of region:`               | Start coordinate of the target genomic window*           | SINGER input            |
| `Enter END bp of region:`                 | End coordinate of the target genomic window              | SINGER input            |
| `Effective population size (Ne):`         | Effective population size parameter for CLUES2 inference  | CLUES2 `--N`           |

*Users should note that the start position should not be 0, but instead it should be the actual start position (i.e. genomic coordinate) for the region of interest.

Phase 2 will then:  
1. Call **SingerToCLUES.py** use the information in the `.trees` output from Phase 1 to generate per-SNP `*_times.txt`.  
2. Run **inference.py**, passing `--N` (Ne) instead of a `.coal` file.  
3. Summarize results in a `.tsv` file and generate plots as for Relate.

---

### Optional parameters (for both Relate and SINGER):
| Prompt                                  | Effect                                                                 |
| --------------------------------------- | ---------------------------------------------------------------------- |
| `AncientSamps / AncientHaps`            | Include ancient samples if available                                   |
| `Disable allele trajectory? (y/N):`     | Answer **y** to skip posterior trajectories (speeds up inference)      |
| `Confidence intervals (e.g. 0.95):`     | Compute CI for the selection coefficient (s); Enter = no CI          |
| `TimeBins (e.g., 200 300):`             | Split 0–tCutoff into custom epochs (e.g. 0-200, 200-300, 300-500)      |
| `Dominance coefficient (default 0.5):`  | Define dominance (h); default is additive (h=0.5)                      |

---

### Internal housekeeping:

NOTE: All output files will be saved to `~/Output_C2Companion/phase2/{prefix}_chr{N}/`. In this example, output files are saved in `~/Output_C2Companion/phase2/Finnish_MCM6_chr2/)`

### The resulting plot will contain:

SNPs analyzed in Phase 2 \
Confidence intervals \ 
Color intensity bar `[−log₁₀(p)]` \ 
Astericks indicating significance above each SNP (specifically, * (0.05), ** (0.01), *** (0.001)) 


### The *.tsv file contains SNP information, including related statistics calculated by CLUES2. 

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

**Additional readings**

CLUES2 manual: https://github.com/avaughn271/CLUES2#command-line-arguments \

Relate SampleBranchLengths: https://myersgroup.github.io/relate/modules.html#SampleBranchLengths \

SINGER manual: https://github.com/popgenmethods/SINGER
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

---

<a name="cliversion"></a>

### Additional notes on available scripts for expert users:

In addition to the main Bash pipeline (`CLUES2Companion.sh`), we also provide fully equivalent **Python CLI scripts** for each phase. These scripts allow more advanced users to bypass interactive prompts and run the pipeline directly from the command line.

The Python versions (Phase 1, Phase 2, and Phase 3) are:

- `CLUES2Companion-Phase1.py`

To access help for CLI Phase 1, please type the following:
```
python CLUES2Companion-Phase1.py --help
```
Example of Phase 1 usage and parameters:
```
usage: CLUES2Companion-Phase1.py [-h] --engine {Relate,Singer} --chr CHR --vcf_prefix VCF_PREFIX --out_prefix OUT_PREFIX --start START --end END
                                           [--poplabels POPLABELS] [-m M] [-N N] [--mu MU] [--Ne NE] [--ratio RATIO] [--recomb RECOMB] [--mutmap MUTMAP] [-n N]
                                           [--thin THIN] [--polar POLAR] [--seed SEED]

CLUES2Companion Phase-1 CLI

options:
  -h, --help            show this help message and exit
  --engine {Relate,Singer}
                        Inference engine
  --chr CHR             Chromosome (e.g. 2, X) [Relate/Singer]
  --vcf_prefix VCF_PREFIX
                        Prefix of VCF (without _chrN.vcf/.vcf.gz)[Relate/Singer]
  --out_prefix OUT_PREFIX
                        Output prefix name [Relate/Singer]
  --start START         Start bp of region [Relate/Singer]
  --end END             End bp of region [Relate/Singer]
  --poplabels POPLABELS
                        Population labels file (*.poplabels) [Relate]
  -m M                  Mutation rate (Relate)
  -N N                  Ne (Relate)
  --mu MU               Mutation rate (Singer)
  --Ne NE               Effective population size (Singer)
  --ratio RATIO         Recomb/mutation ratio (Singer)
  --recomb RECOMB       Recombination map (Singer)
  --mutmap MUTMAP       Mutation map (Singer)
  -n N                  MCMC samples (Singer)
  --thin THIN           Thinning interval (Singer)
  --polar POLAR         Site flip probability (Singer)
  --seed SEED           Random seed (Singer)
```
To access help for CLI Phase 2, please type the following:
```
python CLUES2Companion-Phase2.py --help

```
Example of Phase 2 usage and parameters:
```
python CLUES2Companion-Phase2.py --help
usage: CLUES2Companion-Phase2.py [-h] --engine {Relate,Singer} --chr CHR --out_prefix OUT_PREFIX [--phase1_dir PHASE1_DIR] --tcutoff TCUTOFF --df DF
                                           [--h H] [--anc_samps ANC_SAMPS] [--anc_haps ANC_HAPS] [--noTraj] [--CI CI] [--timeBins TIMEBINS]
                                           [--num_samples NUM_SAMPLES] [--mu MU] [--region_start REGION_START] [--region_end REGION_END] [--Ne NE]

CLUES2Companion Phase-2 CLI (Relate/SINGER to CLUES2)

options:
  -h, --help            show this help message and exit
  --engine {Relate,Singer}
                        Continue with Relate or Singer branch
  --chr CHR             Chromosome (e.g. 2, X)
  --out_prefix OUT_PREFIX
                        Population/output prefix used in Phase-1
  --phase1_dir PHASE1_DIR
                        Override Phase-1 output directory (auto-detected by default)
  --tcutoff TCUTOFF     tCutoff (e.g. 1000)
  --df DF               df (e.g. 600)
  --h H                 Dominance coefficient (default: 0.5, additive model)
  --anc_samps ANC_SAMPS
                        AncientSamps file (optional)
  --anc_haps ANC_HAPS   AncientHaps file (optional)
  --noTraj              Disable allele trajectory
  --CI CI               Confidence interval (e.g. 0.95)
  --timeBins TIMEBINS   Space/comma separated epoch breakpoints (optional)
  --num_samples NUM_SAMPLES
                        Importance sampling of branch lengths (Relate)
  --mu MU               Mutation rate for SampleBranchLengths (Relate)
  --region_start REGION_START
                        START bp of region (Singer)
  --region_end REGION_END
                        END bp of region (Singer)
  --Ne NE               Effective population size for inference.py (Singer)
```
To access help for CLI Phase 3, please type the following:
```
python CLUES2Companion-Phase3.py --help

```
Example of Phase 3 usage and parameters:

```
python CLUES2Companion-Phase3.py --help
usage: CLUES2Companion-Phase3.py [-h] --chr CHR --out_prefix OUT_PREFIX --rsid RSID --df DF --g_start G_START --g_end G_END [--step STEP]
                                           [--phase1_dir PHASE1_DIR] [--phase2_dir PHASE2_DIR] [--bootstrap] [--boot_start BOOT_START] [--boot_end BOOT_END]
                                           [--boot_step BOOT_STEP] [--nboot NBOOT] [--boot_df BOOT_DF] [--num_samples NUM_SAMPLES] [--mu MU] [--h H]

CLUES2Companion Phase-3 CLI (Relate only) – Sliding-window dating + optional bootstrap

options:
  -h, --help            show this help message and exit
  --chr CHR             Chromosome (e.g. 2, X)
  --out_prefix OUT_PREFIX
                        Population/output prefix used in Phase-1/2
  --rsid RSID           Target SNP rsID to date (e.g. rs123)
  --df DF               df score for CLUES2 (e.g. 600)
  --g_start G_START     Initial epoch to scan (generations ago)
  --g_end G_END         Final epoch to scan (generations ago)
  --step STEP           Window size (non-overlapping), default 50
  --phase1_dir PHASE1_DIR
                        Override Phase-1 directory (auto-detected by default)
  --phase2_dir PHASE2_DIR
                        Override Phase-2 directory (auto-detected by default)
  --bootstrap           Enable bootstrap dating
  --boot_start BOOT_START
                        Bootstrap scan START (generations ago)
  --boot_end BOOT_END   Bootstrap scan END (generations ago)
  --boot_step BOOT_STEP
                        Bootstrap window size (default 25)
  --nboot NBOOT         Number of bootstrap replicates (default 100)
  --boot_df BOOT_DF     df score for CLUES2 in bootstrap (default 450)
  --num_samples NUM_SAMPLES
                        Importance sampling for SampleBranchLengths (required if --bootstrap)
  --mu MU               Mutation rate for SampleBranchLengths (default 1.25e-8)
  --h H                 Dominance coefficient (default: 0.5, additive model)


```

For users who already have already existing ARGs and wish to skip Phase 1 execution, they will use the `CLUES2Companion-Linux-ExtractFromPreProcessed.py` script:

To access help for the `CLUES2Companion-Linux-ExtractFromPreProcessed.py` script, please type the following:
```
python CLUES2Companion-Linux-ExtractFromPreProcessed.py --help

```
User should note that they will still need to extract the SNP and derived allele/ALT frequency information from the ARGs. 

Example of `CLUES2Companion-Linux-ExtractFromPreProcessed.py` usage and parameters:

```
python CLUES2Companion-ExtractFromPreProcessed.py --help
usage: CLUES2Companion-ExtractFromPreProcessed.py [-h] --engine {Relate,Singer} --chr CHR --vcf_prefix VCF_PREFIX --out_prefix OUT_PREFIX --start START
                                                        --end END [--haps HAPS] [--mut MUT]

options:
  -h, --help            show this help message and exit
  --engine {Relate,Singer}
  --chr CHR
  --vcf_prefix VCF_PREFIX
                        Prefix of VCF (without _chrN.vcf.gz)
  --out_prefix OUT_PREFIX
  --start START
  --end END
  --haps HAPS           HAPS file (Relate)
  --mut MUT             MUT file (Relate)
```

Again, this tool is designed for users who already have pre-existing ARGs and wish to skip full Phase 1 execution.
 
After using the `CLUES2Companion-Linux-ExtractFromPreProcessed.py` script, the extracted files and the pre-existing ARGs are ready to be used as inputs for **Phase 2** (without the need to re-run the entire Phase 1 pipeline).


---

<a name="license2"></a>

License

MIT License

<a name="contact"></a>

Contact

Alessandro Lisi alisi@usc.edu (alisi1989)
Michael C. Campbell mc44680@usc.edu (mc44680)
