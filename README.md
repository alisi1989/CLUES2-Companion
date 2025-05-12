# CLUES2-Companion  
CLUES-Companion pipeline – **v2025-04-16**

*A streamlined pipeline for estimating, visualizing, and dating selection coefficients using Relate and CLUES v2 for multiple SNPs*

---

**Overview**

**CLUES2Companion.sh** is a comprehensive three-phase Bash pipeline that now includes:

- **Phase-1**:  
  - Convert VCF → `.haps` & `.sample`  
  - Prepare input (masking, flipping, filtering)  
  - Run Relate (two passes → `.anc`, `.mut`, `.coal`)  
  - Extract SNPs & polarize derived alleles → per-SNP files & frequency table  
  - **Automatic cleanup** of intermediate files  

- **Phase-2**:  
  - Sample branch lengths → Newick trees  
  - Run **RelateToCLUES.py** → per-SNP `*_times.txt`  
  - Run **CLUES v2** (`inference.py`) → per-SNP `_inference.txt` (+ optional `_CI.txt`)  
  - **Merge** all SNP inferences → combined TSV  
  - **Integrated plotting** in Python (error bars, −log₁₀(p), significance stars)  

- **Phase-3**:  
  - **Sliding-window dating** of a target rsID → initial onset estimate (generations & years)  
  - **Optional bootstrap** with user-defined settings → confidence intervals  
  - **JSON summary** of dating + CI  

Outputs per Phase-1/2/3 → output_C2Companion/phase{1,2,3}/<PREFIX>_chr<CHR>/.
All phases are driven through a single menu.  Logs are written under `output_C2Companion/log/`.
In addition, the script relies on *relative paths*.  
Please keep the following folders exactly where they are.

Moving or renaming **any** of these directories – or this script itself –
will break Phase1, Phase2 and Phase-3.
---

## Table of Contents

- [Installation](#installation)  
- [Dependencies](#dependencies)  
- [Input files & folder layout](#inputs)  
- [Quick start](#quick-start)  
- [Phase-1 (details)](#phase-1)  
- [Phase-2 (details)](#phase-2)  
- [Phase-3 (details)](#phase-3)  
- [Tips & troubleshooting](#tips)  
- [License](#license)  
- [Contact](#contact)  

---

<a name="installation"></a>
## Installation

```bash
git clone https://github.com/alisi1989/CLUES2-Companion.git
cd CLUES-Companion
chmod +x CLUES2Companion.sh
```
Ensure your working dir contains sub-folders:

```
Relate/     ← Relate binaries  
CLUES2/     ← CLUES v2 scripts  
required_files/  
  ├ ancestor/homo_sapiens_ancestor_chrN.fa  
  ├ mask/PilotMask_chrN.fasta  
  └ map/genetic_map_chrN.txt  
```

<a name="dependencies"></a>

## Dependencies

Bash ≥ 4.0 (for mapfile)
Python 3.8+
Relate v1.2.2
CLUES v2 (master branch)
GNU parallel (optional)
Python packages:

```
pip3 install cyvcf2 numpy pandas matplotlib adjustText
```

<a name="inputs"></a>

Input files & folder layout

```
CLUES-Companion/
├── CLUES2Companion.sh         # main driver
├── Relate/                    # Relate bin/, scripts/
├── CLUES2/                    # CLUES v2 scripts
├── required_files/            # reference data
│   ├ ancestor/…
│   ├ mask/…
│   └ map/…
└── your_data/
    ├ POP_chrN.vcf.gz          # phased & indexed VCF
    └ POP.poplabels            # sampleID, population, group, SEX
```
The user must use a fully phased vcf and the index file (*.tbi). 
The vcf file must contain only one population and one chromosome (e.g see example/FIN_chr2.vcf.gz and example/FIN_chr2.vcf.gz.tbi)

The *.poplabels file must contain 4 columns (sampleID, population, group, SEX). 

e.g. for Diploid organisms:
```
sample population group sex
UNR1 FIN EUR NA
UNR2 FIN EUR NA
UNR3 FIN EUR NA
UNR4 FIN EUR NA
```
for more information please refer to: `https://myersgroup.github.io/relate/input_data.html#Prepare`

---

<a name="quick-start"></a>

## Quick start:
```
./CLUES2Companion.sh
```

Menu prompts: choose 1, 2 or 3.


<a name="phase-1"></a>

## Phase-1 – Relate & SNP preparation

1 - Convert VCF → .haps, .sample \
2 - PrepareInputFiles (mask, flip, filter) \
3 - Relate mode All (run-1 → random Ne) → .anc, .mut \
4 - EstimatePopulationSize → .coal \
5 - Relate mode All (run-2 with --coal) → final .anc, .mut \
6 - Extract SNPs via cyvcf2 within user region \
7 - Polarize derived alleles + compute frequency → ${PREFIX}_Derived_<rs>.txt & ${PREFIX}_Frequency_chr<CHR>_<start>_<end>.txt \
8 - Cleanup of .rate, .pairwise.*, .annot, intermediate .haps/.sample, run-1 files. \

## Example of usage for Phase-1

```
./CLUES2Companion.sh
******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******
Choose phase to run
  1) Phase-1  : Relate (*.mut, *.anc, *.coal files) and SNP/Derived/DAF
  2) Phase-2  : Relate (BranchLengths) → RelateToCLUES.py → inference.py → outputs
  3) Phase-3  : Dating target SNP(s)
Enter option (1/2/3):
```
Here the user has to make his choice. There are 3 options available: Phase-1, Phase-2 and Phase-3.

`then`

```
Enter option (1/2/3): 1   


          ╔═════════════════════════════════════════════════════════╗
          ║        🚀  PHASE 1: RELATE & SNP EXTRACTION  🚀        ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose the chromosome to analyze (e.g. 2, 17, X): 2
Prefix of phased vcf/bcf file (e.g. example/FIN without _chrN.vcf.gz): example/FIN
Path to population‑labels file (*.poplabels)(e.g. example/FIN.poplabels): example/FIN.poplabels
Start bp of target region: 135839626
End bp of target region: 135876443
Prefix of output name (e.g. Finnish): FIN_MCM6
```

In the example above for Phase-1, the user is asked to specify the chromosome to be analyzed, the prefix of the vcf file without `_chrN` and without the `.vcf.gz` extension.
Then the user is asked to specify the start and end genomic position of the target region from which to extract the snps and the information present in the fully phased vcf.
Finally the user has to choose a prefix name for the output which will then be the prefix name to be used for the next phases.

For each step executed by Phase-1 the user will see a `INFO` message about the location of the output files generated by Phase-1, in addition to this the user will see the progress of each step inside Phase-1 (any errors or warnings and the message `"done"` for completed steps).

Upon completion of Phase-1, the user will receive an INFO about the number of SNPs, the location and number of completed *_derived.txt files and the location of the file with the frequencies of the snps)

<a name="phase-2"></a>

### Phase-2 – Selection coefficient inference

1 - SampleBranchLengths → Newick trees \
2 - RelateToCLUES.py → *_times.txt \
3 - inference.py (CLUES v2) → _inference.txt (+ *_CI.txt if requested) \
4 - Merge into <PREFIX>_merged_inference_chr<CHR>.tsv \
5 - Plot \

---

## Example of usage for Phase-2

```
******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******
Choose phase to run
  1) Phase-1  : Relate (*.mut, *.anc, *.coal files) and SNP/Derived/DAF
  2) Phase-2  : Relate (BranchLengths) → RelateToCLUES.py → inference.py → outputs
  3) Phase-3  : Dating target SNP(s)
Enter option (1/2/3): 2

```

`then`

```
******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******
Choose phase to run
  1) Phase-1  : Relate (*.mut, *.anc, *.coal files) and SNP/Derived/DAF
  2) Phase-2  : Relate (BranchLengths) → RelateToCLUES.py → inference.py → outputs
  3) Phase-3  : Dating target SNP(s)
Enter option (1/2/3): 2


          ╔═════════════════════════════════════════════════════════╗
          ║        🧬  PHASE-2 – (requires Phase‑1 outputs)  🧬     ║
          ║   Please read the manual carefully before proceeding    ║
          ╚═════════════════════════════════════════════════════════╝

Choose the chromosome to analyze (e.g. 2, 17, X): 2
Enter population prefix used in Phase‑1 (e.g. Finnish): FIN_MCM6
Phase‑1 auto-detect directory [ENTER = /Users/alessandrolisi1989/desktop/CLUES2Companion2/output_C2Companion/phase1/FIN_MCM6_chr2] or provide a different folder with phase1 outputs :

→ Using frequency file: /Users/alessandrolisi1989/desktop/CLUES2Companion2/output_C2Companion/phase1/FIN_MCM6_chr2/FIN_MCM6_Frequency_chr2_135839626_135876443.txt
→ Using SNPs file: /Users/alessandrolisi1989/desktop/CLUES2Companion2/output_C2Companion/phase1/FIN_MCM6_chr2/FIN_MCM6_SNPs_chr2_135839626_135876443.txt

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

### The script prompts:
1 - Select the chromosome:
Chromosome to process (e.g. 2, 17, X):
If you ran Phase-1 on several chromosomes, the pipeline automatically looks for the Phase-1 output that matches the chromosome you enter.

### 2 - Enter the population prefix:
Provide the same prefix you used in Phase-1 (e.g. Finnish).
The script now knows the exact folder structure
output_C2Companion/phase1/<prefix>_chr<chr>/.

### 3 - Confirm detected Phase-1 folders:
The script shows the paths it has found ( .anc, .mut, .coal, SNP list, frequency table ).
Press Enter to accept.
If you manually moved Phase-1 results (not recommended), type the new directory and press Enter.

### 4 - Mandatory CLUES-2 parameters: (please see CLUES2 and Relate manual, see the links below)

| prompt                                    | description                                          | reference              |
| ----------------------------------------- | ---------------------------------------------------- | ---------------------- |
| `tCutoff (e.g. 1000):`                    | analyse selection back to this many generations      | CLUES-2 `--tCutoff`    |
| `df (e.g. 600):`                          | number of degrees of freedom in the Gaussian process | CLUES-2 `--df`         |
| `Importance-sampling for branch lengths:` | how many trees to sample from *Relate*               | Relate `--num_samples` |

None of these can be left blank.

### 5 - Optional files: (please see CLUES2 and Relate manual, see the links below)

AncientSamps — table of sample ages
AncientHaps — ancient haplotypes
Press Enter to skip either file.
(See CLUES-2 README for exact format.)

### 6 - Optional run-time switches" (please see CLUES2 and Relate manual, see the links below)

| prompt                                  | effect                                                              |
| --------------------------------------- | ------------------------------------------------------------------- |
| `Disable allele trajectory? (y/N):`     | answer **y** to *skip* posterior trajectories – speeds up inference |
| `Confidence interval (0-1, e.g. 0.95):` | compute CI for the selection MLE; **Enter** = no CI                 |
| `Time-bins (e.g. 200 300):`             | split 0-*tCutoff* into custom intervals; *Enter* = single epoch     |

Example 200 300 on a tCutoff 500 gives three epochs:
0-200, 200-300, 300-500.

### 7 - Internal bookkeeping:

NOTE: All output files will be local to `~/Output_C2Companion/phase2/{prefix}_chr{N}/`
(e.g in this example output are located in `~/Output_C2Companion/phase2/FIN_MCM6_original_chr2/)`

### In this output folder the user will find:

### 1 - The plot made with:

All the snps analyzed in Phase-2 
Bars from CI columns \ 
Intensity color bar = −log₁₀(p) \ 
Significance stars based on the p-value above each SNP (* (0.05), ** (0.01), *** (0.01)) 


### 2 - An Excell chart table file with all the SNPs and the related statistics calculated by CLUES2. 

**Example of merged chart table**
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
[INFO] In addition to the final outputs, the user will find 3 folders (`{prefix_inference_chr{N}/` `prefix_times_chr{N}/` `prefix_trees_chr{N}`/ in which are stored respectively: the outputs of inference.py for each snp, the *_times.txt files of RelateToClues.py for each snps, and the *.newick files for each snps). 

**Further reading**

CLUES-2 manual: https://github.com/avaughn271/CLUES2#command-line-arguments \
Relate SampleBranchLengths: https://myersgroup.github.io/relate/modules.html#SampleBranchLengths \

---

<a name="phase-3"></a>

Phase-3 – Dating selective sweeps

Initial scan over fixed windows (0–2000 gen) → onset estimate
Save initial JSON:

`{ "rsID": "...", "onset_gen": 100, "onset_years": 2800, … }`


Optional bootstrap: user-defined epochs, window size, replicates
Compute per-replicate onsets → summary statistics (median, 95% CI)
Save JSON with CI:
rsID_Boostraps_onset_Dating+CI.json

<a name="tips"></a>

Tips & troubleshooting

Logging: check output_C2Companion/log/ for debug
Missing deps: verify RelateFileFormats, PrepareInputFiles.sh, Python packages
Spinner issues: ensure bash -euo pipefail is supported
Plot errors: install adjustText for label adjustment
Bootstrap: skip by answering “N” when prompted
<a name="license"></a>

License

MIT License

<a name="contact"></a>

Contact

Alessandro Lisi alisi@usc.edu (alisi1989)
Michael C. Campbell mc44680@usc.edu (mc44680)
