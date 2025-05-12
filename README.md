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

All phases are driven through a single menu.  Logs are written under `output_C2Companion/log/`.

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
git clone https://github.com/your-org/CLUES-Companion.git
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

To symlink existing installations:
```
ln -s /path/to/Relate   ./Relate
ln -s /path/to/CLUES2   ./CLUES2
```

<a name="dependencies"></a>

Dependencies

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

<a name="quick-start"></a>

Quick start
```
./CLUES2Companion.sh
```


Menu prompts: choose 1, 2 or 3.
Logs → output_C2Companion/log/.
Outputs per phase → output_C2Companion/phase{1,2,3}/<PREFIX>_chr<CHR>/.

<a name="phase-1"></a>

Phase-1 – Relate & SNP preparation

Convert VCF → .haps, .sample
PrepareInputFiles (mask, flip, filter)
Relate mode All (run-1 → random Ne) → .anc, .mut
EstimatePopulationSize → .coal
Relate mode All (run-2 with --coal) → final .anc, .mut
Extract SNPs via cyvcf2 within user region
Polarize derived alleles + compute frequency →
${PREFIX}_Derived_<rs>.txt
${PREFIX}_Frequency_chr<CHR>_<start>_<end>.txt
Cleanup of .rate, .pairwise.*, .annot, intermediate .haps/.sample, run-1 files.

<a name="phase-2"></a>

Phase-2 – Selection coefficient inference

SampleBranchLengths → Newick trees
RelateToCLUES.py → *_times.txt
inference.py (CLUES v2) → _inference.txt (+ *_CI.txt if requested)
Merge into <PREFIX>_merged_inference_chr<CHR>.tsv
Python plot:
Error bars from CI columns
Color = −log₁₀(p)
Significance stars (*, **, ***)
Outputs: *_singleplot_ci.pdf & .png

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
