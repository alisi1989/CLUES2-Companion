# CLUES2-Companion
CLUES‑Companion pipeline

# CLUES‑Companion  
*A streamlined pipeline for estimating and visualising single‑locus
selection coefficients with Relate + CLUES v2*

---

**Overview**

**CLUES‑Companion** is a two‑phase command‑line workflow that

1. runs **Relate** to build local genealogies and sample branch lengths  
2. feeds those branch lengths to **CLUES v2** to infer the selection
   coefficient (*s*) for every SNP in a user‑defined region  
3. merges the per‑SNP results into a single table and generates
   publication‑ready plots (selection coefficient ± CI,
   colour = −log<sub>10</sub> *p*, asterisks for significance, etc.).

Everything is controlled from one Bash script
`clues_pipeline.sh`; an optional Python helper (`plot_CLUES2.py`)
creates the final figure.

The code has been tested on macOS 14 and Ubuntu 22 with
**Relate v1.3.3** and the current master branch of **CLUES v2**.

---

## Table of Contents

- [Installation](#installation)
- [Dependencies](#dependencies)
- [Input files & folder layout](#inputs)
- [Quick start](#quick-start)
- [Phase‑1 (details)](#phase-1)
- [Phase‑2 (details)](#phase-2)
- [Plotting script](#plotting-script)
- [Tips & troubleshooting](#tips)
- [License](#license)
- [Contact](#contact)

---

<a name="installation"></a>
## Installation

```bash
git clone https://github.com/your‑org/CLUES‑Companion.git
cd CLUES‑Companion
chmod +x clues_pipeline.sh                # main driver
conda env create -f env.yml               # optional: creates py env with cyvcf2 etc.
```
Relate must be present as the sub‑folder Relate/ (pre‑compiled binaries work fine).
If you already have Relate elsewhere, simply symlink it:
ln -s /path/to/Relate ./Relate.

<a name="dependencies"></a> ## Dependencies


Type	Tool	Tested version
Core	Relate	v1.3.3
      	CLUES v2 (inference.py)	commit bd5d2bc
Bash	GNU parallel (optional)	20240222
Python	cyvcf2, numpy, pandas, matplotlib	see env.yml

<a name="inputs"></a> ## Input files & folder layout
CLUES‑Companion/
├── Relate/                  ←  Relate executable     (fixed path)
├── mandatory/               ←  helper references
│   ├── ancestor/homo_sapiens_ancestor_chrN.fa
│   ├── mask/PilotMask_chrN.fasta
│   └── map/genetic_map_chrN.txt
├── your_data/
│   ├── POP_chrN.vcf.gz      ←  phased & indexed VCF
│   └── POP.poplabels        ←  2‑column pop‑label file
└── clues_pipeline.sh

> VCFs must be phased, gzipped and tabix‑indexed.
> Chromosomes should be chr1 … chr22 (or 1 … 22) consistently.
> poplabels → two columns: sample ID & population.

 <a name="quick-start"></a> ## Quick start (Chromosome 2, gene MCM6)
 # Phase‑1   (Relate → .anc/.mut → .coal → SNP extraction …)
./clues_pipeline.sh
# choose option 1 and follow prompts
#
#   ‑ chromosome: 2
#   ‑ VCF prefix:   your_data/Douz
#   ‑ poplabels:    your_data/Douz.poplabels
#   ‑ region start / end (bp): 135 800 000 – 135 860 000
#   ‑ output prefix: DOUZ
#
# Outputs land in   working‑dir/phase1/

# Phase‑2   (branch lengths → CLUES → merge + plot)
./clues_pipeline.sh
# choose option 2 and answer:
#
#   ‑ chromosome: 2
#   ‑ population prefix: DOUZ
#   ‑ accept auto‑detected Frequency & SNP files
#   ‑ Relate prefix: working‑dir/phase1/DOUZ_GS_chr2
#   ‑ CI: 0.98
#
# Final table:  working‑dir/phase2/DOUZ_merged_inference_chr2.tsv
# Figure:       working‑dir/phase2/DOUZ_singleplot_ci.pdf

<a name="phase-1"></a> ## Phase‑1 – Relate pipeline + SNP / freq files

Convert‑from‑VCF → *.haps + *.sample
PrepareInputFiles (adds ancestor, mask, dist)
Relate (run‑1, N = 30 000)
EstimatePopulationSize → *.coal
Relate (run‑2, with --coal) → final *.anc/*.mut
SNP extraction (cyvcf2 slice of user region)
Derived‑allele vectors from .mut + .haps → one ${prefix}_Derived_<rs>.txt per SNP
Frequency table ${prefix}_Frequency_chrN_start_end.txt
Intermediate run‑1 files, pairwise .coal and big temporary *.rate files are deleted automatically to save space.

<a name="phase-2"></a> ## Phase‑2 – branch lengths → CLUES inference

SampleBranchLengths (Relate --format n, 200 samples)
RelateToCLUES.py → <rs>_times.txt
CLUES v2 (inference.py) with options set by user
optional arguments:
--CI, --ancientSamps, --ancientHaps, --noAlleleTraj
Per‑SNP outputs are merged into a TSV; if --CI was requested, the lower / upper bounds are appended.

<a name="plotting-script"></a> ## Plotting script

```
python plot_CLUES2.py  working-dir/phase2/DOUZ_merged_inference_chr2.tsv \
       -o DOUZ_chr2_MCM6
```

Creates:

DOUZ_chr2_MCM6_singleplot_ci.pdf/png – selection coefficient with
– CI bars,
– colour = −log<sub>10</sub>(P),
– asterisks for significance (* < 0.05, ** < 0.01, *** < 0.001)
– rsID labels.
All thresholds and colours can be tuned via CLI flags (--p1 0.05 --p2 0.01 …).
See python plot_CLUES2.py --help.

<a name="contact"></a> ## Contact

Alessandro Lisi  alisi@usc.edu   (alisi1989)
Michael C. Campbell  mc44680@usc.edu   (mc44680)
