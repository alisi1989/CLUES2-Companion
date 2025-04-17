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

The code has been tested on macOS 14, Linux and Ubuntu 22 with
**Relate v1.2.2** and the current master branch of **CLUES v2 (CLUES2)**.

---

 ## Table of Contents

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

## Installation

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
```
CLUES‑Companion/
├── Relate/                  ←  Relate executable (bin/, script/)
├── mandatory/               ←  mandatory files already formatted
│   ├── ancestor/homo_sapiens_ancestor_chrN.fa
│   ├── mask/PilotMask_chrN.fasta
│   └── map/genetic_map_chrN.txt
├── your_data/
│   ├── POP_chrN.vcf.gz      ←  phased & indexed VCF
│   └── POP.poplabels        ←  4‑column pop‑label file
└── ./CLUES2Companion.sh
```

- VCFs must be phased, gzipped and tabix‑indexed.
- Chromosomes should be chr1 … chr22 (or 1 … 22) consistently.
- poplabels → four columns: sample ID, population, group and SEX. for more information please visit "https://myersgroup.github.io/relate/input_data.html"

 <a name="quick-start"></a> 
 
 ## Quick start (Chromosome 2, example gene MCM6)

 ## Phase‑1   (Relate → .anc/.mut → .coal → SNP extraction …)

example:
```
./CLUES2Companion.sh

******  CLUES‑Companion – please cite CLUES2 and this helper  ******
Choose phase to run
  1) Phase‑1  : Relate + SNP/Derived/DAF
  2) Phase‑2  : BranchLengths → CLUES → merge   (requires Phase‑1 outputs)
Enter option (1/2): 1

→ Running phase1 … (Relate + SNP extraction + Derived/DAF)
Chromosome (e.g. 2, 17, X): 2
Prefix to phased VCF/BCF (without _chrN.vcf.gz): example/Finnish
Path to population‑labels file (.poplabels): .example/poplabels/Finnish.poplabels
Start bp of target region: 135839626
End   bp of target region: 135876443
Output base prefix (e.g. Finnish): Finnish
```
output dir `working-dir` and subfolder for each phase are automatically created and the outputs land automatically in `working‑dir/phase1/`

 ## Phase‑2   (branch lengths → CLUES → merge + plot)
 
 example:
```
./CLUES2Companion.sh
choose option 2 and answer:

‑ chromosome: 2
‑ population prefix: Finnish
‑ accept auto‑detected Frequency & SNP files
```

Final table:  working‑dir/phase2/DOUZ_merged_inference_chr2.tsv
Figure:       working‑dir/phase2/DOUZ_singleplot_ci.pdf

output dir `working-dir` and subfolder for each phase are automatically created and the outputs land automatically in `working‑dir/phase2/`

<a name="phase-1"></a> ## Phase‑1 – Relate pipeline + SNP / freq files

```
1 → Convert‑from‑VCF → *.haps + *.sample
2 → PrepareInputFiles (flipping snps if necessary using ancestor reference alleles, filtering snps, if necessary, using pilot mask allele, removing non-biallelic snps)
3 → Relate `--mode All` (run‑1, with random N = 30000)
4 → EstimatePopulationSize → to create the *.coal file (coalescence file that serve as input for .anc and .mut file into next step)
5 → Relate `--mode All` (run‑2, with --coal flag specified) → final *.anc/*.mut based on coalescence rates
6 → SNP extraction (cyvcf2 slice of user region)
7 → Derived‑allele polarization using .mut (that contains information about polarization) + .haps.gz (that contains allele coded as 0 and 1 non polarized) → one ${prefix}_Derived_<rs>.txt per SNP
8 → Frequency table of the pre-calculated derived alleles →  ${prefix}_Frequency_chrN_start_end.txt
9 → Intermediate run‑1 files, pairwise .coal and big temporary *.rate files are deleted automatically to save space.
```

<a name="phase-2"></a> ## Phase‑2 – branch lengths → CLUES inference

```
1 → SampleBranchLengths (Relate --format n, 200 importance-samples (upper limit to use)
1 → RelateToCLUES.py → population_<rs>_times.txt
1 → CLUES v2 (inference.py) with options set by user
1 → optional arguments:
1 → --CI, --ancientSamps, --ancientHaps, --noAlleleTraj
1 → Per‑SNP outputs are merged into a TSV; if --CI was requested, the lower / upper bounds are appended.
```

<a name="plotting-script"></a> ## Plotting script

example:
```
python plot_CLUES2.py  working-dir/phase2/Finnish_merged_inference_chr2.tsv \
       -o Finnish_chr2_MCM6 (optional --db)
```
the optional flag `--db`, if provided, apply the correction for multiple test

the `plot_CLUES2.py creates:

`DOUZ_chr2_MCM6_singleplot_ci.pdf/png – selection coefficient with`

where are reported:
`CI bars(if specified in phase2)`
`Bar intenisity colour = −log<sub>10</sub>(P)`
`Asterisks for significance (* < 0.05, ** < 0.01, *** < 0.001)`
`rsID labels`
`All thresholds and colours can be tuned via CLI flags (--p1 0.05 --p2 0.01 …)`
`See python plot_CLUES2.py --help`

<a name="contact"></a> ## Contact

Alessandro Lisi  alisi@usc.edu   (alisi1989)
Michael C. Campbell  mc44680@usc.edu   (mc44680)
