#!/usr/bin/env bash

# ------------------------------------------------------------------
#  CLUES2Companion pipeline
#  • Phase‑1  : Relate + SNP extraction + Derived files + frequency
#  • Phase‑2  : BranchLengths → RelateToCLUES → inference → merge
#  • Phase‑3  : Dating target SNP(s)
# ------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# ---------- helpers ----------------------------------------------------------
abs_path() {
  local p="$1"
  [[ -z "$p" ]] && { echo ""; return; }
  case "$p" in
    /*) echo "$p" ;;
    ~*) echo "${p/#\~/${HOME}}" ;;
     *) echo "$(pwd)/$p" ;;
  esac
}

SCRIPT_NAME=$(basename "$0")
BASE_DIR=$(pwd)
WORK_BASE_RELATE_DEFAULT="${BASE_DIR}/output_CLUES2Companion-Relate"
WORK_BASE_SINGER_DEFAULT="${BASE_DIR}/output_CLUES2Companion-Singer"

LOG_ROOT="${BASE_DIR}/logs"
mkdir -p "$LOG_ROOT"
LOG_FILE="${LOG_ROOT}/${SCRIPT_NAME%.*}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# colours
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  BOLD="\e[1m"; GREEN="\e[32m"; YELLOW="\e[33m"; RED="\e[31m"; RESET="\e[0m"
else
  BOLD=""; GREEN=""; YELLOW=""; RED=""; RESET=""
fi
info()  { echo -e "${BOLD}${GREEN}→ $*${RESET}"; }
warn()  { echo -e "${YELLOW}[WARN] $*${RESET}"; }
error() { echo -e "${RED}[ERROR] $*${RESET}"; exit 1; }

mark_done(){ touch "$1/.done"; }
phase_done(){ [[ -f "$1/.done" ]]; }

# ---------- menu -------------------------------------------------------------
echo
echo "******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******"
cat <<EOF
Choose phase to run
  1) Phase 1  : Apply Relate/SINGER (*.mut, *.anc, *.coal, and derived allele frequency file for Relate or *.trees files for SINGER)
  2) Phase 2  : Apply Relate or SINGER, and then run CLUES2 to infer selection coefficient (s) along with table and figure
  3) Phase 3  : Date onset of selective sweeps on target SNP(s) using Relate
EOF
read -rp "Enter option (1/2/3): " OPTION
echo

# ---------- PHASE 1  helpers ---------------------------------------------------------

run_phase1_relate() {
  PH="phase1"
  ENGINE="Relate"
  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║         🚀  PHASE 1: RELATE & SNP EXTRACTION  🚀        ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"

  # 1) Input utente
  read -rp  "Choose the chromosome to analyze (e.g. 2, 17, X): " chr
  read -e -rp "Prefix of phased vcf/bcf file (e.g. example/Finnish without _chrN.vcf.gz): " vcf_prefix_in
  read -e -rp "Path to population-labels file (*.poplabels)(e.g. example/Finnish.poplabels): "  pop_labels_in
  read -rp  "Start bp of target region: " start_bp
  read -rp  "End bp of target region: " end_bp
  read -e -rp "Prefix of output name (e.g. Finnish): " out_prefix
  export OUT_PREFIX="$out_prefix"
  
  WORK_BASE="$WORK_BASE_RELATE_DEFAULT"
  RUN_ID="${out_prefix}_chr${chr}"
  WORK_DIR="${WORK_BASE}/${PH}/${RUN_ID}"
  mkdir -p "$WORK_DIR"
  echo
  echo -e "[INFO] Working dir: $WORK_DIR"

  # 2) Percorsi Relate e required_files
  relate_root="${BASE_DIR}/Relate-Linux"
  [[ -x "${relate_root}/bin/Relate" ]] || error "Relate binary not found in ${relate_root}/bin"

  required_files_dir="${BASE_DIR}/required_files"
  ancestral="${required_files_dir}/ancestor/homo_sapiens_ancestor_chr${chr}.fa"
  mask="${required_files_dir}/mask/PilotMask_chr${chr}.fasta"
  genmap="${required_files_dir}/map/genetic_map_chr${chr}.txt"
  for f in "$ancestral" "$mask" "$genmap"; do [[ -f $f ]] || error "Missing $f"; done
  pop_labels=$(abs_path "$pop_labels_in")

  # 3) Prefix interni
  vcf_gz="${vcf_prefix_in}_chr${chr}"   # senza estensione
  haps_path="${WORK_DIR}/${out_prefix}_chr${chr}.haps"
  sample_path="${WORK_DIR}/${out_prefix}_chr${chr}.sample"
  pif_base="${WORK_DIR}/${out_prefix}_PIF_chr${chr}"
  gs_base="${WORK_DIR}/${out_prefix}_GS_chr${chr}"
  gs_coal="${WORK_DIR}/${out_prefix}_GS+COAL_chr${chr}"
  coal_base="${WORK_DIR}/${out_prefix}_EPS4COAL_chr${chr}"

  # 4) Convert VCF → haps/sample
  echo
  info "Convert vcf → *.haps and *.sample files"
  "${relate_root}/bin/RelateFileFormats" --mode ConvertFromVcf \
      -i "$vcf_gz" \
      --haps "$haps_path" --sample "$sample_path" \
      -o "${WORK_DIR}/${out_prefix}_chr${chr}" 2>/dev/null \
      || error "RelateFileFormats failed"
  echo "✓ Convert vcf file done"

  # 5) PrepareInputFiles
  echo
  info "Preparing input files → *.haps.gz, *.sample.gz, *.dist, *.annot"
  rand_seed() { echo $(( ( RANDOM << 15 ) | RANDOM )); }
  "${relate_root}/scripts/PrepareInputFiles/PrepareInputFiles.sh" \
      --haps      "$haps_path" \
      --sample    "$sample_path" \
      --ancestor  "$ancestral" \
      --mask      "$mask" \
      --poplabels "$pop_labels" \
      --seed "$(rand_seed)" \
      -o          "$pif_base" 2>/dev/null \
      || error "PrepareInputFiles failed"
  echo "✓ Preparing input files done"

  pif_haps="${pif_base}.haps.gz"
  pif_sample="${pif_base}.sample.gz"
  pif_annot="${pif_base}.annot"

  # 6) Relate run-1 (MCMC branch lengths iniziali)
  echo
  info "Running Relate mode All with random Ne → *.anc and *.mut files"
  pushd "$WORK_DIR" >/dev/null
  "${relate_root}/bin/Relate" --mode All \
      --haps  "$pif_haps" --sample "$pif_sample" \
      --map   "$genmap"   --annot "$pif_annot" \
      -N 30000 -m 1.25e-8 -o GS_tmp --seed "$(rand_seed)" 2>/dev/null \
      || { popd >/dev/null; error "Relate --mode All failed"; }
  for f in "$WORK_DIR"/GS_tmp.*; do mv "$f" "${gs_base}_run1${f##*/GS_tmp}"; done
  echo "✓ Relate mode All with random Ne done"
  popd >/dev/null

  # 7) EstimatePopulationSize → .coal
  echo
  info "Estimating population size with Relate → *.coal"
  pushd "$WORK_DIR" >/dev/null
  if ! "${relate_root}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh" \
          -i "${gs_base}_run1" \
          -o "${coal_base}" \
          --noanc 1 --poplabels "$pop_labels" \
          -m 1.25e-8 --years_per_gen 28 --seed "$(rand_seed)" 2>/dev/null; then
    # molte installazioni falliscono solo nella fase di plot R/ggplot: continuiamo
    warn "EstimatePopulationSize completed without plotting (likely missing R/ggplot) – continuing."
  fi
  popd >/dev/null
  echo "✓ Estimating population size done"

  # 8) Re-estimate branch lengths con la .coal
  echo
  info "Re-estimating branch lengths (MCMC) → final *.anc/*.mut"
  pushd "$WORK_DIR" >/dev/null
  if ! "${relate_root}/scripts/SampleBranchLengths/ReEstimateBranchLengths.sh" \
        -i  "${gs_base}_run1" \
        -o  "${gs_coal}" \
        -m  1.25e-8 \
        --coal "${coal_base}.coal" \
        --seed "$(rand_seed)" 2>/dev/null; then
    popd >/dev/null
    error "ReEstimateBranchLengths failed – check log."
  fi
  popd >/dev/null
  echo "✓ Branch lengths re-estimated"

  # Esponi percorsi finali per gli step comuni
  mut_file="${gs_coal}.mut"
  haps_file="${pif_base}.haps.gz"
}

run_phase1_singer() {
  PH="phase1"
  ENGINE="Singer"
  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║     🧬  PHASE 1: SINGER – ARGs INFERENCE FROM VCF  🧬     ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"

  # --- Input parameters ---
  read -rp  "Choose the chromosome to analyze (e.g. 2, 17, X): " chr
  read -e -rp "Prefix of phased VCF file (without _chrN.vcf): " vcf_prefix_in
  read -e -rp "Output prefix name (e.g. Bedouin_MCM6): " out_prefix
  read -rp  "Start bp of target region: " start_bp
  read -rp  "End bp of target region: " end_bp
  read -rp  "Mutation rate (-m, e.g. 1.25e-8): " mu
  read -rp  "Effective population size (-Ne) [Optional. press ENTER to skipp]: " Ne
  read -rp  "Recombination/mutation ratio (-ratio) [Default: 1, press ENTER to skipp]: " ratio
  read -e -rp "Recombination map file (-recomb_map, press ENTER to skipp): " recomb_map
  read -e -rp "Mutation rate map file (-mut_map) [Optional, press ENTER to skipp]: " mut_map
  read -rp  "Number of MCMC samples (-n) [Default: 100, press ENTER to skipp]: " n_mcmc
  # Normalizza n_mcmc e rendilo disponibile fuori
if [[ -z "$n_mcmc" ]]; then
  N_MCMC=100
elif [[ "$n_mcmc" =~ ^[0-9]+$ ]]; then
  N_MCMC="$n_mcmc"
else
  warn "Invalid -n MCMC samples '$n_mcmc'; falling back to 100"
  N_MCMC=100
fi
export N_MCMC

  read -rp  "Thinning interval (-thin) [Default: 20]: " thin
  read -rp  "Site flip probability (-polar) [Default: 0.5]: " polar
  read -rp  "Random seed (-seed) [Default: 42]: " seed
  export OUT_PREFIX="$out_prefix"

  WORK_BASE="$WORK_BASE_SINGER_DEFAULT"
  RUN_ID="${out_prefix}_chr${chr}"
  WORK_DIR="${WORK_BASE}/${PH}/${RUN_ID}"
  mkdir -p "$WORK_DIR"
  echo
  echo -e "[INFO] Working dir: $WORK_DIR"

  # --- Run SINGER ---
  info "Running SINGER ARG inference..."

  singer_exec="${BASE_DIR}/Singer-Linux/singer_master"

  [[ -x "$singer_exec" ]] || error "SINGER executable not found at $singer_exec"

  #vcf_gz="${vcf_prefix_in}_chr${chr}"
  vcf_prefix="${vcf_prefix_in}_chr${chr}"
  vcf_plain="${vcf_prefix}"
  vcf_gz="${vcf_prefix}.vcf.gz"   # servirà dopo per cyvcf2 (estrazione SNP)


  # Build command
  singer_cmd="$singer_exec -vcf $vcf_plain -output $WORK_DIR/$out_prefix -m $mu -start $start_bp -end $end_bp"
 
  [[ -n "$Ne" ]] && singer_cmd="$singer_cmd -Ne $Ne"
  [[ -n "$recomb_map" ]] && singer_cmd="$singer_cmd -recomb_map $recomb_map"
  [[ -n "$mut_map" ]]    && singer_cmd="$singer_cmd -mut_map $mut_map"
  [[ -n "$ratio" ]]      && singer_cmd="$singer_cmd -ratio $ratio"
  [[ -n "$n_mcmc" ]]     && singer_cmd="$singer_cmd -n $n_mcmc"
  [[ -n "$thin" ]]       && singer_cmd="$singer_cmd -thin $thin"
  [[ -n "$polar" ]]      && singer_cmd="$singer_cmd -polar $polar"
  [[ -n "$seed" ]]       && singer_cmd="$singer_cmd -seed $seed"

  echo -e "[COMMAND] $singer_cmd"
  eval "$singer_cmd"

  mark_done "$WORK_DIR"
  echo
  info "✓ SINGER inference completed — outputs in ${WORK_DIR}"
  printf "\n"
}
# ---------- PHASE 1  --------------------------------------------------------

phase1() {
  PH="phase1"
  echo
  echo "Which method do you want to use for genealogical inference?"
  echo "  1) Relate"
  echo "  2) SINGER"
  read -rp "Enter option (1/2): " PHASE1_METHOD

  case "$PHASE1_METHOD" in
    1) run_phase1_relate ;;
    2) run_phase1_singer ;;
    *) error "Invalid method selection" ;;
  esac

# Dopo la scelta del metodo, abbiamo vcf_prefix definito nei run_* (variabile globale in bash se non "local")
  vcf_prefix="${vcf_prefix_in}_chr${chr}"   # ridichiariamo esplicitamente per chiarezza
  VCF_GZ_PATH="${vcf_prefix}.vcf.gz"
  [[ -f "$VCF_GZ_PATH" ]] || error "Missing compressed VCF: $VCF_GZ_PATH (bgzip + tabix required)"

  # -------------------------------------------------------------------------
  # 9. SNP extraction
echo
info "Extracting SNPs using cyvcf2"

export OUT_PREFIX="$OUT_PREFIX"
export OUTPUT_DIR="$WORK_DIR" CHR="$chr" VCF_GZ="$VCF_GZ_PATH"
export START_BP="$start_bp" END_BP="$end_bp"
snp_file="${WORK_DIR}/${OUT_PREFIX}_SNPs_chr${chr}_${start_bp}_${end_bp}.txt"
export SNP_FILE="$snp_file"

python - <<'PY'
import os
from cyvcf2 import VCF

vcf_path = os.environ["VCF_GZ"]
chr_in   = os.environ["CHR"].strip()
start_bp = int(os.environ["START_BP"])
end_bp   = int(os.environ["END_BP"])
outp     = os.environ["SNP_FILE"]

def norm_chr(s: str) -> str:
    s = s.strip()
    return s[3:] if s.lower().startswith("chr") else s

chr_norm = norm_chr(chr_in)

vcf = VCF(vcf_path)
cnt = 0
with open(outp, "w") as out:
    for v in vcf:
        if norm_chr(str(v.CHROM)) != chr_norm:
            continue
        pos = int(v.POS)
        if pos < start_bp or pos > end_bp:
            continue
        if v.ID and v.ID != ".":
            out.write(f"{v.ID}\t{pos}\n")
            cnt += 1

print(f"Extracted {cnt} SNPs → {outp}")
PY

[[ -s "$snp_file" ]] || error "No SNPs extracted."

  # -------------------------------------------------------------------------
  # 10. Derived alleles & frequency
  # -------------------------------------------------------------------------

if [[ "$PHASE1_METHOD" == "1" ]]; then
  echo
  info "Polarizing derived alleles and calculating derived allele frequency"
  export HAPS_FILE="$haps_file" MUT_FILE="$mut_file"
  freq_file="${WORK_DIR}/${OUT_PREFIX}_Frequency_chr${chr}_${start_bp}_${end_bp}.txt"
  export FREQ_FILE="$freq_file"

  python - <<'PY'
import os, gzip
op = lambda p,m="rt": gzip.open(p,m) if p.endswith(".gz") else open(p,m)
prefix=os.environ["OUT_PREFIX"]; outdir=os.environ["OUTPUT_DIR"]

# snp list
snps=[ln.split() for ln in open(os.environ["SNP_FILE"])]
# mut dict
mut={}
with op(os.environ["MUT_FILE"]) as f:
    next(f)
    for ln in f:
        fs=ln.strip().split(";")
        if len(fs)<11 or fs[3]==".":continue
        rs=fs[3]; a,d=fs[10].split("/")
        if fs[7]=="1": a,d=d,a
        mut[rs]=d.upper()
# haps dict
haps={}
with op(os.environ["HAPS_FILE"]) as f:
    for ln in f:
        p=ln.split();
        if len(p)<6: continue
        haps[p[1]]=(p[3].upper(),p[4].upper(),p[5:])

out=open(os.environ["FREQ_FILE"],"w")
out.write("rsID\tPOS\tallele_der\tfreq_der\n")
for rs,pos in snps:
    if rs not in mut or rs not in haps: continue
    der=mut[rs]; ref,alt,vec=haps[rs]
    bits=[int(b) if alt==der else 1-int(b) for b in vec]
    freq=sum(bits)/len(bits)
    dpath=f"{outdir}/{prefix}_Derived_{rs}.txt"
    with open(dpath,"w") as d: d.write("\n".join(map(str,bits)))
    out.write(f"{rs}\t{pos}\t{der}\t{freq:.4f}\n")
    print(f"  ↳ {dpath}")
print(f"✓ Frequency table → {os.environ['FREQ_FILE']}")
PY

else

# --- SINGER: usa solo il VCF per le frequenze ALT (no polarizzazione) ---
# Usando SEMPRE il .vcf.gz per cyvcf2
freq_file="${WORK_DIR}/${OUT_PREFIX}_Frequency_chr${chr}_${start_bp}_${end_bp}.txt"
export FREQ_FILE="$freq_file"
export VCF_GZ="$VCF_GZ_PATH"

python - <<'PY'
import os
from cyvcf2 import VCF

prefix   = os.environ["OUT_PREFIX"]
outdir   = os.environ["OUTPUT_DIR"]
vcf_path = os.environ["VCF_GZ"]
chr_in   = os.environ["CHR"].strip()
start_bp = int(os.environ["START_BP"])
end_bp   = int(os.environ["END_BP"])

def norm_chr(s: str) -> str:
    s = s.strip()
    return s[3:] if s.lower().startswith("chr") else s

chr_norm = norm_chr(chr_in)

# Carica lista SNPs (rsID -> pos)
snps_set = set()
pos_by_rs = {}
with open(os.environ["SNP_FILE"]) as f:
    for ln in f:
        rs, pos = ln.split()[:2]
        snps_set.add(rs)
        pos_by_rs[rs] = int(pos)

vcf = VCF(vcf_path)
out = open(os.environ["FREQ_FILE"], "w")
out.write("rsID\tPOS\tallele_ALT\tfreq_ALT\n")

n_written = 0
for v in vcf:
    vid = v.ID
    if not vid or vid == "." or vid not in snps_set:
        continue
    if norm_chr(str(v.CHROM)) != chr_norm:
        continue
    pos = int(v.POS)
    if pos < start_bp or pos > end_bp:
        continue

    gt = v.genotypes  # [[a1,a2,phased], ...]
    if not gt:
        continue

    alleles = []
    for a1, a2, *_ in gt:
        if a1 in (0,1): alleles.append(a1)
        if a2 in (0,1): alleles.append(a2)
    if not alleles:
        continue

    alt_freq = sum(1 for a in alleles if a == 1) / len(alleles)
    alt_allele = v.ALT[0] if v.ALT else "."
    out.write(f"{vid}\t{pos}\t{alt_allele}\t{alt_freq:.4f}\n")
    n_written += 1

out.close()
print(f"✓ Frequency table → {os.environ['FREQ_FILE']}")
print(f"✓ SNPs written: {n_written}")
PY

  # --- Converti gli output di SINGER in tskit (.trees) ----------------------
  convert_exec="${BASE_DIR}/Singer-Linux/convert_to_tskit"
  if [[ -x "$convert_exec" ]]; then
    # Input/output prefix: usiamo lo stesso prefix passato a singer_master
    input_prefix="${WORK_DIR}/${OUT_PREFIX}"
    output_prefix="${WORK_DIR}/${OUT_PREFIX}"

    # Intervallo campioni: 0..N_MCMC (inclusivo come richiesto)
    start_idx=0
    end_idx="${N_MCMC:-100}"

    info "Converting SINGER outputs to tskit: ${start_idx}..${end_idx}"
    conv_cmd=( "$convert_exec" -input "$input_prefix" -output "$output_prefix" \
               -start "$start_idx" -end "$end_idx" )
    echo "[COMMAND] ${conv_cmd[*]}"
    if ! "${conv_cmd[@]}"; then
      warn "convert_to_tskit failed – check logs. Continuing workflow."
    else
      echo "✓ convert_to_tskit completed → ${output_prefix}.trees (and related)"
    fi
  else
    warn "convert_to_tskit not found at $convert_exec – skipping tskit conversion."
  fi

fi

# --- clean‑up: remove intermediate files we no longer need -------------------
rm -f "${WORK_DIR}"/*.rate \
      "${WORK_DIR}"/*.pairwise.coal \
      "${WORK_DIR}"/*.pairwise.bin \
      "${WORK_DIR}"/*.annot \
      "${WORK_DIR}"/*.haps \
      "${WORK_DIR}"/*.sample \
      "${WORK_DIR}"/*_run1.anc \
      "${WORK_DIR}"/*_run1.mut
      
  mark_done "$WORK_DIR"
  echo
  info "✓ Phase 1 completed — outputs in ${WORK_DIR}"
  printf "\n"
}


# ---------- PHASE 2  (BranchLengths → RelateToCLUES → CLUES2 → merge) ----------
phase2() {

  PH="phase2"

  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║       🧬  PHASE 2 – (requires Phase-1 outputs)  🧬      ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"

  # --- Scelta metodo per Phase-2 (Relate vs SINGER) -------------------------
  echo "Which method do you want to continue with?"
  echo "  1) Relate"
  echo "  2) SINGER"
  read -rp "Enter option (1/2): " PHASE2_METHOD
  [[ "$PHASE2_METHOD" =~ ^[12]$ ]] || error "Invalid method selection for Phase-2"

  # -------------------------------------------------------------------------
  # 1. INPUTS COMUNI
  # -------------------------------------------------------------------------
  read -rp  "Choose the chromosome to analyze (e.g. 2, 17, X): " chr
  read -e -rp "Enter population prefix used in Phase-1 (e.g. Finnish): " out_prefix
  export OUT_PREFIX="$out_prefix"

  RUN_ID="${OUT_PREFIX}_chr${chr}"

  # Base differente per Relate/SINGER (coerente con Phase-1)
  if [[ "$PHASE2_METHOD" == "1" ]]; then
    WORK_BASE="${BASE_DIR}/output_CLUES2Companion-Relate"
  else
    WORK_BASE="${BASE_DIR}/output_CLUES2Companion-Singer"
  fi

  WORK_DIR="${WORK_BASE}/${PH}/${RUN_ID}"
  mkdir -p "$WORK_DIR"

  # where to find Phase-1 output (auto-detect coerente col metodo scelto)
  phase1_dir_default="${WORK_BASE}/phase1/${out_prefix}_chr${chr}"
  read -e -rp "Phase-1 auto-detect directory [ENTER = ${phase1_dir_default}] or provide a different folder with phase1 outputs: " phase1_dir
  phase1_dir=${phase1_dir:-$phase1_dir_default}
  [[ -d "$phase1_dir" ]] || error "Phase-1 directory not found: $phase1_dir"

  # -------------------------------------------------------------------------
  # 2. files comuni: freq e snps (servono a entrambi i rami)
  # -------------------------------------------------------------------------
  freq_file=$(ls "${phase1_dir}/${OUT_PREFIX}_Frequency_chr${chr}_"*.txt 2>/dev/null | head -1)
  [[ -f "$freq_file" ]] || error "Frequency file not found in ${phase1_dir}"
  echo; info "Using frequency file: $freq_file"

  snp_file=$(ls "${phase1_dir}/${OUT_PREFIX}_SNPs_chr${chr}_"*.txt 2>/dev/null | head -1)
  [[ -f "$snp_file" ]] || error "SNP coordinate file not found in ${phase1_dir}"
  info "Using SNPs file: $snp_file"
  echo

  # -------------------------------------------------------------------------
  # 3. script/tool paths
  # -------------------------------------------------------------------------
  relate_root="$(dirname "$0")/Relate-Linux"
  clues_dir="$(dirname "$0")/CLUES2"
  rtc_py="${clues_dir}/RelateToCLUES.py"
  inf_py="${clues_dir}/inference.py"
  stc_py="${clues_dir}/SingerToCLUES.py"

  # Verifiche minime
  [[ -f "$inf_py" ]] || error "inference.py not found in $clues_dir"
  if [[ "$PHASE2_METHOD" == "1" ]]; then
    [[ -x "${relate_root}/bin/Relate" ]] || error "Relate binaries not found"
    [[ -f "$rtc_py" ]] || error "RelateToCLUES.py not found in $clues_dir"
  else
    [[ -f "$stc_py" ]] || error "SingerToCLUES.py not found in $clues_dir"
  fi

  # -------------------------------------------------------------------------
  # 4. CLUES2 parameters (comuni)
  # -------------------------------------------------------------------------
  read -rp "tCutoff (e.g. 1000): " tcutoff
  read -rp "df (e.g. 600): "      df_score
  read -rp "AncientSamps file (optional, press ENTER to skip): " anc_samps
  read -rp "AncientHaps  file (optional, press ENTER to skip): " anc_haps
  read -rp "Disable allele trajectory? (y/N): " no_traj
  read -rp "Confidence interval (e.g. 0.95) [press ENTER to skip]: " ci_val
  read -rp "TimeBins (a list of epoch breakpoints; optional, press ENTER to skip): " time_bins
  read -rp "Dominance coefficient (default: 0.5, additive model, press ENTER to skip): " h_val

  time_bins_args=()
  if [[ -n $time_bins ]]; then
    IFS=', ' read -r -a tb_array <<< "$time_bins"
    tb_array=( "${tb_array[@]}" )
    [[ ${#tb_array[@]} -gt 0 ]] && time_bins_args=( --timeBins "${tb_array[@]}" )
  fi

  #############################################################################
  # 5. Directory output Phase-2
  #############################################################################
  TREE_DIR="${WORK_DIR}/${OUT_PREFIX}_trees_chr${chr}"
  TIMES_DIR="${WORK_DIR}/${OUT_PREFIX}_times_chr${chr}"
  INFER_DIR="${WORK_DIR}/${OUT_PREFIX}_inference_chr${chr}"
  mkdir -p "$TREE_DIR" "$TIMES_DIR" "$INFER_DIR"

  # ------------ spinner + runner (invariati) --------------------------------
  spinner_bar() {
    local pid=$1 delay=${2:-1} max=${3:-30} count=0
    echo "THIS MAY TAKE A WHILE, PLEASE WAIT."
    printf "["
    while kill -0 "$pid" 2>/dev/null; do
      (( ++count > max )) && { printf "\r[\033[K"; count=1; }
      printf "▊"; sleep "$delay"
    done
    printf "] done!\n"
  }
  run_step() {
    local func=$1 delay=${2:-0.5} max=${3:-30}
    local tmpwarn; tmpwarn=$(mktemp)
    (
      warn() { echo "[WARN] $*" >>"$tmpwarn"; }
      "$func"
    ) &
    local step_pid=$!
    spinner_bar "$step_pid" "$delay" "$max"
    wait "$step_pid" || true
    [[ -s $tmpwarn ]] && cat "$tmpwarn"
    rm -f "$tmpwarn"
  }

  if [[ "$PHASE2_METHOD" == "1" ]]; then
    #############################################################################
    # ==========================  RELATE BRANCH  ================================
    #############################################################################
    read -rp "Importance sampling of branch lengths: " num_samples
    [[ "$num_samples" =~ ^[0-9]+$ ]] || error "must be a positive integer"
    
    # 3R. Relate artifacts da Phase-1
    gs_prefix="${phase1_dir}/${OUT_PREFIX}_GS+COAL_chr${chr}"
    [[ -f "${gs_prefix}.anc" || -f "${gs_prefix}.anc.gz" ]] \
          || error "*.anc file not found for prefix $gs_prefix"

    coal_file="${phase1_dir}/${OUT_PREFIX}_EPS4COAL_chr${chr}.coal"
    [[ -f "$coal_file" ]] || error "*.coal file not found: $coal_file"

    # 4R. Derived files glob (Phase-1 Relate)
    derived_glob="${phase1_dir}/${OUT_PREFIX}_Derived_*.txt"
    compgen -G "$derived_glob" >/dev/null || error "Derived files not found in ${phase1_dir}"
    

    #############################################################################
    # --- Step A: SampleBranchLengths (Relate) ----------------------------------
    #############################################################################
    stepA() {
      grep -v -E '^(rsID|rsid|#)' "$snp_file" |
      while read -r rsid pos; do
        [[ -z $rsid || -z $pos ]] && continue
        rand_seed() { echo $(( ( RANDOM << 15 ) | RANDOM )); }
        "${relate_root}/scripts/SampleBranchLengths/SampleBranchLengths.sh" \
            --input "$gs_prefix" \
            --output "${TREE_DIR}/${rsid}" \
            --first_bp "$pos" --last_bp "$pos" \
            --format n --num_samples "$num_samples" \
            --coal "$coal_file" -m 1.25e-8 --seed "$(rand_seed)" &>/dev/null \
            || warn "BranchLengths failed for $rsid"
      done
    }
    echo; info "Step A: SampleBranchLengths"
    stepA
    echo "✓ Step A completed"

    #############################################################################
    # --- Step B: RelateToCLUES -------------------------------------------------
    #############################################################################
    stepB() {
      for dfile in $derived_glob; do
        rsid=$(basename "$dfile" | sed "s/^${OUT_PREFIX}_Derived_//; s/\.txt$//")
        nwx="${TREE_DIR}/${rsid}.newick"
        if [[ -f $nwx ]]; then
          python "$rtc_py" --DerivedFile "$dfile" --RelateSamples "$nwx" \
                           --out "${TIMES_DIR}/${rsid}" &>/dev/null \
                           || warn "RTC failed for $rsid"
        else
          warn "newick not found for $rsid"
        fi
      done
    }
    echo; info "Step B: RelateToCLUES"
    stepB
    echo "✓ Step B completed"

  else
    #############################################################################
    # ==========================   SINGER BRANCH  ===============================
    #############################################################################

    # Regione di interesse
    read -rp "Enter START bp of region: " singer_start
    read -rp "Enter END bp of region: " singer_end
    [[ "$singer_start" =~ ^[0-9]+$ && "$singer_end" =~ ^[0-9]+$ ]] || error "Start/End must be integers"

    # Ne richiesto da inference.py in assenza di .coal
    read -rp "Effective population size (Ne) for inference.py: " singer_ne
    [[ "$singer_ne" =~ ^[0-9]+$ ]] || error "Ne must be a positive integer"

    # `.trees` creato in Phase-1 SINGER via convert_to_tskit
    trees_path="${phase1_dir}"
    [[ -d "$trees_path" ]] || error "TSKIT trees directory not found: $trees_path"

    # Estrai tutti gli SNP nell’intervallo [singer_start, singer_end]
    mapfile -t snp_list < <(awk -v s="$singer_start" -v e="$singer_end" '$2>=s && $2<=e {print $1"\t"$2}' "$snp_file")
    [[ ${#snp_list[@]} -gt 0 ]] || error "No SNPs found in interval $singer_start-$singer_end"

    info "Found ${#snp_list[@]} SNPs in interval $singer_start-$singer_end"

    # Prepara freq file ridotto con solo SNPs selezionati
    freq_one="${WORK_DIR}/__freq_selected.txt"
    head -1 "$freq_file" > "$freq_one"   # header
    for line in "${snp_list[@]}"; do
        rs=$(echo "$line" | cut -f1)
        grep -w "$rs" "$freq_file" >> "$freq_one"
    done
    freq_file="$freq_one"

    # Loop SingerToCLUES su ogni SNP dell’intervallo
    for line in "${snp_list[@]}"; do
        rs=$(echo "$line" | cut -f1)
        pos=$(echo "$line" | cut -f2)
        echo; info "SingerToCLUES for $rs at POS=$pos"
        python "$stc_py" \
          --position "$pos" \
          --tree_path "$trees_path" \
          --output "${TIMES_DIR}/${rs}" \
          || warn "SingerToCLUES failed for $rs ($pos)"
    done

    # Saltiamo Step A/B (Relate). Il resto (Step C + Merge + Plot) resta identico.
    # In stepC ricordati di usare:
    #   cmd+=( --Ne "$singer_ne" )
    # invece di --coal
  fi

  #############################################################################
  # --- Step C – inference.py (comune a entrambi i rami) ----------------------
  #############################################################################
  stepC() {
    grep -v -E '^(rsID|rsid|#)' "$freq_file" |
    while read -r rs pos _ freq; do
      tfile="${TIMES_DIR}/${rs}_times.txt"
      if [[ -f $tfile ]]; then
        cmd=( python "$inf_py" --times "$tfile" --popFreq "$freq" \
              --out "${INFER_DIR}/${rs}" --tCutoff "$tcutoff" --df "$df_score" )
        # Per Relate abbiamo $coal_file; in SINGER di solito non serve/è assente.
        if [[ "$PHASE2_METHOD" == "1" ]]; then
          # Relate → usa coal se esiste
          if [[ -n ${coal_file:-} && -f ${coal_file:-/dev/null} ]]; then
            cmd+=( --coal "$coal_file" )
          fi
        else
          # Singer → obbligatorio specificare Ne
          cmd+=( --N "$singer_ne" )
        fi

        [[ -n ${ci_val:-}       ]] && cmd+=( --CI "$ci_val" )
        [[ -n ${h_val:-}       ]] && cmd+=( --h "$h_val" )
        [[ -f ${anc_samps:-}    ]] && cmd+=( --ancientSamps "$anc_samps" )
        [[ -n ${time_bins_args+x} ]] && cmd+=( "${time_bins_args[@]}" )
        [[ -f ${anc_haps:-}     ]] && cmd+=( --ancientHaps "$anc_haps" )
        [[ $no_traj =~ ^[yY]$   ]] && cmd+=( --noAlleleTraj )
        "${cmd[@]}" || warn "inference failed for $rs"
      else
        warn "times not found for $rs"
      fi
    done
  }

  echo; info "Step C: inference.py"
  stepC
  echo "✓ Step C completed"

  # -------------------------------------------------------------------------
  # Step D: merge inference outputs (uno o molti rs; funziona per entrambi)
  # -------------------------------------------------------------------------
  echo
  info "Step D: merge inference outputs"
  merged="${WORK_DIR}/${OUT_PREFIX}_merged_inference_chr${chr}.tsv"

  # inizializza una tabella vuota
  printf "rsID\tPOS\tder_freq" >"$merged"

  # trova un primo inference per ricavare le colonne (gestisce anche il caso singolo)
  first_inf=$(ls "${INFER_DIR}"/*_inference.txt 2>/dev/null | head -1) || true
  [[ -f "$first_inf" ]] || error "No inference files found in ${INFER_DIR}"

  cols=$(head -1 "$first_inf")
  printf "\tlogLR\t-log10(p)" >>"$merged"

  n_sel=$(echo "$cols" | tr '\t' '\n' | grep -c '^SelectionMLE')
  for k in $(seq "$n_sel"); do
      printf "\tSelectionMLE${k}\tCI${k}_lower\tCI${k}_upper" >>"$merged"
  done
  printf "\n" >>"$merged"

  grep -v -E '^(rsID|rsid|#)' "$freq_file" |
  while read -r rs pos _ freq; do
      inf_file="${INFER_DIR}/${rs}_inference.txt"
      [[ -f $inf_file ]] || { warn "inference missing for $rs"; continue; }

      IFS=$'\t' read -r -a hdr  <<<"$(head -1  "$inf_file")"
      IFS=$'\t' read -r -a data <<<"$(sed -n '2p' "$inf_file")"

      col_idx() {
        local needle=$1
        for i in "${!hdr[@]}"; do
          [[ ${hdr[i]} == *"$needle"* ]] && { echo "$i"; return; }
        done
        echo -1
      }
      loglr_idx=$(col_idx "logLR")
      logp_idx=$(col_idx "-log10(")
      if [[ $loglr_idx -lt 0 || $logp_idx -lt 0 ]]; then
        warn "Header mismatch in $inf_file"; continue
      fi

      out=("${rs}" "${pos}" "${freq}" "${data[loglr_idx]}" "${data[logp_idx]}")

      ci_lo=()  ci_up=()
      ci_file="${INFER_DIR}/${rs}_CI.txt"
      if [[ -f $ci_file ]]; then
        while IFS=$'\t' read -r epoch lo up; do
            ci_lo+=("$lo") ; ci_up+=("$up")
        done < <(tail -n +2 "$ci_file")
      fi

      for k in $(seq "$n_sel"); do
        sel_idx=$(col_idx "SelectionMLE${k}")
        sel_val=${data[sel_idx]:-"NA"}
        if [[ -f $ci_file ]]; then
            lo=${ci_lo[k-1]:-NA}
            up=${ci_up[k-1]:-NA}
        else
            lo="NA"; up="NA"
        fi
        out+=("$sel_val" "$lo" "$up")
      done

      printf "%s\n" "$(IFS=$'\t'; echo "${out[*]}")" >>"$merged"
  done

  info "Merged file → $merged"

  mark_done "$WORK_DIR"
  echo
  info "Generating integrated plot"
  export MERGED_TSV="$merged"
  export PLOT_PREFIX="${WORK_DIR}/${OUT_PREFIX}_clues"

  python - <<'PYCODE'
import os, pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from adjustText import adjust_text
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import offset_copy

merged_tsv = os.environ["MERGED_TSV"]
prefix     = os.environ["PLOT_PREFIX"]

df = pd.read_csv(merged_tsv, sep="\t")
df["p_raw"] = 10**(-df["-log10(p)"])
star = lambda p: "***" if p<1e-3 else ("**" if p<1e-2 else ("*" if p<5e-2 else ""))

epochs = [c for c in df.columns if c.startswith("SelectionMLE")]
n_ep   = len(epochs)

for k in range(1, n_ep+1):
    sel   = df[f"SelectionMLE{k}"]
    loCI  = df[f"CI{k}_lower"]
    upCI  = df[f"CI{k}_upper"]
    yerr  = np.vstack([np.abs(sel-loCI), np.abs(upCI-sel)])

    fig, ax = plt.subplots(figsize=(11,4))
    sc = ax.scatter(df["POS"], sel, c=-np.log10(df["p_raw"]),
                    cmap="Reds", vmin=0, vmax=10, s=40, zorder=3)
    ax.errorbar(df["POS"], sel, yerr=yerr,
                fmt='none', ecolor="grey", elinewidth=1,
                capsize=1.5, zorder=2)

    texts = []
    for i, r in df.iterrows():
        sig = star(r["p_raw"])
        if sig:
            trans_star = offset_copy(ax.transData, fig, y=4, units='points')
            ax.text(r["POS"], upCI[i], sig, ha='center', va='bottom',
                    transform=trans_star, color="firebrick", fontsize=7)
            trans_rs = offset_copy(ax.transData, fig, y=12, units='points')
            texts.append(ax.text(r["POS"], upCI[i], r["rsID"],
                                 ha='center', va='bottom',
                                 transform=trans_rs, fontsize=7))

    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylabel("Selection coefficient (s)")
    ax.set_title(f"Epoch {k}   (* P<0.05  ** P<0.01  *** P<0.001)")
    ax.grid(alpha=.4)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_: f"{int(x)}"))

    cbar = fig.colorbar(sc, ax=ax, pad=0.01)
    cbar.set_label("−log10 P  (0–10)")

    fn = f"{prefix}_epoch{k}.png"
    fig.tight_layout()
    fig.savefig(fn, dpi=300)
    fig.savefig(fn.replace(".png", ".pdf"))
    plt.close(fig)
    print("Wrote:", Path(fn).name)
PYCODE

  printf "\n"
  info "Phase-2 completed."
  echo
}

# ---------- PHASE-3  (re-inference + sweep-dating for ONE SNP) --------------
phase3() {
  PH="phase3"
  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║        ⏱️  PHASE 3 – Dating a selective sweep  ⏱️         ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"
  # ------------------------ user input --------------------------------------
  read -rp "Choose chromosome to analyze (e.g. 2, 17, X): "  chr
  read -rp "Same prefix population name as used in Phase 2 (e.g. Finnish): "    pop
  read -rp "rsID of SNP to date   (e.g. rs123): "  rsid
  read -rp "df score for CLUES2 (e.g. 600): "      df_score
  read -rp "Initial epoch to scan for initial onset (e.g. 0 or 50): "   G_START
  read -rp "Final epoch to scan for initial onset (e.g. 500 or 1000): " G_END
  read -rp "Non overlapping windows size  [default is 50]: "       STEP ; STEP=${STEP:-50}
  # ------------------------ paths produced in phase-1/2 ---------------------
  P1_DIR="${WORK_BASE}/phase1/${pop}_chr${chr}"
  P2_DIR="${WORK_BASE}/phase2/${pop}_chr${chr}"
  [[ -d "$P1_DIR" && -d "$P2_DIR" ]] || \
        error "phase1/phase2 output not found for ${pop}_chr${chr}"
  
  # Phase-3 working directory
  RUN_ID="${pop}_chr${chr}"
  WORK_DIR="${WORK_BASE}/phase3/${RUN_ID}"
  mkdir -p "$WORK_DIR"
  printf "\n"
  info "[INFO]Phase 3 working dir: $WORK_DIR"

  COAL="${P1_DIR}/${pop}_EPS4COAL_chr${chr}.coal"
  [[ -f "$COAL" ]] || error "*.coal file missing"

  DERIVED="${P1_DIR}/${pop}_Derived_${rsid}.txt"
  [[ -f "$DERIVED" ]] || error "Derived file missing ($DERIVED)"

  # ----- derived-allele frequency (needed by inference.py) ------------------
  FREQ_FILE=$(ls "${P1_DIR}/${pop}_Frequency_chr${chr}_"*.txt | head -1)
  [[ -f "$FREQ_FILE" ]] || error "frequency file not found"
  freq=$(awk -v RSID="$rsid" '$1==RSID{print $4; exit}' "$FREQ_FILE")
  [[ -n "$freq" ]] || error "derived allele frequency for $rsid not found"

  # ------------------------ output directories ------------------------------
  DAT_DIR="${WORK_DIR}/Dating";              mkdir -p "$DAT_DIR"
  INF_OUT="${DAT_DIR}/${rsid}"             # prefix for first run
  OUT_JSON="${DAT_DIR}/${rsid}_Dating.json"

  RELATE_ROOT="$(dirname "$0")/Relate-Linux"
  CLUES_DIR="$(dirname "$0")/CLUES2"
  inf_py="${CLUES_DIR}/inference.py"
  rtc_py="${CLUES_DIR}/RelateToCLUES.py"

# -------------------------------------------------------------------------
# (A)  Sliding-window scan: one CLUES2 run per interval
# -------------------------------------------------------------------------
TIMES0="${P2_DIR}/${pop}_times_chr${chr}/${rsid}_times.txt"
[[ -f "$TIMES0" ]] || error "times file missing ($TIMES0)"

# Copiamo il file–times nella cartella di Phase-3 (così rimane tutto lì)
TIMES_FILE="${WORK_DIR}/${rsid}_times.txt"
cp "$TIMES0" "$TIMES_FILE"

# ---------------- build BREAKS from user input -----------------
if (( STEP <= 0 )) || (( G_START >= G_END )); then
    error "STEP must be >0 and G_START < G_END"
fi

# array BREAKS: G_START, G_START+STEP, …, G_END
if declare -F mapfile >/dev/null; then      # Bash ≥ 4
    mapfile -t BREAKS < <( seq "$G_START" "$STEP" "$G_END" )
else                                        # fallback Bash 3.x
    IFS=$'\n' read -r -d '' -a BREAKS < <( seq "$G_START" "$STEP" "$G_END"; printf '\0' )
fi

# assicura che l'ultima boundary sia ESATTAMENTE G_END
last_break=${BREAKS[${#BREAKS[@]}-1]}
(( last_break != G_END )) && BREAKS+=( "$G_END" )

# ---------------- stampa riepilogo BREAKS -----------------------------------
printf "\nBreak-points (STEP (%s-%s)):\n" "left" "right"
for (( i=0; i<${#BREAKS[@]}-1; i++ )); do
    left=${BREAKS[i]}
    right=${BREAKS[i+1]}
    printf "%-4d (%d-%d)\n"  "$STEP"  "$left" "$right"
done
printf "tCutoff      : %d\n"  "${BREAKS[${#BREAKS[@]}-1]}"
printf "Window size  : %d\n\n"  "$STEP"

# ---------------- loop sulle finestre ---------------------------------------
echo "▶  Scanning windows:"
for (( i=0; i<${#BREAKS[@]}-1; i++ )); do
    left=${BREAKS[i]}
    right=${BREAKS[i+1]}

    win_pref="${INF_OUT}_${left}_${right}"
    info "    CLUES2 • window [${left}-${right}]"

    cmd=( python "$inf_py"
          --times   "$TIMES_FILE"
          --popFreq "$freq"
          --out     "$win_pref"
          --df      "$df_score"
          --CI      0.95
          --tCutoff "$right"
          --coal    "$COAL"
          --noAlleleTraj )

    (( left > 0 )) && cmd+=( --timeBins "$left" )

    "${cmd[@]}" >/dev/null 2>&1 || {
        warn "      ➜ CLUES2 failed – window skipped"
        continue
    }
done

# -------------------------------------------------------------------------
#  (A-2)  Collect inference results and determine onset
# -------------------------------------------------------------------------
# Expects the following environment variables to be set by the shell layer:
#   INF_OUT – root prefix to inference files (e.g.  .../rs4988235)
#   rsid    – SNP identifier
#   pop     – population label
#   chr     – chromosome
# -------------------------------------------------------------------------
export INF_OUT rsid pop chr # let's read the variable for PY block
python - <<'PY'
import glob, os, re, pandas as pd, json, pathlib, numpy as np, sys

# ------------------------- gather files ----------------------------------
pref_root = os.environ['INF_OUT']                 # e.g.  /path/to/rs4988235
def win_start(path: str) -> int:
    """
    Estrae la generazione iniziale della finestra dal nome file
    (…_400_500_inference.txt → 400).
    """
    m = re.search(r'_(\d+)_(\d+)_inference\.txt$', os.path.basename(path))
    if not m:
        raise ValueError(f"Nome file non atteso: {path}")
    return int(m.group(1))

files = sorted(
    glob.glob(f"{pref_root}_*_inference.txt"),
    key=win_start,
    reverse=True) 

windows = []  # will hold tuples: (start, end, s_MLE, p_val)
for f in files:
    m = re.search(r'_(\d+)_(\d+)_inference\.txt$', f)
    if not m:
        continue
    start, end = map(int, m.groups())
    row = pd.read_csv(f, sep="\t").iloc[0]
    # take the *last* SelectionMLE column (highest k)
    last_k = max(int(c.split("SelectionMLE")[1])
                 for c in row.index if c.startswith("SelectionMLE"))
    s_val = row[f"SelectionMLE{last_k}"]
    p_val = 10**(-row["-log10(p-value)"]) if "-log10(p-value)" in row else np.nan
    windows.append((start, end, s_val, p_val))

if not windows:
    sys.exit("[ERROR] No inference files found – check INF_OUT prefix and path.")

# ----------------------- helper to find consecutive -----------------------

def find_consecutive(pos_windows, k):
    """Return (window, k) if k consecutive windows have s>0, else None."""
    for i in range(len(pos_windows) - k + 1):
        segment = pos_windows[i:i + k]      # ancient → recent order
        if all(w[2] > 0 for w in segment):
            return segment[0], k            # left‑most (oldest) window and k
    return None

# ----------------------- onset detection logic ---------------------------
criterion_messages = []
result = None
for k in (4, 3, 2):
    res = find_consecutive(windows, k)
    if res:
        result = res
        method_used = f"{k} consecutive windows > 0"
        if k == 4:
            print("Detected onset using criterion: 4 consecutive windows > 0\n")
        elif k == 3:
            print("Criterion 4 consecutive s > 0 not satisfied; using 3 consecutive windows > 0\n")
        else:  # k == 2
            print("[WARNING] Criterion 4 consecutive s > 0 not satisfied and 3 consecutive s > 0 not satisfied; using 2 consecutive windows > 0\n")
        break

# If none of the consecutive‑positive criteria hit, fall back to previous heuristics
if result is None:
    print("[WARNING] No block of 2 consecutive positive‑s windows found – falling back to maximum s estimate\n")
    result = (next((w for w in windows if w[2] > 0), None) or max(windows, key=lambda t: t[2]), 1)
    method_used = "fallback: first positive or max s"

(onset_win, k_used) = result
st, en, s_mle, _ = onset_win
median_gen   = int(round((st + en) / 2))
median_years = median_gen * 28

# ----------------------- screen output -----------------------------------
print(f"\nInitial onset (median of window) ≈ {median_gen} generations  (≈ {median_years} years)")
print(f"   window   : {st} – {en} generations")
print(f"   s(MLE)   : {s_mle:.5f}")
print(f"   method   : {method_used}\n")

# ----------------------- JSON output -------------------------------------
json_out = pathlib.Path(f"{pref_root}_InitialOnset_Dating.json")
json_out.write_text(json.dumps(dict(
    rsID                 = os.environ["rsid"],
    population           = os.environ["pop"],
    chromosome           = os.environ["chr"],
    median_onset_gen     = median_gen,
    median_onset_years   = median_years,
    epoch_start          = st,
    epoch_end            = en,
    epoch_start_years    = st * 28,
    epoch_end_years      = en * 28,
    method               = method_used,
    s_MLE                = round(s_mle, 5)
), indent=2))

print(f"[INFO] JSON written to: {json_out}\n")
PY

# -------------------------------------------------------------------------
#  (A-2)  initial onset done – JSON for first dating saved
# -------------------------------------------------------------------------

# >>> ask if the user want proceed with the bootstraps, otherwise exit <<<
read -rp "Proceed with bootstrap dating? [y/N]: " REPLY
if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
    echo -e "\nBootstrap step skipped – Phase 3 completed."
    return            # esce dalla funzione phase3 (o 'exit 0' se non è in funzione)
fi

# -----------------------------------------------------------------
#  (B)  bootstrap parameters
# -----------------------------------------------------------------
echo -e "\n▶  Bootstrap settings"
read -rp "Start point to scan in generations ago before initial onset (e.g. 0 or 100): "   G_START
read -rp "End point to scan in generations ago after initial onset (e.g. 500): " G_END
read -rp "Non-overlapping time bin size  [default is 25]: "       STEP ; STEP=${STEP:-25}
read -rp "Number of bootstrap replicates  [default is 100]: "       NBOOT; NBOOT=${NBOOT:-100}
read -rp "df score for CLUES2  [default is 450]: "      DF   ; DF=${DF:-450}
read -rp "Importance sampling of branch lengths: " num_samples

# ---------- build BREAKS ----------------------------------------------------
if declare -F mapfile >/dev/null; then          # Bash ≥ 4
    mapfile -t BREAKS < <( seq "$G_START" "$STEP" "$G_END" )
else                                            # Bash 3.x fallback
    IFS=$'\n' read -r -d '' -a BREAKS < <( seq "$G_START" "$STEP" "$G_END"; printf '\0' )
fi

# assicura che l’ultima boundary sia esattamente G_END
last_break=${BREAKS[${#BREAKS[@]}-1]}
(( last_break != G_END )) && BREAKS+=( "$G_END" )

# sanity-check
(( ${#BREAKS[@]} == 0 )) && error "seq produced no break-points"

TCUT=$(( G_END + STEP ))        # --tCutoff passato a CLUES

# ---------------- stampa riepilogo ------------------------------------------
printf "\nBreak-points (STEP (%s-%s)):\n" "left" "right"
for (( i=0; i<${#BREAKS[@]}-1; i++ )); do
    left=${BREAKS[i]}
    right=${BREAKS[i+1]}
    printf "%-4d (%d-%d)\n"  "$STEP"  "$left" "$right"
done
printf "tCutoff      : %d\n"   "$TCUT"
printf "Replicates   : %d\n\n" "$NBOOT"

##############################################################################
#  C)  N BOOTSTRAP REPLICATES                           #
##############################################################################
BOOT_DIR="${DAT_DIR}/bootstrap_${rsid}"; mkdir -p "$BOOT_DIR"
GS_PREFIX="${P1_DIR}/${pop}_GS+COAL_chr${chr}"
POS=$(grep -m1 -w "$rsid" "${P1_DIR}/${pop}_SNPs_chr${chr}_"*.txt | cut -f2) \
     || error "bp position for $rsid not found"

echo -e "\n▶  Generating $NBOOT bootstrap trees & CLUES2 runs"

for (( rep=1; rep<=NBOOT; rep++ )); do
    seed=$(( ( RANDOM << 15 ) ^ RANDOM ))
    info "[bootstrap $rep]  SampleBranchLengths (seed=$seed)"
    nw_pref="${BOOT_DIR}/nw_${rep}"

    "${RELATE_ROOT}/scripts/SampleBranchLengths/SampleBranchLengths.sh" \
        --input "$GS_PREFIX" --output "$nw_pref" \
        --first_bp "$POS" --last_bp "$POS" \
        --format n --num_samples "$num_samples" --coal "$COAL" \
        --seed "$seed" -m 1.25e-8 &>/dev/null

    info "  RelateToCLUES"
    python "$rtc_py" --DerivedFile "$DERIVED" \
                     --RelateSamples "${nw_pref}.newick" \
                     --out "$nw_pref" &>/dev/null

# ---- CLUES2 runs each windows --------------------------------------
for (( i=0; i<${#BREAKS[@]}-1; i++ )); do
    left=${BREAKS[i]} ; right=${BREAKS[i+1]}
    win_pref="${nw_pref}_${left}_${right}"
    #log="${win_pref}_clues.log"
    
    info "    CLUES2 inference • window [$left - $right]"

    cmd=( python "$inf_py" --times "${nw_pref}_times.txt"
          --popFreq "$freq" --out "$win_pref"
          --df "$DF" --CI 0.95 --tCutoff "$right"
          --coal "$COAL" --noAlleleTraj )
    (( left > 0 )) && cmd+=( --timeBins "$left" )

    "${cmd[@]}" >/dev/null 2>&1 || {
        warn "      ➜ CLUES2 failed – window skipped"
        continue
    }
done

# ---------- AFTER all windows: concat Epoch-2 -------------------------
combo="${BOOT_DIR}/bootstrap_${rep}_${rsid}.txt"
: > "$combo"                            # svuota/crea

# header fisso
printf "logLR\t-log10(p-value)\tEpoch2_start\tEpoch2_end\tSelectionMLE2\n" \
       >> "$combo"

for f in "${nw_pref}"_*"_inference.txt"; do
    [[ -f $f ]] || continue

    awk -F'\t' '
        NR==1 {
            for(i=1;i<=NF;i++){
                if($i=="logLR")            c_l=i;
                if($i=="-log10(p-value)")  c_p=i;
                if($i=="Epoch2_start")     c_s=i;
                if($i=="Epoch2_end")       c_e=i;
                if($i=="SelectionMLE2")    c_m=i;
            }
            next
        }
        # stampa solo se tutte le colonne esistono nella riga
        (c_l && c_p && c_s && c_e && c_m && NF>=c_m) {
            printf("%s\t%s\t%s\t%s\t%s\n",
                   $c_l, $c_p, $c_s, $c_e, $c_m);
        }
    ' "$f" >> "$combo"
done

rm "${nw_pref}"_*"_inference.txt"



done   # -------------------------end bootstrap cycle----------------------


##############################################################################
#  D)  CALCULATE ONSET FOR ALL BOOTSTRAP                                    #
##############################################################################
output="onset_bootstraps.txt"
> "$output"   # svuota o crea il file

echo -e "\n▶  Computing onset for each bootstrap…"
onset_vals=()
method_vals=()

for (( rep=1; rep<=NBOOT; rep++ )); do
    combo="${BOOT_DIR}/bootstrap_${rep}_${rsid}.txt"
    if [[ ! -s $combo ]]; then
        warn "Bootstrap $rep – file $combo empty or missing, skip"
        continue
    fi

onset_and_method=$(python3 - <<'PY' "$combo"
import sys, pandas as pd, math

wins = []
df = pd.read_csv(sys.argv[1], sep="\t")

for _, row in df.iterrows():
    st, en = int(row["Epoch2_start"]), int(row["Epoch2_end"])
    s      = float(row["SelectionMLE2"])
    wins.append((st, en, s))

wins.sort(key=lambda w: w[0], reverse=True)

def median(st, en):
    return int(round((st + en) / 2))

def find_onset(wlist, k):
    for i in range(len(wlist) - k + 1):
        seg = wlist[i:i + k]
        if all(w[2] > 0 for w in seg):
            st, en, _ = seg[0]
            return median(st, en), f"{k} consecutive s>0"
    return None, None

onset, method = None, None
for k in (4, 3, 2):
    onset, method = find_onset(wins, k)
    if onset is not None:
        break

if onset is None:
    for st, en, s in wins:
        if s > 0:
            onset  = median(st, en)
            method = "first positive s"
            break

if onset is None:
    st, en, _ = max(wins, key=lambda w: w[2])
    onset  = median(st, en)
    method = "max s"

print(f"{onset}\t{method}")
PY
)

IFS=$'\t' read -r onset method <<< "$onset_and_method"
method_vals+=("$method")
line=$(printf "Bootstrap %2d : %s gen  (%s)" "$rep" "$onset" "$method")

    echo  "$line"        # stampa a terminale
    echo  "$line" >> "$output"  # salva in append

    [[ $onset =~ ^[0-9]+$ ]] && onset_vals+=("$onset")
done



##############################################################################
#  E)  STATISTICHE                                                           #
##############################################################################
echo -e "\n▶  Summary statistics"
echo "Valid onsets collected: ${#onset_vals[@]}"

if (( ${#onset_vals[@]} < 2 )); then
    error "Fewer than 2 valid onsets – CI cannot be computed"
fi

stats=$(python3 - <<'PY' "${onset_vals[@]}"
import sys, numpy as np, math
vals = np.array(list(map(int, sys.argv[1:])), dtype=int)
print(int(np.median(vals)),
      int(math.floor(np.percentile(vals, 2.5))),
      int(math.ceil(np.percentile(vals, 97.5))))
PY
)

# --- split robusto, forziamo IFS a “spazio” / “newline” --------------------
IFS=$' \t\n' read -r MEDIAN CI_LOW CI_HIGH <<< "$stats"

echo
#echo "Median onset : $MEDIAN gen (~$(( MEDIAN * 28 )) yr)"
echo "95 % CI      : $CI_LOW – $CI_HIGH gen "\
     "(~$(( CI_LOW * 28 )) – $(( CI_HIGH * 28 )) yr)"
     

# ------------------------------------------------------------------
#  Metodo prevalente tra i bootstrap
# ------------------------------------------------------------------
if ((${#method_vals[@]})); then
    prevailing_method=$(printf '%s\n' "${method_vals[@]}" \
                        | sort | uniq -c | sort -nr | head -1 \
                        | awk '{print $2" ("$1"×)"}')
else
    prevailing_method="NA"
fi


##############################################################################
#  F)  JSON OUTPUT                                                           #
##############################################################################
export PREV_METHOD="$prevailing_method"

python3 - <<PY
import json, pathlib, os
out = dict(
    rsID              = os.environ["rsid"],
    population        = os.environ["pop"],
    chromosome        = os.environ["chr"],
    #onset_median_gen  = int("$MEDIAN"),
    #onset_median_year = int("$MEDIAN")*28,
    CI95_low_gen      = int("$CI_LOW"),
    CI95_high_gen     = int("$CI_HIGH"),
    CI95_low_year     = int("$CI_LOW")*28,
    CI95_high_year    = int("$CI_HIGH")*28,
    bootstraps        = int("$NBOOT"),
    step              = int("$STEP"),
    method            = os.environ["PREV_METHOD"]
)

out_dir = pathlib.Path("$BOOT_DIR").parent          # ← usa la variabile bash
out_path = out_dir / f"{os.environ['rsid']}_Boostraps_onset_Dating+CI.json"
out_path.write_text(json.dumps(out, indent=2))

print(f"JSON saved → {out_path}")
PY

}
# ---------- dispatch ----------------------------------------------------------
case "$OPTION" in
  1) phase1 ;;
  2) phase2 ;;
  3) phase3 ;;
  *) error "Invalid option" ;;
esac
