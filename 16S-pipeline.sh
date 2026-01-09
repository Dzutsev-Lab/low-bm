#!/usr/bin/env bash
# =============================================================================
# RNA 16S (nested-like PCR; Nextera UMIs on R2) — Refactored, Safe, Reproducible
# =============================================================================
# 
# - UMI parsing ( extracts 16 nt immediately upstream of primer GGACTAC).
#
# 
# 
#
#
# Inputs (Biowulf/Helix paths)
#   RAW_DIR: directory with *_R1_001.fastq[.gz] and *_R2_001.fastq[.gz]
#   OUT_ROOT: working/output folder that will contain original.fastq/ and clean.fastq/
#   HUMAN_REF: bwa index prefix for human/cDNA (clean.human.fa)
#   VIRAL_REF: bwa index prefix for viral fasta (idexes/all.viral.fna)
#   BACT16S_REF: bwa index prefix for bacterial rRNA set (all.rrna.bacteria)
#   MOTHUR_TEMPLATE / MOTHUR_TAX: mothur ncbi20 references
#
# What you get
#   OUT_ROOT/
#     original.fastq/      (original *.fastq)
#     clean.fastq/         (final FASTA, UMI maps, per-sample files)
#     taxonomy/            (mothur outputs + count tables)
#     logs/                (per-step logs)
#     xin1.stats.tsv       (ID, Total_reads, Human_mapped_reads, Clean_reads)
#
# Usage
#   1) Edit CONFIG section below.
#   2) module load seqkit seqtk bwa samtools trimmomatic mothur pigz
#   3) bash xin16S_refactored_pipeline.sh 2> logs/run.stderr.log | tee logs/run.stdout.log
#
# Notes
#   - Assumes UMI is 16 nt immediately before the primer motif GGACTAC on R2.
#   - All mapping steps are single-end and use FASTA to mirror your original.
#   - For large runs, wrap the per-sample loops with swarm if desired (hooks provided).
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

module load trimmomatic seqkit seqtk bwa samtools mothur pigz java

# ------------------------------
# CONFIG — EDIT ME
# ------------------------------
RAW_DIR="./Select_Trial_Data"               # input FASTQs
OUT_ROOT="./Select_Trial_Output"          # output root
PRIMER_MOTIF="GGACTAC"                                          # 16nt UMI is immediately upstream of this
UMI_LEN="16"

HUMAN_REF="./Ref_Data/Mus_musculus.GRCm38.cdna.all.fa" 
VIRAL_REF="./Ref_Data/all.viral.fna"
BACT16S_REF="./Ref_Data/all.rrna.bacteria"

# mothur references
MOTHUR_TEMPLATE="./Ref_Data/ncbi20.fasta"
MOTHUR_TAX="./Ref_Data/ncbi20.tax"

# Trimming (match your historical settings)
TRIMMOMATIC_JAR="$TRIMMOJAR"                                   # set by module trimmomatic
TRIM_ARGS=("AVGQUAL:30")                                        # keep semantics; optionally add ILLUMINACLIP/MINLEN/CROP

# Negative control filtering (optional)
NEG_DB_TOP_K="10000"                                            # top N negative-control uniques to include
NEG_AS_KEEP_LT="160"                                            # keep if AS < this (i.e., not similar to neg DB)

THREADS="2"

# Swarm toggle: if set to 1, writes swarm files alongside loops (does not submit)
EMIT_SWARM="0"

# ------------------------------
# Helpers
# ------------------------------

# LOG: prints a formatted string including date/time and status update on current processes
# Input: some record to be logged 
#
log() { printf "[log][%s] %s\n" "$(date +"%F %T")" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

# Check tools early
for tool in seqkit seqtk bwa samtools mothur pigz java; do need "$tool"; done
[[ -f "$TRIMMOMATIC_JAR" ]] || die "Trimmomatic jar not found via \$TRIMMOJAR ($TRIMMOMATIC_JAR). Load module trimmomatic."

# Prepare folders
mkdir -p "$OUT_ROOT"/{original.fastq/,clean.fastq/,taxonomy/,logs/}
LOG_DIR="$OUT_ROOT/logs"
CLEAN_DIR="$OUT_ROOT/clean.fastq"
ORIG_DIR="$OUT_ROOT/original.fastq"
TAX_DIR="$OUT_ROOT/taxonomy"

# ------------------------------
# 0) Gather & normalize inputs
# ------------------------------
shopt -s nullglob
#cd "$RAW_DIR" || die "Cannot cd RAW_DIR: $RAW_DIR"

# DECOMPRESSION
# checks for compressed files '*.gz' (discarding the output by funneling to /dev/null)
# compgen -G exit status will be 0 if at least one is found, executing if statement block
# Ensures pigz decompression does not error out if there are no compressed files
if compgen -G "$RAW_DIR/*.gz" > /dev/null; then
  log "Decompressing .gz with pigz…"
  # pigz is parallel implementation of gzip that is faster on multi-core systems
  # '--' ensures that filenames starting with '-' are not treated as flags on execution
  pigz -d -- $RAW_DIR/*.gz
fi

# SAMPLE CATALOGING
# Constructs arrays of fastq file names and saves them to array variables R1S and R2S
# '2>/dev/null' supresses errors if no matching file names are found
# file names are sorted alphabetically before being saved to array variable 
mapfile -t R1S < <(ls $RAW_DIR/*_R1_001.fastq 2>/dev/null | sort)
mapfile -t R2S < <(ls $RAW_DIR/*_R2_001.fastq 2>/dev/null | sort)
# Checks that length of R1S > 0 (gt = greather than), else dies
[[ ${#R1S[@]} -gt 0 ]] || die "No R1 fastqs found in $RAW_DIR"
# Checks that length RS1 = length RS2 (-eq = equals)
[[ ${#R1S[@]} -eq ${#R2S[@]} ]] || die "R1/R2 count mismatch"

# Derive sample basenames by stripping file suffixes, stored in 'samples' array
samples=()
for r1 in "${R1S[@]}"; do
  # strip to just file name (remove directory path)
  fname=${r1##*/}
  # removes filename suffix
  base=${fname%_R1_001.fastq}
  # checks that correlating sample file is in R2
  [[ -f "$RAW_DIR/${base}_R2_001.fastq" ]] || die "Missing R2 for ${base}"
  samples+=("$base")
  # archive originals
  # creates a soft link for easy access of original data files
  ln -sf "$RAW_DIR/${base}_R1_001.fastq" "$ORIG_DIR/"
  ln -sf "$RAW_DIR/${base}_R2_001.fastq" "$ORIG_DIR/"
done

log "Found ${#samples[@]} samples"

# Record sample names to sample.names in output directory
printf "%s\n" "${samples[@]}" > "$OUT_ROOT/sample.names"

# ------------------------------
# 1) Extract R2 UMIs and select corresponding R1 reads
#       GTP's Why: Parse UMI as 16 nt immediately upstream of primer on R2, 
#       avoid brittle cut/sed.
# ------------------------------

extract_umi_and_select_r1() {
  # sets the variable 'base' to the first positioned parameter passed
  local base="$1"
  log "[$base] UMI extraction from R2"

  # construct r1 and r2 filenames from base
  # subsitutes the string stored in base into fastq file names a
  local r1="$RAW_DIR/${base}_R1_001.fastq"
  local r2="$RAW_DIR/${base}_R2_001.fastq"

  # output file path contruction
  local umi_map="$CLEAN_DIR/UMI.${base}.tsv"        # NAME\tUMI
  local sel_names="$CLEAN_DIR/selected.${base}.names" # R1 read names to keep
  local sel_r1="$CLEAN_DIR/selected.${base}.R1.fastq"

  # extracting UMI using awk utlitiy (separate language = awk sepcializing in text parsing) 
  # Passed Variables: Primer Motif and UMI Length
  # Outputs: umi_map.tsv (2 columns = NAME<tab>UMI) and selected.bases.names (1 column = NAME)
  # Name Values = name of reads with detectable UMI+motif pattern in R2 
  #     **looks like name is just be header line from read file (following '@')**
  # UMI Values = unique molecular identifiers of length $UMI_LEN in alphabet [ACGTN]
  # **Utlizing UMI to collapse reads results in a low read count, making analysis difficult**
  awk -v motif="$PRIMER_MOTIF" -v L=$UMI_LEN '
    NR%4==1 { hdr=$0; name=substr($0,2); next } 
    NR%4==2 { seq=$0; 
              if (match(seq, "([ACGTN]{"L"})" motif, m)) {
                print name"\t"m[1]
              }
            }
    ' "$r2" | tee "$umi_map" | cut -f1 > "$sel_names"

  log "[$base] Subsetting R1 to names selected based on motif and UMI identification."

  # Subsets R1 to those entries for which we successfully extracted a UMI in R2
  #     **Why are we not also subsetting R2???**
  #         Amiran said that R2 are only usable as a means of detecting UMI and
  #         are unusable otherwise. (Why??)
  seqtk subseq "$r1" "$sel_names" > "$sel_r1"
}

# SWARM FILE CONSTRUCTION (optional, based on EMIT_SWARm global variable)
#   **construct swarm files to execute up through R1 subsetting based on UMIs**
#   **return to this for understanding**
if [[ "$EMIT_SWARM" == "1" ]]; then
  SW1="$OUT_ROOT/01_umi_select.swarm"
  : > "$SW1"
  for s in "${samples[@]}"; do
    echo "awk -v motif=$PRIMER_MOTIF -v L=$UMI_LEN 'NR%4==1{hdr=\$0; name=substr(\$0,2); next} NR%4==2{seq=\$0; if(match(seq, \"([ACGTN]{\"L\"})\"motif, m)){print name\"\\t\"m[1]}}' ${s}_R2_001.fastq | tee $CLEAN_DIR/UMI.${s}.tsv | cut -f1 > $CLEAN_DIR/selected.${s}.names; seqtk subseq ${s}_R1_001.fastq $CLEAN_DIR/selected.${s}.names > $CLEAN_DIR/selected.${s}.R1.fastq" >> "$SW1"
  done
  log "Wrote swarm file: $SW1"
else
  for s in "${samples[@]}"; do extract_umi_and_select_r1 "$s"; done
fi

# ------------------------------
# 2) Trim selected R1 (single-end)
#       GTP's Why: Matches your AVGQUAL:30 logic; adapters optional.
# ------------------------------
trim_selected_r1() {
  local base="$1"
  local in_fq="$CLEAN_DIR/selected.${base}.R1.fastq"
  # TODO: change name of trimmed fastq to match previous naming convention
  local out_fq="$CLEAN_DIR/${base}.trimmed.R1.fastq"

  # Check to ensure input fastq exists
  #   Changed to exists and file flag -f (was -s which checks if exists and is nonempty)
  [[ -f "$in_fq" ]] || die "Missing FASTQ File for Trimming: $in_fq"


  log "[$base] Trimming with Trimmomatic"
  # Trimming R1 fastq's using trimmomatic (java based processing and trimming tool for Illumina tech)
  #     Flags:
  #       SE = single end reads
  #       -phred33 = phred+33 encoding of quality (most likely with modern Illumina systems)
  #     Trimmomatic Filtering Arguments:
  #       AVGQUAL:30 (given as global variable)
  #     
  #     TODO: Confirm that our fastq's use Phred+33 encoding
  java -jar "$TRIMMOMATIC_JAR" SE -threads "$THREADS" -phred33 "$in_fq" "$out_fq" "${TRIM_ARGS[@]}"
}

# SWARM FILE CONSTRUCTION (optional, based on EMIT_SWARM global variable)
#   **return to this for understanding**
if [[ "$EMIT_SWARM" == "1" ]]; then
  SW2="$OUT_ROOT/02_trim.swarm"; : > "$SW2"
  for s in "${samples[@]}"; do
    echo "java -jar $TRIMMOMATIC_JAR SE -threads $THREADS $CLEAN_DIR/selected.${s}.R1.fastq $CLEAN_DIR/${s}.trimmed.R1.fastq ${TRIM_ARGS[*]}" >> "$SW2"
  done
  log "Wrote swarm file: $SW2"
else
  for s in "${samples[@]}"; do trim_selected_r1 "$s"; done
fi

# ------------------------------
# 3) FASTQ -> FASTA for mapping (to mirror your original behavior)
#           **Why doing this??**
# ------------------------------
log "Converting trimmed FASTQ to FASTA"
for s in "${samples[@]}"; do
  seqkit fq2fa "$CLEAN_DIR/${s}.trimmed.R1.fastq" > "$CLEAN_DIR/${s}.trimmed.R1.fasta"
done



# ------------------------------
# 4) Remove human, then remove viral; keep unmapped names each time
#    Why: Use SAM flags, not "grep NC_". Reliable across references.
# ------------------------------

# REFERENCE INDEXING
# BWA index prefixes (must exist: *.bwt etc.)
REFS=(
  "$HUMAN_REF"
  "$VIRAL_REF"
  "$BACT16S_REF"
)
for ref in "${REFS[@]}"; do
  ref_name="${ref##*/}"

  if [[ ! -f "${ref}.bwt" ]]; then
    log "[$ref_name] indexing with bwa"
    bwa index "$ref"
  fi
done


map_and_keep_unmapped_names() {
  local base="$1" ref="$2" in_fa="$3" tag="$4"  # tag: hum|viral
  local sam="$CLEAN_DIR/${base}.${tag}.sam"
  local unmapped_names="$CLEAN_DIR/${base}.${tag}.unmapped.names"

  log "[$base] bwa mem vs $tag"

  bwa mem -t "$THREADS" "$ref" "$in_fa" > "$sam"

  # Filters reads in the output SAM file to those that did not align based on bwa mem flagging
  # Flag = 4 -> the query sequence itself is unmapped
  # Save names of unmapped reads to text file (specific to unmapped with human or virus)
  samtools view -f 4 "$sam" | cut -f1 | sort -u > "$unmapped_names"
  printf "%s\n" "$sam" "$unmapped_names"
}

# Filter reads from input fastq to those unmapped after alignment to human genome
nonhuman_fa() {
  local base="$1"
  local in_fa="$CLEAN_DIR/${base}.trimmed.R1.fasta"
  map_and_keep_unmapped_names "$base" "$HUMAN_REF" "$in_fa" hum >/dev/null
  seqtk subseq "$in_fa" "$CLEAN_DIR/${base}.hum.unmapped.names" > "$CLEAN_DIR/${base}.nonhuman.fasta"
}

# Further filter reads from previously filtered nonhuman reads down to those also unmapped after alignment to viral genomes
nonviral_fa() {
  local base="$1"
  local in_fa="$CLEAN_DIR/${base}.nonhuman.fasta"
  map_and_keep_unmapped_names "$base" "$VIRAL_REF" "$in_fa" viral >/dev/null
  log "[$base] merging nonviral-nonhuman reads)"
  seqtk subseq "$in_fa" "$CLEAN_DIR/${base}.viral.unmapped.names" > "$CLEAN_DIR/${base}.nonhuman.nonviral.fasta"
}

# TODO: onstruct order-independent means of merging nonhuman and nonviral read sets
# **Has to be a more efficient way of forming union between the two unmapped read sets**
for s in "${samples[@]}"; do
  nonhuman_fa "$s"
  nonviral_fa "$s"
done

# ------------------------------
# 5) Map to bacterial 16S rRNA reference and keep mapped names (positives)
# ------------------------------
for s in "${samples[@]}"; do
  log "Mapping $s to bacterial 16S rRNA reference"
  bwa mem -t "$THREADS" "$BACT16S_REF" "$CLEAN_DIR/${s}.nonhuman.nonviral.fasta" > "$CLEAN_DIR/${s}.bact.sam"
  
  # SAMs are tab delimited 
  # $1 statement skips header lines (first character == @)
  # $3 statement selects reads that are unmapped (RNAME in alignment records != *)
  # prints $1 = QNAME (read name in alignment records), stores in awk output
  # sorts those read names output piped from awk
  # saves sorted read names
  awk '
    ($1 !~ /^@/) && ($3 != "*") {
      print $1
    }
  ' "$CLEAN_DIR/${s}.bact.sam" \
  | sort -u > "$CLEAN_DIR/${s}.bact.pos.names"

  # Filtering all reads down to those identified via 16S alignment
  #   TODO: need to change output to fit better naming convention
  seqtk subseq "$CLEAN_DIR/${s}.nonhuman.nonviral.fasta" "$CLEAN_DIR/${s}.bact.pos.names" > "$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
done

  ######### PIPELINE CURRENTLY FUNCTION UP TO THIS POINT #########
# ------------------------------
# 6) Stats (robust counts)
# ------------------------------
log "Computing count stats"
{
  echo -e "ID\tTotal_reads\tHuman_mapped_reads\tClean_reads"
  for s in "${samples[@]}"; do
    # Determining Number of reads in original fastq 
    # (each read entry has 4 lines: @Header, Sequence, +Header, Quality), hence divide by 4
    # TODO: ensure that raw data is not being altered after original data link was created
    #       originally used the original data root but was having issues with linking
    total=$(($(wc -l < "$RAW_DIR/${s}_R1_001.fastq")/4))

    # Counts number of reads mapped to human genome
    # TODO: make switchable to other reference organisms (namely mouse)
    # TODO: look into possibility that secondary, supplementary, split reads possibly inflating this value
    hum_mapped=$(samtools view -F 4 "$CLEAN_DIR/${s}.hum.sam" | wc -l)

    # Counts number of reads in cleaned nonhuman and nonviral fasta
    #   If the fasta file is not found or there are no reads in the fasta file,
    #   still exits with true status (clean = 0) avoiding erroring out
    clean=$(grep -c '^>' "$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta" || true)

    # Pringts a tab-delimited entry for the read counts for the sample to output
    echo -e "${s}\t${total}\t${hum_mapped}\t${clean}"
  done
} > "$OUT_ROOT/xin1.stats.tsv"  # redirects all stdout from this block (header + loop) into a TSV file


# ------------------------------
# 7) OPTIONAL: Negative control DB & filtering (keeps reads dissimilar to neg DB)
# ------------------------------
# Inputs: any sample name containing strings NEGATIVECONTROL|CELLSCONTROL|NEGWATER|BLANK|NEGcDNA|NEGPCR (adjust below)

NEG_TAGS=(
  "NEGATIVECONTROL"
  "CELLSCONTROL"
  "NEGWATER"
  "BLANK"
  "NEGcDNA"
  "NEGPCR"
  "NEG"
  )

build_neg_db_and_filter() {
  log "Building negative-control DB (top $NEG_DB_TOP_K uniques)"

  # enabling extended globbing patterns for ls to recoginze negative control file names
  shopt -s extglob

  local -a neg_all=()
  local -a neg_tag=()

  #collects all candidate neg files (across all tags)
  for tag in "${NEG_TAGS[@]}"; do

    # TODO: makesure that formatted string works as intended
    # CURRENTLY RUNNING INTO SYNTAX ERROR HERE
    neg_tag=( "$CLEAN_DIR"/clean.*_"${tag}"?([0-9])_*.R1.nonhuman.nonviral.fasta )
    log "DEBUG: files for tag '$tag':"
    printf '  - %q\n' "${neg_tag[@]}" >&2

    #builds a per-tag composite fasta if any are found for this tag
    if ((${#neg_tag[@]} > 0)); then
      # deduplicating in case ls found duplicate files (consider removing)
      mapfile -t neg_tag < <(printf "%s\n" "${neg_tag[@]}" | sort -u)

      log "Found ${#neg_tag[@]} negative-control FASTA(s) for conrtol tag '$tag'"



      
      # concatenate all negative control files found with this tag
      cat "${neg_tag[@]}" > "$CLEAN_DIR/combined.${tag}.neg.fasta"

      # also add all of the files found with this tag to list of all negative control fasta's
      neg_all+=( "${neg_tag[@]}" )

    else
      log "No FASTA matched negative-control tag: '$tag'"
    fi
  done

  # disabling extended globbing patterns for ls
  shopt -u extglob

  # deduplicates the overall negative control file list (consider removing)
  if ((${#neg_all[@]} > 0)); then
    mapfile -t neg_all < <(printf "%s\n" "${neg_all[@]}" | sort -u)
  
  #if no negative controls are found -> skip
  else
    log "No negative control FASTA found (skipping negative control filtering)"
    return 0
  fi

  # Build overal negative control composite (all tags)
  log "Building combined negative-control fasta from ${#neg_all[@]} file(s)"
  log "DEBUG: files for combined.neg.fasta:"
  printf '  - %q\n' "${neg_all[@]}" >&2

  cat "${neg_all[@]}" > "$CLEAN_DIR/combined.neg.fasta"

  # uses Mothur to collapse duplicate reads in the combined negative control fasta
  #   this generates a new fasta with the name convention <original_name>.unique.fasta (stored in pwd)
  # redirects the stderr from the mothur command (mothur's log) to the log directory
  mothur "#unique.seqs(fasta=$CLEAN_DIR/combined.neg.fasta)" 2>"$LOG_DIR/mothur.neg.unique.log"

  # takes the combined negative control
  cut -f1 "$CLEAN_DIR/combined.neg.count_table" | tail -n +2 \
    | paste - <(cut -f2 "$CLEAN_DIR/combined.neg.count_table" | tail -n +2) \
    | sort -k2,2nr | head -n "$NEG_DB_TOP_K" | cut -f1 > "$CLEAN_DIR/top.$NEG_DB_TOP_K.neg.names"
  
  seqtk subseq "$CLEAN_DIR/combined.neg.unique.fasta" "$CLEAN_DIR/top.$NEG_DB_TOP_K.neg.names" > "$CLEAN_DIR/top.neg.db.fasta"
  bwa index "$CLEAN_DIR/top.neg.db.fasta"

  for s in "${samples[@]}"; do
    in_fa="$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
    [[ -s "$in_fa" ]] || continue
    bwa mem -t "$THREADS" "$CLEAN_DIR/top.neg.db.fasta" "$in_fa" > "$CLEAN_DIR/${s}.neg.sam"
    # keep names with low AS (not similar to negative DB)
    awk -v th="$NEG_AS_KEEP_LT" 'BEGIN{FS="\t"} $0 !~ /^@/ { m=match($0, /AS:i:([0-9]+)/, a); if(m && a[1] < th) print $1 }' "$CLEAN_DIR/${s}.neg.sam" | sort -u > "$CLEAN_DIR/${s}.selected.names"
    seqtk subseq "$in_fa" "$CLEAN_DIR/${s}.selected.names" > "$CLEAN_DIR/${s}.selected.fasta"
  done
}

OLD_build_neg_db_and_filter() {
  ###### OLD VERSION OF THIS FUNCTION ##########
  # finds and lists negative control files that match the "$CLEAN_DIR"/CLEAN.*NEG*.R1.nonhumna,nonviral.fasta
  # stores these final names in the new array 'neg-fa'
  # 2>/dev/null || true ensure that this continues even if no negative control files are found
  #   TODO: check to see if this is exclusive to literal 'NEG' instead of all possible negative strings
  #         like those listed above in NEG_STRINGS (not seeing any collapsed fasta after mothur unique.seq)
  mapfile -t neg_fa < <(ls "$CLEAN_DIR"/clean.*NEG*.R1.nonhuman.nonviral.fasta 2>/dev/null || true)

  # checks if indeed, any negative control files were round and exits if not (len neg_fa > 0)
  [[ ${#neg_fa[@]} -gt 0 ]] || { log "No negative-control FASTA found (skipping)"; return 0; }

  # combines all negative control fasta files into a single combined fasta file
  # TODO: built combined fasta on per negative control type and combine all negative control for full purity filters
  cat "${neg_fa[@]}" > "$CLEAN_DIR/combined.neg.fasta"

  # uses Mother to collapse duplicate reads in the combined negative control fasta
  #   this generates a new fasta with the name convention <original_name>.unique.fasta (stored in pwd)
  # redirects the stderr from the mothur command (mothur's log) to the log directory
  mothur "#unique.seqs(fasta=$CLEAN_DIR/combined.neg.fasta)" 2>"$LOG_DIR/mothur.neg.unique.log"

  # takes the combined negative control
  cut -f1 "$CLEAN_DIR/combined.neg.count_table" | tail -n +2 \
    | paste - <(cut -f2 "$CLEAN_DIR/combined.neg.count_table" | tail -n +2) \
    | sort -k2,2nr | head -n "$NEG_DB_TOP_K" | cut -f1 > "$CLEAN_DIR/top.$NEG_DB_TOP_K.neg.names"
  
  seqtk subseq "$CLEAN_DIR/combined.neg.unique.fasta" "$CLEAN_DIR/top.$NEG_DB_TOP_K.neg.names" > "$CLEAN_DIR/top.neg.db.fasta"
  bwa index "$CLEAN_DIR/top.neg.db.fasta"

  for s in "${samples[@]}"; do
    local in_fa="$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
    [[ -s "$in_fa" ]] || continue
    bwa mem -t "$THREADS" "$CLEAN_DIR/top.neg.db.fasta" "$in_fa" > "$CLEAN_DIR/${s}.neg.sam"
    # keep names with low AS (not similar to negative DB)
    awk -v th="$NEG_AS_KEEP_LT" 'BEGIN{FS="\t"} $0 !~ /^@/ { m=match($0, /AS:i:([0-9]+)/, a); if(m && a[1] < th) print $1 }' "$CLEAN_DIR/${s}.neg.sam" | sort -u > "$CLEAN_DIR/${s}.selected.names"
    seqtk subseq "$in_fa" "$CLEAN_DIR/${s}.selected.names" > "$CLEAN_DIR/${s}.selected.fasta"
  done
}

# execute negative control filtering
build_neg_db_and_filter

# ------------------------------
# 8) Mothur classify (genus/species/class/phylum summaries)
# ------------------------------
classify_and_summarize() {
  log "mothur classify.seqs on final FASTA (or .selected.fasta if present)"
  pushd "$TAX_DIR" >/dev/null
  for s in "${samples[@]}"; do
    src="$CLEAN_DIR/${s}.selected.fasta"
    [[ -s "$src" ]] || src="$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
    [[ -s "$src" ]] || { log "[$s] no FASTA to classify"; continue; }

    outprefix="${s}.final"
    mothur "#classify.seqs(fasta=$src, template=$MOTHUR_TEMPLATE, taxonomy=$MOTHUR_TAX, method=wang, cutoff=0, processors=$THREADS)" \
      2>"$LOG_DIR/mothur.${s}.log"

    tax="$(basename "$src").ncbi20.wang.taxonomy"
    [[ -f "$tax" ]] || tax="${outprefix}.ncbi20.wang.taxonomy" # fallback if mothur renames
    [[ -f "$tax" ]] || { log "[$s] taxonomy file not found"; continue; }

    # strip confidences, make rank-specific counts
    sed -E 's/\([^()]*\)//g' "$tax" | awk -F"\t" '{print $2}' > "${outprefix}.lineages"

    for level in genus species class phylum; do
      case $level in
        genus)   cutf='-d";" -f2,3,4,5,6' ; min=10 ;;
        species) cutf='-d";" -f2,3,4,5,6,7' ; min=5 ;;
        class)   cutf='-d";" -f2,3' ; min=50 ;;
        phylum)  cutf='-d";" -f2' ; min=50 ;;
      esac
      eval "cut $cutf \"${outprefix}.lineages\"" | sort | uniq -c | awk -v m="$min" '$1>m{print $1, $2}' | sed -E 's/^ +//' | sort -k2 > "${outprefix}.${level}.fin.count"
    done
  done

  # Combine across samples: genus/species/class/phylum
  combine_level() {
    local level="$1"
    # collect IDs
    ls *.${level}.fin.count | head -n1 | xargs -I{} awk '{print $2}' {} > "taxa.${level}.names"
    echo ID | cat - "taxa.${level}.names" > "header.${level}"
    # columns
    tmpcols=()
    for s in "${samples[@]}"; do
      f="${s}.final.${level}.fin.count"
      [[ -s "$f" ]] || continue
      echo "$s" > ".id.${level}.${s}"; awk '{print $1}' "$f" > ".col.${level}.${s}"
      paste ".id.${level}.${s}" ".col.${level}.${s}" > ".named.${level}.${s}"
      tmpcols+=(".named.${level}.${s}")
    done
    paste "header.${level}" "${tmpcols[@]}" > "xin1.final.${level}.count.table"
  }
  combine_level genus
  combine_level species
  combine_level class
  combine_level phylum
  popd >/dev/null
}

#execute classification and summarization script
#classify_and_summarize

log "Done. Outputs in: $OUT_ROOT"
