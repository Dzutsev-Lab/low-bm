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
# enabling extended globbing patterns for ls to recoginze negative control file names
shopt -s extglob
IFS=$'\n\t'

module load trimmomatic seqkit seqtk bwa samtools mothur pigz java

# ------------------------------
# CONFIG — EDIT ME
# ------------------------------
RAW_DIR="./Select_Trial_Data"               # input FASTQs
OUT_ROOT="./Select_Trial_Output"          # output root
THREADS="2"
PRIMER_MOTIF="GGACTAC"                                          # 16nt UMI is immediately upstream of this
UMI_LEN="16"

# BWA index prefixes (must exist: *.bwt etc.)
HUMAN_REF="/data/dzutseva/refdata/refseq/clean.human.fa"        # or your Mus reference for mouse sets
VIRAL_REF="/data/dzutseva/refdata/idexes/all.viral.fna"
BACT16S_REF="/data/Trinchieri_lab/dzutseva/genomes/all.rrna.bacteria/all.rrna.bacteria"

# mothur references
MOTHUR_TEMPLATE="/data/Trinchieri_lab/jenny/ncbi20.fasta"
MOTHUR_TAX="/data/Trinchieri_lab/jenny/ncbi20.tax"

# Trimming (match your historical settings)
TRIMMOMATIC_JAR="$TRIMMOJAR"                                   # set by module trimmomatic
TRIM_ARGS=("AVGQUAL:30")                                        # keep semantics; optionally add ILLUMINACLIP/MINLEN/CROP

# Negative control filtering (optional)
NEG_DB_TOP_K="10000"                                            # top N negative-control uniques to include
NEG_AS_KEEP_LT="160"                                            # keep if AS < this (i.e., not similar to neg DB)

# Swarm toggle: if set to 1, writes swarm files alongside loops (does not submit)
EMIT_SWARM="0"

# ------------------------------
# Helpers
# ------------------------------

# LOG: prints a formatted string including date/time and status update on current processes
# Input: some record to be logged 
#
log() { printf "[%s] %s\n" "$(date +"%F %T")" "$*"; }
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
  # WHY ARE WE DOING THIS AND DO I NEED TO BE MORE AWARE OF 
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
    NR%4==1 { 
        hdr=$0; 
        name=substr($0,2)
        sub(/ .*/, "", name)
        next 
    }
    NR%4==2 { 
        seq=$0
        if (match(seq, "([ACGTN]{" L "})" motif, m)) {
                    print name"\t"m[1]
        }
    }
    ' "$r2" | tee "$umi_map" | cut -f1 > "$sel_names"

  log "[$base] Subsetting R1 to names selected based on motif and UMI identification."

  # Subsets R1 to those entries for which we successfully extracted a UMI in R2
  #     **Why are we not also subsetting R2???**
  #         Amiran
  seqtk subseq "$r1" "$sel_names" > "$sel_r1"
}

# SWARM FILE CONSTRUCTION (optional, based on EMIT_SWARm global variable)
#   **seems to only construct swarm files to execute up through R1 subsetting based on UMIs**
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

  ######### PIPELINE CURRENTLY FUNCTION UP TO THIS POINT #########

# ------------------------------
# 4) Remove human, then remove viral; keep unmapped names each time
#    Why: Use SAM flags, not "grep NC_". Reliable across references.
# ------------------------------

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
  awk \
  '($1 !~ /^@/) && '\
  '($3 != "*") '\
  '{print $1}' \
  "$CLEAN_DIR/${s}.bact.sam" \
  | sort -u > "$CLEAN_DIR/${s}.bact.pos.names"

  # Filtering all reads down to those identified via 16S alignment
  #   TODO: need to change output to fit better naming convention
  seqtk subseq "$CLEAN_DIR/${s}.nonhuman.nonviral.fasta" "$CLEAN_DIR/${s}.bact.pos.names" > "$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
done

# ------------------------------
# 6) Stats (robust counts)
# ------------------------------
log "Computing count stats"
{
  echo -e "ID\tTotal_reads\tHuman_mapped_reads\tClean_reads"
  for s in "${samples[@]}"; do
    total=$(($(wc -l < "$ORIG_DIR/${s}_R1_001.fastq")/4))
    hum_mapped=$(samtools view -F 4 "$CLEAN_DIR/${s}.hum.sam" | wc -l)
    clean=$(grep -c '^>' "$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta" || true)
    echo -e "${s}\t${total}\t${hum_mapped}\t${clean}"
  done
} > "$OUT_ROOT/xin1.stats.tsv"

# ------------------------------
# 7) OPTIONAL: Negative control DB & filtering (keeps reads dissimilar to neg DB)
# ------------------------------
# Inputs: any sample name containing strings NEGATIVECONTROL|CELLSCONTROL|NEGWATER|BLANK|NEGcDNA|NEGPCR (adjust below)

NEG_STRINGS='NEGATIVECONTROL|CELLSCONTROL|NEGWATER|BLANK|NEGcDNA|NEGPCR|NEG'

build_neg_db() {
  log "Building negative-control DB (top $NEG_DB_TOP_K uniques)"
  mapfile -t neg_fa < <(ls "$CLEAN_DIR"/clean.*NEG*.R1.nonhuman.nonviral.fasta 2>/dev/null || true)
  [[ ${#neg_fa[@]} -gt 0 ]] || { log "No negative-control FASTA found (skipping)"; return 0; }

  cat "${neg_fa[@]}" > "$CLEAN_DIR/combined.neg.fasta"
  mothur "#unique.seqs(fasta=$CLEAN_DIR/combined.neg.fasta)" 2>"$LOG_DIR/mothur.neg.unique.log"
  cut -f1 "$CLEAN_DIR/combined.neg.count_table" | tail -n +2 \
    | paste - <(cut -f2 "$CLEAN_DIR/combined.neg.count_table" | tail -n +2) \
    | sort -k2,2nr | head -n "$NEG_DB_TOP_K" | cut -f1 > "$CLEAN_DIR/top.$NEG_DB_TOP_K.neg.names"
  seqtk subseq "$CLEAN_DIR/combined.neg.unique.fasta" "$CLEAN_DIR/top.$NEG_DB_TOP_K.neg.names" > "$CLEAN_DIR/top.neg.db.fasta"
  # Indexes the top negative control reads fasta
  bwa index "$CLEAN_DIR/top.neg.db.fasta"

  # filters the reads in other samples based on the top reads found in the negative controls
  # skipping for now (will do negative control filtering downstream)
  # for s in "${samples[@]}"; do
  #   in_fa="$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
  #   [[ -s "$in_fa" ]] || continue
  #   bwa mem -t "$THREADS" "$CLEAN_DIR/top.neg.db.fasta" "$in_fa" > "$CLEAN_DIR/${s}.neg.sam"
  #   # keep names with low AS (not similar to negative DB)
  #   awk -v th="$NEG_AS_KEEP_LT" 'BEGIN{FS="\t"} $0 !~ /^@/ { m=match($0, /AS:i:([0-9]+)/, a); if(m && a[1] < th) print $1 }' "$CLEAN_DIR/${s}.neg.sam" | sort -u > "$CLEAN_DIR/${s}.selected.names"
  #   seqtk subseq "$in_fa" "$CLEAN_DIR/${s}.selected.names" > "$CLEAN_DIR/${s}.selected.fasta"
  # done
}

# Uncomment to enable negative-control filtering
# build_neg_db_and_filter

# ------------------------------
# 8) Mothur classify (genus/species/class/phylum summaries)
# ------------------------------
classify_and_summarize() {
  log "mothur classify.seqs on final FASTA (or .selected.fasta if present)"

  mkdir -p "$TAX_DIR"

  for s in "${samples[@]}"; do
    # input fasta (cleaned)
    src="$CLEAN_DIR/clean.${s}.R1.nonhuman.nonviral.fasta"
    [[ -s "$src" ]] || { log "[$s] no FASTA to classify"; continue; }

    #set output file name prefix
    outprefix="$TAX_DIR/${s}.final"

    #classify reads in sample file using wang method 
    #   cutoff set to 0 = including all possible assignments (should include confidence score for each assignment)
    mothur "#set.dir(output=$TAX_DIR); classify.seqs(fasta=$src, template=$MOTHUR_TEMPLATE, taxonomy=$MOTHUR_TAX, method=wang, cutoff=0, processors=$THREADS)" \
      2>"$LOG_DIR/mothur.${s}.log"

    # constructs expected taxonomy file for the sample
    tax="$TAX_DIR/$(basename "$src").ncbi20.wang.taxonomy"
    [[ -f "$tax" ]] || tax="${outprefix}.ncbi20.wang.taxonomy" # fallback if mothur renames
    [[ -f "$tax" ]] || { log "[$s] taxonomy file not found"; continue; }

    # extracts taxonomic lineage strings from Mothur output 
    #   strips confidence scores and sequence names, single column
    sed -E 's/\([^()]*\)//g' "$tax" | awk -F"\t" '{print $2}' > "${outprefix}.lineages"

    for level in genus species class phylum; do
      # sets specific minimum count cutoffs for different levels of classification
      case $level in
        genus)   cutf='-d";" -f2,3,4,5,6' ; min=10 ;;
        species) cutf='-d";" -f2,3,4,5,6,7' ; min=5 ;;
        class)   cutf='-d";" -f2,3' ; min=50 ;;
        phylum)  cutf='-d";" -f2' ; min=50 ;;
      esac
      # consolidates identical lineage strings and counts the number of occurances of each lineage
      # filters out taxa not meeting count cut off for the level
      # totals output to level-specific count file
      eval "cut $cutf \"${outprefix}.lineages\"" \
      | sort \
      | uniq -c \
      | awk -v m="$min" '$1>m{print $1, $2}' \
      | sed -E 's/^ +//' \
      | sort -k2 \
      > "${outprefix}.${level}.fin.count"
    done
  done

  # Combine across samples: genus/species/class/phylum
  combine_level() {
    local level="$1"
    # collect IDs
    ls "$TAX_DIR"*.${level}.fin.count | head -n1 | xargs -I{} awk '{print $2}' {} > "taxa.${level}.names"
    echo ID | cat - "taxa.${level}.names" > "header.${level}"
    # columns
    tmpcols=()
    for s in "${samples[@]}"; do
      f="$TAX_DIR/${s}.final.${level}.fin.count"
      [[ -s "$f" ]] || continue
      echo "$s" > "$TAX_DIR/.id.${level}.${s}"
      awk '{print $1}' "$f" > "$TAX_DIR/.col.${level}.${s}"
      paste "$TAX_DIR/.id.${level}.${s}" "$TAX_DIR/.col.${level}.${s}" > "$TAX_DIR/.named.${level}.${s}"
      tmpcols+=("$TAX_DIR/.named.${level}.${s}")
    done
    paste "$TAX_DIR/header.${level}" "${tmpcols[@]}" > "$TAX_DIR/xin1.final.${level}.count.table"
  }
  combine_level genus
  combine_level species
  combine_level class
  combine_level phylum
}

# Uncomment to classify now
# classify_and_summarize

log "Done. Outputs in: $OUT_ROOT"
