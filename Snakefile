import json
import os
import sys
from pathlib import Path

# Snakemake's --directory changes the process working directory, and SLURM
# execution can evaluate sources from Snakemake's runtime cache. The low-bm
# launcher injects the real checkout path so config paths stay repo-relative.
REPO_ROOT = Path(config.get("low_bm_repo_root", Path.cwd())).resolve()


def repo_path(value):
    text = os.path.expandvars(os.path.expanduser(str(value)))
    path = Path(text)
    if path.is_absolute():
        return str(path)
    return str(REPO_ROOT / path)


# Input and Output Directories
TRIAL_ID = config["trialID"]
TRIAL_NAME = str(TRIAL_ID) + "_" + config["trial_descript"]
EXP_DIR = config["exp_dir"]

IN_ROOT = repo_path(config["in_root"])
IP_ROOT = repo_path(config["ip_root"])
OUT_ROOT = repo_path(config["out_root"])

IN_DIR =    os.path.join(IN_ROOT, EXP_DIR)
IP_DIR =    os.path.join(IP_ROOT, EXP_DIR)
OUT_DIR =   os.path.join(OUT_ROOT, TRIAL_NAME)
METADATA =  os.path.join(IN_ROOT, config["metadata"])

SCRIPTS = repo_path(config["script_dir"])
REF_DIR = repo_path(config.get("ref_dir", "Ref_Data"))
CONDA_ENV_DIR = repo_path(config.get("conda_env_dir", "workflow/envs"))
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)
from SamplePrep import build_fastq_manifest


def conda_env(name):
    return f"{CONDA_ENV_DIR}/{name}.yaml"


RAW = f"{IN_DIR}/OUTPUT/Data/fastq"
NORM_RAW_DIR = f"{IP_DIR}/00_RawNorm"
UMI_SELECT_DIR = f"{IP_DIR}/01_UMISelection"
UMI_DEDUP_DIR = f"{IP_DIR}/02_UMIDeduplication"
DADA_DENOISE_DIR = f"{IP_DIR}/03_Dada_Denoising"
NEG_ALIGNMENT_DIR = f"{IP_DIR}/04_Neg_Alignment"
NEG_ALIGN_FILT_DIR = f"{IP_DIR}/05_Neg_Alignment_Filter"
POS_ALIGNMENT_DIR = f"{IP_DIR}/06_Pos_Alignment"
MOTHUR_TAX_DIR = f"{IP_DIR}/08.1_Mothur_Taxonomy"
KRAKEN_TAX_DIR = f"{IP_DIR}/08.2_Kraken_Taxonomy"
BLAST_TAX_DIR = f"{IP_DIR}/08.3_BLAST_Taxonomy"
RECONCILED_TAX_DIR = f"{IP_DIR}/08.4_Reconciled_Taxonomy"
NORM_COUNT_DIR = f"{OUT_DIR}/CountNormalization"
PHYLOSEQ_DIR = f"{OUT_DIR}"
DIFF_ABUND_DIR = f"{OUT_DIR}"

TRACK_DIR = f"{OUT_DIR}/Tracking"
LOG_DIR = f"{OUT_DIR}/Logs"

# Sequencing metadata
def parse_bool(value, default=True, name="value"):
    if value is None:
        return default
    if isinstance(value, bool):
        return value

    text = str(value).strip().lower()
    if text == "":
        return default
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    raise ValueError(f"Invalid boolean for {name}: {value!r}")

PROCESS_UMIS = parse_bool(config.get("process_umis"), default=True, name="process_umis")
config["process_umis"] = PROCESS_UMIS

R1_PRIMER = config["r1_primer"]
R2_PRIMER = config["r2_primer"]
R1_PRIMER_MOTIF_LEN = config["r1_primer_motif_len"]
R2_PRIMER_MOTIF_LEN = config["r2_primer_motif_len"]
R2_PRIMER_SKIP = config["r2_primer_skip"]
POLY_G_THRESHOLD = config["poly_G_threshold"]
UMI_LEN = config["umi_len"]
MAX_OFFSET = config["max_offset"]

# Reference Data
HOST = str(config.get("host", "")).strip().lower()
if HOST not in {"human", "mouse"}:
    raise ValueError("Config key 'host' must be either 'human' or 'mouse'. Set it in the per-batch row config.")
config["host"] = HOST
if HOST == "human":
    HOST_REF = repo_path(config["human_ref"])
elif HOST == "mouse":
    HOST_REF = repo_path(config["mouse_ref"])

VIRAL_REF = repo_path(config["viral_ref"])
BACT16S_REF = repo_path(config["bact16s_ref"])
BWA_INDEX_EXTENSIONS = [".amb", ".ann", ".bwt", ".pac", ".sa"]


def bwa_index_files(reference):
    return [f"{reference}{extension}" for extension in BWA_INDEX_EXTENSIONS]

# mothur references
MOTHUR_REFERENCE = repo_path(config["mothur_reference"])
MOTHUR_TAX_FILE = repo_path(config["mothur_tax_file"])
# mothur params
MOTHUR_CUTOFF = config["mothur_cutoff"]
MOTHUR_METHOD = config["mothur_method"]

# Kraken Parameters
KRAKEN_DB = config["kraken_db"]


# Trimming Parameters
TRIMMOMATIC_JAR = config["trimmomatic_jar"]
TRIM_ARGS = config["trim_args"]
TRIM_ARGS_STR =  " ".join(TRIM_ARGS)

# Dada2 Parameters
CHUNK_SIZE = config["chunk_size"]
TRUNC_LEN = config["truncLen"]
MAX_N = config["maxN"]
MAX_EE = config["maxEE"]
TRUNC_Q = config["truncQ"]

# BLAST Parameters
CONF_CUTOFF = config["conf_cutoff"]


# Phyloseq Parameters
ADD_UNCLASSIFIED_PREFIX = config["add_unclassified_prefix"]

#------------------------------------
# Helper Functions
#------------------------------------

def dada_input_reads():
    if PROCESS_UMIS:
        return UMI_DEDUP_READS
    return NORM_R1_READS


def sample_prep_done():
    if PROCESS_UMIS:
        return UMI_DEDUP_DONE
    return NO_UMI_COUNT_DONE


def count_summary_done():
    if PROCESS_UMIS:
        return UMI_SELECTION_DONE
    return NO_UMI_COUNT_DONE

#----------------------------
# All Rule
#----------------------------
rule all:
    input:
        phyloseq_image = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_physeq.RData",
        asv_fasta = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_ASV.fasta"

rule write_effective_config:
    output:
        effective_config = f"{OUT_DIR}/effective_config.yaml"
    run:
        Path(output.effective_config).parent.mkdir(parents=True, exist_ok=True)
        # JSON is valid YAML 1.2 and records Snakemake's fully merged config state.
        Path(output.effective_config).write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n"
        )



#----------------------------
# Validate FASTQ Manifest
#----------------------------
FASTQ_MANIFEST_RESULT = build_fastq_manifest(RAW)
SAMPLES = [row["SampleID"] for row in FASTQ_MANIFEST_RESULT.rows]
SAMPLE_NAMES = f"{OUT_DIR}/sample.names"
FASTQ_MANIFEST = f"{TRACK_DIR}/fastq_manifest.tsv"
FASTQ_VALIDATION_OK = f"{TRACK_DIR}/fastq_validation.ok"
FASTQ_VALIDATION_REPORT = f"{TRACK_DIR}/fastq_validation_report.txt"

NORM_R1_READS = expand(f"{NORM_RAW_DIR}/{{s}}_R1_001.fastq", s=SAMPLES)
NORM_R2_READS = expand(f"{NORM_RAW_DIR}/{{s}}_R2_001.fastq", s=SAMPLES)
NORM_FASTQ_DONE = f"{NORM_RAW_DIR}/.norm_fastq.done"

UMI_SELECTED_R1S = expand(f"{UMI_SELECT_DIR}/Selected.{{s}}.UMI_R1.fastq", s=SAMPLES)
UMI_COUNT_SUMMARIES = expand(f"{UMI_SELECT_DIR}/CountSummary.{{s}}.tsv", s=SAMPLES)
UMI_SELECTION_DONE = f"{UMI_SELECT_DIR}/.umi_selection.done"
UMI_DEDUP_READS = expand(f"{UMI_DEDUP_DIR}/Deduped.{{s}}.fastq", s=SAMPLES)
UMI_DEDUP_DONE = f"{UMI_DEDUP_DIR}/.umi_dedup.done"
NO_UMI_COUNT_DONE = f"{UMI_SELECT_DIR}/.no_umi_count_summary.done"


rule validate_fastqs:
    output:
        sample_names = SAMPLE_NAMES,
        manifest = FASTQ_MANIFEST,
        validation_ok = FASTQ_VALIDATION_OK
    log:
        report = FASTQ_VALIDATION_REPORT
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.sample_names}")" "$(dirname "{output.manifest}")" "$(dirname "{log.report}")"
        python3 "{SCRIPTS}/SamplePrep.py" validate \
            --raw-dir "{RAW}" \
            --manifest "{output.manifest}" \
            --sample-names "{output.sample_names}" \
            --report "{log.report}" \
            --ok "{output.validation_ok}"
        """

#-------------------------------------
# 00 Normalizing and Unzipping FASTQs
#-------------------------------------
rule norm_fastq:
    input:
        manifest = FASTQ_MANIFEST,
        validation_ok = FASTQ_VALIDATION_OK,
    output:
        r1_norm = NORM_R1_READS,
        r2_norm = NORM_R2_READS,
        done = NORM_FASTQ_DONE
    threads: 4
    log: f"{LOG_DIR}/00_norm/00_norm.log"
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "{NORM_RAW_DIR}"
        exec > "{log}" 2>&1
        python3 "{SCRIPTS}/SamplePrep.py" norm-fastq \
            --manifest "{input.manifest}" \
            --out-dir "{NORM_RAW_DIR}" \
            --threads {threads} \
            --done "{output.done}"
        """

#-----------------------------------------------
# 01 UMI extraction from R2 and Select R1 Reads
#-----------------------------------------------
if PROCESS_UMIS:
    rule umi_selection:
        input:
            manifest = FASTQ_MANIFEST,
            norm_done = NORM_FASTQ_DONE,
        output:
            sel_umi_r1 = UMI_SELECTED_R1S,
            count_summaries = UMI_COUNT_SUMMARIES,
            done = UMI_SELECTION_DONE
        params:
            r2_primer_motif = R2_PRIMER[:R2_PRIMER_MOTIF_LEN],
            r2_primer_skip_flag = "--r2-primer-skip" if R2_PRIMER_SKIP else "",
            poly_G_threshold = POLY_G_THRESHOLD,
            sample_log_dir = f"{LOG_DIR}/01_umi_select"
        threads: 4
        log: f"{LOG_DIR}/01_umi_select/01_umi_select.log"
        conda: conda_env("bio-tools-env")
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname "{log}")" "{UMI_SELECT_DIR}" "{params.sample_log_dir}"
            exec > "{log}" 2>&1
            python3 "{SCRIPTS}/SamplePrep.py" umi-selection \
                --manifest "{input.manifest}" \
                --norm-dir "{NORM_RAW_DIR}" \
                --out-dir "{UMI_SELECT_DIR}" \
                --sample-log-dir "{params.sample_log_dir}" \
                --umi-selection-script "{SCRIPTS}/UMISelection.py" \
                --r2-primer-motif "{params.r2_primer_motif}" \
                {params.r2_primer_skip_flag} \
                --poly-G-threshold {params.poly_G_threshold} \
                --umi-len "{UMI_LEN}" \
                --max-offset "{MAX_OFFSET}" \
                --threads {threads} \
                --done "{output.done}"
            """


    #----------------------------
    # 02 UMI Deduplication
    #----------------------------
    rule umi_dedup:
        input:
            manifest = FASTQ_MANIFEST,
            umi_selection_done = UMI_SELECTION_DONE,
        output:
            umi_dedup_reads = UMI_DEDUP_READS,
            done = UMI_DEDUP_DONE
        params:
            AmpUMI_regex = "^" + ("I" * UMI_LEN),
            sample_log_dir = f"{LOG_DIR}/02_umi_dedup"
        threads: 8
        log:    f"{LOG_DIR}/02_umi_dedup/02_umi_dedup.log"
        conda: conda_env("AmpUMI-env")
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname "{log}")" "{UMI_DEDUP_DIR}" "{params.sample_log_dir}"
            exec > "{log}" 2>&1
            python3 "{SCRIPTS}/SamplePrep.py" umi-dedup \
                --manifest "{input.manifest}" \
                --selected-dir "{UMI_SELECT_DIR}" \
                --out-dir "{UMI_DEDUP_DIR}" \
                --sample-log-dir "{params.sample_log_dir}" \
                --umi-regex "{params.AmpUMI_regex}" \
                --threads {threads} \
                --done "{output.done}"
            """
else:
    rule no_umi_count_summary:
        input:
            manifest = FASTQ_MANIFEST,
            norm_done = NORM_FASTQ_DONE,
        output:
            count_summaries = UMI_COUNT_SUMMARIES,
            done = NO_UMI_COUNT_DONE
        threads: 4
        log: f"{LOG_DIR}/01_no_umi_count_summary/01_no_umi_count_summary.log"
        conda: conda_env("bio-tools-env")
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname "{log}")" "{UMI_SELECT_DIR}"
            exec > "{log}" 2>&1
            python3 "{SCRIPTS}/SamplePrep.py" no-umi-count-summary \
                --manifest "{input.manifest}" \
                --norm-dir "{NORM_RAW_DIR}" \
                --out-dir "{UMI_SELECT_DIR}" \
                --threads {threads} \
                --done "{output.done}"
            """

#----------------------------
# 03 Dada Denoising
#---------------------------- 
rule dada_denoising:
    input: 
        sample_names = SAMPLE_NAMES,
        sample_prep_done = sample_prep_done()
    output:
        filtered_reads = temp(expand(f"{DADA_DENOISE_DIR}/filteredAndTrimmed/filtered.{{s}}.fastq", s=SAMPLES)),
        seq_err_plot = f"{DADA_DENOISE_DIR}/dada_error_plot.png",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        asv_map = f"{DADA_DENOISE_DIR}/ASVMap.tsv",
        filter_stage_counts = f"{DADA_DENOISE_DIR}/dada_read_counts.tsv",
    params:
        # Processign Params
        chunk_size = CHUNK_SIZE,
        # Filter and Trim Params
        primerLen = R1_PRIMER_MOTIF_LEN,
        truncLen = TRUNC_LEN,
        maxN = MAX_N,
        maxEE = MAX_EE,
        truncQ = TRUNC_Q,
        filtered_reads_dir = f"{DADA_DENOISE_DIR}/filteredAndTrimmed",
        dada_reads = dada_input_reads(),
    threads: 16
    log:    f"{LOG_DIR}/03_dada.log"
    conda: conda_env("R-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.seq_table}")" "{params.filtered_reads_dir}"
        exec > "{log}" 2>&1 

        cd "{REPO_ROOT}"
        Rscript "{SCRIPTS}/DadaASVFilter.R" \
            --fqs {params.dada_reads} \
            --sample-names {input.sample_names} \
            --filtered-fqs {output.filtered_reads} \
            --err-plt {output.seq_err_plot} \
            --filt-counts {output.filter_stage_counts} \
            --asv-fa {output.rep_asv_fasta} \
            --seq-table {output.seq_table} \
            --asv-map {output.asv_map} \
            --chunk-size {params.chunk_size} \
            --truncLen {params.truncLen} \
            --primerLen {params.primerLen} \
            --maxN {params.maxN} \
            --maxEE {params.maxEE} \
            --truncQ {params.truncQ} \
            --threads {threads}
        """

#--------------------------------------------------------
# 04 Host and Viral Mapping Selection
#--------------------------------------------------------
rule host_viral_alignment:
    input:
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        reference_fasta = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag],
        reference_indexes = lambda wc: bwa_index_files({"host": HOST_REF, "viral": VIRAL_REF}[wc.tag])
    output:
        sam = f"{NEG_ALIGNMENT_DIR}/{{tag}}.ASV.sam",
        unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.{{tag}}.ASV.names"
    threads: 8
    log:    f"{LOG_DIR}/04_alignment/04_alignment_{{tag}}.log"
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.sam}")" "$(dirname "{output.unmapped_names}")"
        bwa mem -t {threads} "{input.reference_fasta}" "{input.rep_asv_fasta}" > "{output.sam}" 2> "{log}"
        samtools view -f 4 "{output.sam}" 2> "{log}" | cut -f1 | sort -u > "{output.unmapped_names}"
        """

#-----------------------------------------
# 05 Unmapped Host and Viral Read Filter
#----------------------------------------- 
rule nonhost_nonviral_filter:
    input:
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
        viral_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.viral.ASV.names"
    output:
        nonhost_nonviral_ASVs = f"{NEG_ALIGN_FILT_DIR}/nonhost.nonviral.ASV.fasta"
    threads: 8
    log:    f"{LOG_DIR}/05_nonhost_nonviral.log"
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.nonhost_nonviral_ASVs}")"

        seqtk subseq "{input.rep_asv_fasta}" "{input.host_unmapped_names}" \
         | seqtk subseq - "{input.viral_unmapped_names}" \
         > "{output.nonhost_nonviral_ASVs}"
        """


#-------------------------------------
# 06 Positive Bacterial Alignment
#------------------------------------- 
rule bacterial_alignment:
    input:
        nonhost_nonviral_ASVs = f"{NEG_ALIGN_FILT_DIR}/nonhost.nonviral.ASV.fasta",
        reference_fasta = BACT16S_REF,
        reference_indexes = bwa_index_files(BACT16S_REF)
    output:
        bacterial_alignment = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.sam",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        bacterial_ASVs = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta"
    threads: 8
    log:    f"{LOG_DIR}/06_bact_map.log"
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.bacterial_alignment}")" "$(dirname "{output.bacterial_names}")" "$(dirname "{output.bacterial_ASVs}")"
        bwa mem -t {threads} "{input.reference_fasta}" "{input.nonhost_nonviral_ASVs}" > "{output.bacterial_alignment}" 2> "{log}"
        samtools view -F 4 "{output.bacterial_alignment}" | cut -f1 | sort -u > "{output.bacterial_names}"
        seqtk subseq "{input.nonhost_nonviral_ASVs}" "{output.bacterial_names}" > "{output.bacterial_ASVs}" 2> "{log}"
        """



#-------------------------------------
# 07 Canonical Pre-Decontamination ASV FASTA
#-------------------------------------
# copies the intermediate ASV to output directory to use as pipeline endpoint (used in batch compilation)
rule canonical_asv_fasta:
    input:
        bacterial_ASV_fa = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta"
    output:
        asv_fasta = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_ASV.fasta"
    log: f"{LOG_DIR}/07_canonical_asv_fasta.log"
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.asv_fasta}")"
        exec > "{log}" 2>&1
        cp "{input.bacterial_ASV_fa}" "{output.asv_fasta}"
        """

#--------------------------------------
# 08.2 Kraken2 Taxanomic Classifcation 
#--------------------------------------
rule kraken_classification:
    input:
        bacterial_ASV_fa = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta",
        kraken_database = f"{REF_DIR}/{KRAKEN_DB}"
    output:
        kraken_class_file = f"{KRAKEN_TAX_DIR}/bacterial.ASV.{KRAKEN_DB}.kraken2",
        kraken_report_file = f"{KRAKEN_TAX_DIR}/bacterial.ASV.{KRAKEN_DB}.k2report",
    threads: 16
    log:    f"{LOG_DIR}/08.2_kraken_class.log"
    conda: conda_env("kraken-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.kraken_class_file}")" "$(dirname "{output.kraken_report_file}")"
        exec > "{log}" 2>&1
        kraken2 --threads {threads} \
                --db {input.kraken_database} \
                --report {output.kraken_report_file} \
                --output {output.kraken_class_file} \
                {input.bacterial_ASV_fa}
        """

#--------------------------------------
# 00 BLAST Database Construction
#--------------------------------------
# Not used by pipeline but used later in confirmatory taxanomic analysis
rule blast_db_construction:
    input:
        reference_db_fasta = f"{REF_DIR}/{KRAKEN_DB}.fasta",
        taxid_map = f"{REF_DIR}/{KRAKEN_DB}/seqid2taxid.map"
    output:
        #representative database construction artifact inidcating successful database construction
        blast_db = f"{REF_DIR}/{KRAKEN_DB}/blast_format_db/{KRAKEN_DB}_blast.nsq"
    params: 
        blast_db_base = f"{REF_DIR}/{KRAKEN_DB}/blast_format_db/{KRAKEN_DB}_blast"
    threads: 16
    log:    f"{LOG_DIR}/00_blast_db_construction.log"
    conda: conda_env("bio-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.blast_db}")"
        exec > "{log}" 2>&1

        makeblastdb -in {input.reference_db_fasta} \
            -dbtype nucl \
            -parse_seqids \
            -taxid_map {input.taxid_map} \
            -out {params.blast_db_base}
        """

rule phyloseq_construction:
    input:
        effective_config = f"{OUT_DIR}/effective_config.yaml",
        kraken_class_file = f"{KRAKEN_TAX_DIR}/bacterial.ASV.{KRAKEN_DB}.kraken2",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        raw_seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        sample_names = SAMPLE_NAMES,
        metadata_sheet = METADATA,
        library_counts = f"{TRACK_DIR}/combined_read_counts.tsv"
    output:
        phyloseq_object = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_physeq.RData",
        missing_metadata_report = f"{PHYLOSEQ_DIR}/DroppedSamplesMissingMetadata.tsv"
    params:
        dump_dir = f"{REF_DIR}/taxdump",        
        unclassified_prefix_flag = "--add-unclassified-prefix" if ADD_UNCLASSIFIED_PREFIX else ""
    threads: 8
    log:    f"{LOG_DIR}/10_phyloseq.log"
    conda: conda_env("R-tools-env")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.phyloseq_object}")" "$(dirname "{output.missing_metadata_report}")"
        exec > "{log}" 2>&1
        cd "{REPO_ROOT}"
        Rscript "{SCRIPTS}/PhyloseqConstruction.R" \
            --kraken-file {input.kraken_class_file} \
            --bacterial-names {input.bacterial_names} \
            --raw-seq-table {input.raw_seq_table} \
            --sample-names {input.sample_names} \
            --metadata {input.metadata_sheet} \
            --library-counts {input.library_counts} \
            --dump-dir {params.dump_dir} \
            {params.unclassified_prefix_flag} \
            --trialID {TRIAL_ID} \
            --out {PHYLOSEQ_DIR}
        """


#-----------------------------
# 11 Read Count Calculations
#-----------------------------
rule read_counts:
    input:
        sample_names = SAMPLE_NAMES,
        count_summary_done = count_summary_done(),
        dada_read_counts = f"{DADA_DENOISE_DIR}/dada_read_counts.tsv",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
        viral_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.viral.ASV.names",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
    output:
        library_counts = f'{TRACK_DIR}/combined_read_counts.tsv'
    params:
        selected = UMI_SELECT_DIR,
    threads: 8
    log:    f"{LOG_DIR}/09_read_count.log"
    conda: conda_env("R-tools-env")

    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.library_counts}")"
        exec > "{log}" 2>&1

        cd "{REPO_ROOT}"
        Rscript "{SCRIPTS}/ReadCountCompilation.R" \
            --sample-name-file {input.sample_names} \
            --selected-dir {params.selected} \
            --dada-filter-counts {input.dada_read_counts} \
            --seq-table {input.seq_table} \
            --host-names {input.host_unmapped_names} \
            --viral-names {input.viral_unmapped_names} \
            --bacterial-names {input.bacterial_names} \
            --combined-counts {output.library_counts}
        """
