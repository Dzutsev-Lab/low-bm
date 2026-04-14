import os, re, glob
from pathlib import Path

configfile: "config.yaml"

# Input and Output Directories
TRIAL_ID = config["trialID"]
TRIAL_NAME = TRIAL_ID + "_" + config["trial_descript"]
EXP_NAME = config["exp_name"]

IN_DIR = os.path.join(config["in_root"], f"{EXP_NAME}_Data")
IP_DIR = os.path.join(config["ip_root"], EXP_NAME)
OUT_DIR = os.path.join(config["out_root"], TRIAL_NAME)

SCRIPTS = config["script_dir"]
REF_DIR = "Ref_Data"
METADATA = f"{IN_DIR}/{EXP_NAME}_metadata.xlsx"
CONDA_ENV_DIR = "/vf/users/taylorng/conda/envs"


RAW = f"{IN_DIR}/OUTPUT/Data/fastq"
NORM_RAW_DIR = f"{IP_DIR}/00_RawNorm"
UMI_SELECT_DIR = f"{IP_DIR}/01_UMISelection"
UMI_DEDUP_DIR = f"{IP_DIR}/02_UMIDeduplication"
DADA_DENOISE_DIR = f"{IP_DIR}/03_Dada_Denoising"
NEG_ALIGNMENT_DIR = f"{IP_DIR}/04_Neg_Alignment"
NEG_ALIGN_FILT_DIR = f"{IP_DIR}/05_Neg_Alignment_Filter"
POS_ALIGNMENT_DIR = f"{IP_DIR}/06_Pos_Alignment"
MICROCLEAN_DECONTAM_DIR = f"{IP_DIR}/07_micRoclean_Decontam"
MOTHUR_TAX_DIR = f"{IP_DIR}/08.1_Mothur_Taxonomy"
KRAKEN_TAX_DIR = f"{IP_DIR}/08.2_Kraken_Taxonomy"
NORM_COUNT_DIR = f"{IP_DIR}/09_CountNormalization"
PHYLOSEQ_DIR = f"{OUT_DIR}"

TRACK_DIR = f"{OUT_DIR}/Tracking"
LOG_DIR = f"{OUT_DIR}/Logs"

# Sequencing metadata
R1_PRIMER = config["r1_primer"]
R2_PRIMER = config["r2_primer"]
R1_PRIMER_MOTIF_LEN = config["r1_primer_motif_len"]
R2_PRIMER_MOTIF_LEN = config["r2_primer_motif_len"]
R2_PRIMER_SKIP = config["r2_primer_skip"]
POLY_G_THRESHOLD = config["poly_G_threshold"]
UMI_LEN = config["umi_len"]
MAX_OFFSET = config["max_offset"]

# Reference Data
HOST = config["host"]
if HOST == "human":
    HOST_REF = config["human_ref"]
elif HOST == "mouse":
    HOST_REF = config["mouse_ref"]

VIRAL_REF = config["viral_ref"]
BACT16S_REF = config["bact16s_ref"]

# mothur references
MOTHUR_REFERENCE = config["mothur_reference"]
MOTHUR_TAX_FILE = config["mothur_tax_file"]
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

# Read Count Normalization Parameters
NORM_METHOD = config["norm_method"]
NORM_OFFSET = config["norm_offset"]

# Phyloseq Analysis Parameters
ADD_UNCLASSIFIED_PREFIX = config["add_unclassified_prefix"]

# Negative Control Filtering Parameters
# TODO: figure out how to format for snakemake without config file
NEG_DB_TOP_K="10000"
NEG_AS_KEEP_LT="160"   

# Taxonomy Cutoffs
TAX_MIN = config["tax_min_counts"]
TAX_LEVELS = ['genus', 'species', 'class', 'phylum']
TAX_FIELDS = {
    "phylum":  "2",
    "class":   "2,3",
    "genus":   "2,3,4,5,6",
    "species": "2,3,4,5,6,7",
}


#------------------------------------
# Helper Functions
#------------------------------------

# Sample Discovery
def discover_samples():
    pats = [
        os.path.join(RAW, "*_R1_001.fastq"),
        os.path.join(RAW, "*_R1_001.fastq.gz"),
    ]
    r1s = []
    for p in pats:
        r1s.extend(glob.glob(p))
    samples = []
    for r1 in sorted(r1s):
        name = os.path.basename(r1)
        base = name.replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
        # skip undertermined fastq
        if "Undetermined" in base:
            continue
        # require matching R2 (either gz or not)
        r2a = os.path.join(RAW, f"{base}_R2_001.fastq")
        r2b = os.path.join(RAW, f"{base}_R2_001.fastq.gz")
        if os.path.exists(r2a) or os.path.exists(r2b):
            samples.append(base)
    if not samples:
        raise ValueError(f"No samples found in {RAW} matching *_R1_001.fastq(.gz) with paired R2.")
    return samples


# Input FASTQ Compression check
#   checks if the input raw fastq need to be unzipped for input to raw normalization
def pick_raw_fastq(wc, read):
    fastq_path = os.path.join(RAW, f"{wc.s}_R{read}_001.fastq")
    gz_path = fastq_path + ".gz"
    if os.path.exists(fastq_path):
        return fastq_path
    elif os.path.exists(gz_path):
        return gz_path
    else:
        raise (ValueError(f"Missing raw FASTQ for sample {wc.s} read R{read}: {fastq_path}(.gz)"))

#----------------------------
# All Rule
#----------------------------
rule all:
    input:
        kraken_countXtypeXsample_plot = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_kraken_countXtypeXsample.png",
        kraken_genus_table = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_kraken_genus_table.tsv"


rule copy_config:
    output:
        config_copy = f"{OUT_DIR}/config.yaml"
    shell:
        r"""
        set -euo pipefail
        cp config.yaml {output.config_copy}
        """



#----------------------------
# Record Sample Names
#----------------------------
#Construct global variable of all samples based on available fastq's
SAMPLES = discover_samples()
#makes samples list into an format that works better with bash for-loop in read counts rule
SAMPLES_STR = " ".join(SAMPLES)
rule sample_names:
    output:
        f"{OUT_DIR}/sample.names"
    run:
        Path(output[0]).write_text("\n".join(SAMPLES) + "\n")

#----------------------------
# Reference Indexing
#----------------------------  
rule index_ref_bwa:
    input:
        reference = "{ref}"
    output:
        # indexing produces many files but .bwt is the key one (will use as sentinel)
        bwt = "{ref}.bwt"
    threads: 8
    log: f"{LOG_DIR}/00_bwa_{{ref}}_ref.log"
    conda: f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1
        bwa index {input.reference}
        """

#------------------------------
# Kraken Database Construction
#------------------------------
rule kraken_db_construction:
    output:
        # using database dump files as sentinal 
        # (needed for tax table reconstruction in phyloseq analysis)
        names_dump = f"{REF_DIR}/{KRAKEN_DB}/taxonomy/names.dmp",
        nodes_dump = f"{REF_DIR}/{KRAKEN_DB}/taxonomy/nodes.dmp"
    params:
        kraken_db = KRAKEN_DB
    threads: 16
    log: f"{LOG_DIR}/00_kraken_db.log"
    conda: f"{CONDA_ENV_DIR}/kraken-env"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1
        
        kraken2-build \
            --db {REF_DIR}/{params.kraken_db} \
            --special {params.kraken_db} \
            --threads {threads}
        """


#-------------------------------------
# 00 Normalizing and Unzipping FASTQs
#-------------------------------------
rule norm_fastq:
    input:
        r1_raw = lambda wc: pick_raw_fastq(wc, 1),
        r2_raw = lambda wc: pick_raw_fastq(wc, 2),
    output:
        r1_norm = temp(f"{NORM_RAW_DIR}/{{s}}_R1_001.fastq"),
        r2_norm = temp(f"{NORM_RAW_DIR}/{{s}}_R2_001.fastq")
    threads: 2
    log: f"{LOG_DIR}/00_norm/00_norm.{{s}}.log"
    conda: f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
        exec 2> "{log}"

        # check if input is compressed or not
        # decompress to norm raw directory if it is
        # add link to original raw fastq in norm raw directory
        #   if it is already decompress
        norm_fastq() {{
            local fastq_in="$1"
            local fastq_out="$2"
            if [[ "$fastq_in" == *.gz ]]; then
                pigz -dc "$fastq_in" > "$fastq_out"
            else
                ln -sf "$(realpath "$fastq_in")" "$fastq_out"
            fi
        }}

        norm_fastq "{input.r1_raw}" "{output.r1_norm}"
        norm_fastq "{input.r2_raw}" "{output.r2_norm}"
        """

#-----------------------------------------------
# 01 UMI extraction from R2 and Select R1 Reads
#-----------------------------------------------
rule umi_selection:
    input:
        r1 = f'{NORM_RAW_DIR}/{{s}}_R1_001.fastq',
        r2 = f'{NORM_RAW_DIR}/{{s}}_R2_001.fastq',
    output:
        sel_umi_r1 = temp(f"{UMI_SELECT_DIR}/Selected.{{s}}.UMI_R1.fastq"),
        count_summary = f"{UMI_SELECT_DIR}/CountSummary.{{s}}.tsv"
    params:
        r2_primer_motif = R2_PRIMER[:R2_PRIMER_MOTIF_LEN],
        r2_primer_skip_flag = "--r2-primer-skip" if R2_PRIMER_SKIP else "",
        poly_G_threshold = POLY_G_THRESHOLD
    threads: 1
    log: f"{LOG_DIR}/01_umi_select/01_umi_select.{{s}}.log"
    conda: f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
        exec 2> "{log}"

        python3 scripts/UMISelection.py \
            --sample-name "{wildcards.s}" \
            --r1 "{input.r1}" \
            --r2 "{input.r2}" \
            --r2-primer-motif "{params.r2_primer_motif}" \
            {params.r2_primer_skip_flag} \
            --poly-G-threshold {params.poly_G_threshold} \
            --umi-len "{UMI_LEN}" \
            --max-offset "{MAX_OFFSET}" \
            --out-count-summary "{output.count_summary}" \
            --out-umi-r1 "{output.sel_umi_r1}"
        """


#----------------------------
# 02 UMI Deduplication
#----------------------------
rule umi_dedup:
    input:
        selected_umi_r1 = f"{UMI_SELECT_DIR}/Selected.{{s}}.UMI_R1.fastq",
    output:
        umi_dedup_reads = temp(f"{UMI_DEDUP_DIR}/Deduped.{{s}}.fastq"),
    params:
        AmpUMI_regex = "^" + ("I" * UMI_LEN)
    threads: 8
    log:    f"{LOG_DIR}/02_umi_dedup/02_umi_dedup.{{s}}.log"
    conda:  f"{CONDA_ENV_DIR}/AmpUMI-env"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1

        # ---- Check if Empty FASTQ ----
        n_lines=$(wc -l < {input.selected_umi_r1})

        if [[ "$n_lines" -eq 0 ]]; then
            echo "Input FASTQ ({input.selected_umi_r1}) is empty - skipping AmpUMI"
            : > {output.umi_dedup_reads}
        else
            echo "Running AmpUMI on $n_lines lines from {input.selected_umi_r1}"
            AmpUMI Process\
                --fastq {input.selected_umi_r1} \
                --fastq_out {output.umi_dedup_reads} \
                --umi_regex "{params.AmpUMI_regex}"
        fi
        """

#----------------------------
# 03 Dada Denoising
#---------------------------- 
rule dada_denoising:
    input: 
        sample_names = f"{OUT_DIR}/sample.names",
        umi_dedup_reads = expand(f"{UMI_DEDUP_DIR}/Deduped.{{s}}.fastq", s=SAMPLES)
    output:
        filtered_reads = temp(expand(f"{DADA_DENOISE_DIR}/filteredAndTrimmed/filtered.{{s}}.fastq", s=SAMPLES)),
        seq_err_plot = f"{DADA_DENOISE_DIR}/dada_error_plot.png",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
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
    threads: 16
    log:    f"{LOG_DIR}/03_dada.log"
    conda:  f"{CONDA_ENV_DIR}/R-tools-env"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1 

        Rscript scripts/DadaASVFilter.R \
            --fqs {input.umi_dedup_reads} \
            --sample-names {input.sample_names} \
            --filtered-fqs {output.filtered_reads} \
            --err-plt {output.seq_err_plot} \
            --filt-counts {output.filter_stage_counts} \
            --asv-fa {output.rep_asv_fasta} \
            --seq-table {output.seq_table} \
            --chunk-size {params.chunk_size} \
            --truncLen {params.truncLen} \
            --primerLen {params.primerLen} \
            --maxN {params.maxN} \
            --maxEE {params.maxEE} \
            --truncQ {params.truncQ} 
        """

#--------------------------------------------------------
# 04 Host and Viral Mapping Selection
#--------------------------------------------------------
rule host_viral_alignment:
    input:
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        reference_fasta = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag],
        reference_bwt = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag] + ".bwt"
    output:
        sam = f"{NEG_ALIGNMENT_DIR}/{{tag}}.ASV.sam",
        unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.{{tag}}.ASV.names"
    threads: 8
    log:    f"{LOG_DIR}/04_alignment/04_alignment_{{tag}}.log"
    conda:  f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
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
    conda:  f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail

        seqtk subseq "{input.rep_asv_fasta}" "{input.host_unmapped_names}" \
         | seqtk subseq - "{input.viral_unmapped_names}" \
         > "{output.nonhost_nonviral_ASVs}"
        """


#-------------------------------------
# 06 Positive Bacterial Alignment
#------------------------------------- 
rule bacterial_alignment:
    input:
        nonhost_nonviral_ASVs = f"{NEG_ALIGN_FILT_DIR}/nonhost.nonviral.ASV.fasta"
    output:
        bacterial_alignment = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.sam",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        bacterial_ASVs = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta"
    threads: 8
    log:    f"{LOG_DIR}/06_bact_map.log"
    conda:  f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} "{BACT16S_REF}" "{input.nonhost_nonviral_ASVs}" > "{output.bacterial_alignment}" 2> "{log}"
        samtools view -F 4 "{output.bacterial_alignment}" | cut -f1 | sort -u > "{output.bacterial_names}"
        seqtk subseq "{input.nonhost_nonviral_ASVs}" "{output.bacterial_names}" > "{output.bacterial_ASVs}" 2> "{log}"
        """



#-------------------------------------
# 07 micRoclean Decontamination
#------------------------------------- 
rule micRoclean_decontamination_detection:
    input:
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        metadata_sheet = METADATA
    output:
        decontaminated_seq_table = f"{MICROCLEAN_DECONTAM_DIR}/DecontamSeqTable.tsv",
        decontaminated_names = f"{MICROCLEAN_DECONTAM_DIR}/decontaminated.ASV.names",
        filtering_report = f"{MICROCLEAN_DECONTAM_DIR}/DecontamFilterReport.tsv",
        
    threads: 8
    log:    f"{LOG_DIR}/07.1_micRoclean_decontam.log"
    conda:  f"{CONDA_ENV_DIR}/micRoclean-env-new"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        Rscript scripts/Decontamination.R \
            --seq-table {input.seq_table} \
            --bacterial-names {input.bacterial_names} \
            --metadata {input.metadata_sheet} \
            --trialID {TRIAL_ID} \
            --out {MICROCLEAN_DECONTAM_DIR}
        """

rule decontamination_filter:
    input:
        bacterial_ASV_fa = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta",
        decontaminated_names = f"{MICROCLEAN_DECONTAM_DIR}/decontaminated.ASV.names"
    output:
        decontaminated_ASV_fa = f"{MICROCLEAN_DECONTAM_DIR}/decontaminated.ASV.fasta"
    threads: 8
    log:    f"{LOG_DIR}/07.2_decontam_filter.log"
    conda:  f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        seqtk subseq "{input.bacterial_ASV_fa}" "{input.decontaminated_names}" > "{output.decontaminated_ASV_fa}" 2> "{log}"
        """


#------------------------------
# 08.1 Mothur Taxa Classifcation 
#------------------------------
rule mothur_classify:
    input:
        bacterial_ASVs = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta",
        mothur_ref = MOTHUR_REFERENCE,
        mothur_tax_file = MOTHUR_TAX_FILE
    output:
        taxfile = f"{MOTHUR_TAX_DIR}/bacterial.ASV.ncbi20.wang.taxonomy",
        tax_summary = f"{MOTHUR_TAX_DIR}/bacterial.ASV.ncbi20.wang.tax.summary",
    threads: 8
    log:    f"{LOG_DIR}/08.1_mothur_class.log"
    conda:  f"{CONDA_ENV_DIR}/mothur-env"
    shell:
        r"""
        set -euo pipefail

        # ---- Run Mothur ----
        mothur "#set.dir(output={MOTHUR_TAX_DIR}); 
                set.logfile(name={log});
                classify.seqs(
                    fasta={input.bacterial_ASVs}, 
                    reference={input.mothur_ref}, 
                    taxonomy={input.mothur_tax_file}, 
                    method={MOTHUR_METHOD}, 
                    cutoff={MOTHUR_CUTOFF}, 
                    processors={threads}
                )" >> {log} 2>&1
        
        
        # ---- Remove Redundant Mothur Logs ----
        rm *.logfile
        """

#--------------------------------------
# 08.2 Kraken2 Taxanomic Classifcation 
#--------------------------------------
rule kraken_classification:
    input:
        decontaminated_ASV_fa = f"{MICROCLEAN_DECONTAM_DIR}/decontaminated.ASV.fasta",
        kraken_database = f"{REF_DIR}/{KRAKEN_DB}"
    output:
        kraken_class_file = f"{KRAKEN_TAX_DIR}/bacterial.ASV.{KRAKEN_DB}.kraken2",
        kraken_report_file = f"{KRAKEN_TAX_DIR}/bacterial.ASV.{KRAKEN_DB}.k2report",
    threads: 16
    log:    f"{LOG_DIR}/08.2_kraken_class.log"
    conda:  f"{CONDA_ENV_DIR}/kraken-env"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        kraken2 --threads {threads} \
                --db {input.kraken_database} \
                --report {output.kraken_report_file} \
                --output {output.kraken_class_file} \
                {input.decontaminated_ASV_fa}
        """


#--------------------------
# 09 Count Normalization
#--------------------------
rule count_normalization:
    input:
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
    output:
        norm_seq_table = f"{NORM_COUNT_DIR}/NormSeqTable.tsv"
    params:
        method = NORM_METHOD,
        offset = NORM_OFFSET
    threads: 1
    log:    f"{LOG_DIR}/09_normalization.log"
    conda:  f"{CONDA_ENV_DIR}/bio-tools-env"
    shell:
        r"""
        set -euo pipefail
        python scripts/SequenceTableNormalization.py \
            --in {input.seq_table} \
            --host-names {input.host_unmapped_names} \
            --out {output.norm_seq_table} \
            --method {params.method} \
            --offset {params.offset} \
        > "{log}" 2>&1
        """

#-------------------------------
# 10 Phyloseq Taxonomy Analysis
#-------------------------------
rule phyloseq_analysis:
    input:
        norm_seq_table = f"{NORM_COUNT_DIR}/NormSeqTable.tsv",
        raw_seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        kraken_file = f"{KRAKEN_TAX_DIR}/bacterial.ASV.{KRAKEN_DB}.kraken2",
        names_dump = f"{REF_DIR}/{KRAKEN_DB}/taxonomy/names.dmp",
        nodes_dump = f"{REF_DIR}/{KRAKEN_DB}/taxonomy/nodes.dmp",
        combined_read_counts = f'{TRACK_DIR}/combined_read_counts.tsv',
        metadata_sheet = METADATA
    output:
        kraken_abunXtypeXsample_plot = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_kraken_countXtypeXsample.png",
        kraken_genus_table = f"{PHYLOSEQ_DIR}/{TRIAL_ID}_kraken_genus_table.tsv"
    params:
        kraken_db = KRAKEN_DB,
        unclassified_prefix_flag = "--add-unclassified-prefix" if ADD_UNCLASSIFIED_PREFIX else "",
        norm_method = NORM_METHOD
    threads: 8
    log:    f"{LOG_DIR}/10_phyloseq.log"
    conda:  f"{CONDA_ENV_DIR}/R-tools-env"

    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        Rscript scripts/PhyloseqTaxAnalysis.R \
            --kraken-file {input.kraken_file} \
            --dump-dir {REF_DIR}/{params.kraken_db}/taxonomy \
            --norm-method {params.norm_method} \
            --norm-seq-table {input.norm_seq_table} \
            --raw-seq-table {input.raw_seq_table} \
            --bacterial-names {input.bacterial_names} \
            --trialID {TRIAL_ID} \
            {params.unclassified_prefix_flag} \
            --read-count-file {input.combined_read_counts} \
            --metadata {input.metadata_sheet} \
            --out {PHYLOSEQ_DIR}
        """

#-----------------------------
# -01 Read Count Calculations
#-----------------------------
rule read_counts:
    input:
        sample_names = f"{OUT_DIR}/sample.names",
        
        umi_selection_count_summary_tsv = expand(f"{UMI_SELECT_DIR}/CountSummary.{{s}}.tsv", s=SAMPLES),

        dada_read_counts = f"{DADA_DENOISE_DIR}/dada_read_counts.tsv",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
        viral_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.viral.ASV.names",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
    output:
        combined_read_counts = f'{TRACK_DIR}/combined_read_counts.tsv'
    params:
        selected = UMI_SELECT_DIR,
    threads: 8
    log:    f"{LOG_DIR}/-01_read_count.log"
    conda:  f"{CONDA_ENV_DIR}/R-tools-env"

    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1

        Rscript scripts/ReadCountCompilation.R \
            --sample-name-file {input.sample_names} \
            --selected-dir {params.selected} \
            --dada-filter-counts {input.dada_read_counts} \
            --seq-table {input.seq_table} \
            --host-names {input.host_unmapped_names} \
            --viral-names {input.viral_unmapped_names} \
            --bacterial-names {input.bacterial_names} \
            --combined-counts {output.combined_read_counts}
        """

        



