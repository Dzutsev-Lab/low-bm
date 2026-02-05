# TODO: adjust log file numbering
import os, re, glob
from pathlib import Path

configfile: "config.yaml"

# Input and Output Directories
RAW = config["raw_dir"]
OUT = config["out_root"]
SCRIPTS = config["script_dir"]

ORIG_DIR = f"{OUT}/original.fastq"

IP_DIR = f"{OUT}/InProcess"
NORM_RAW_DIR = f"{IP_DIR}/00_RawNorm"
UMI_SELECT_DIR = f"{IP_DIR}/01_UMISelection"
DADA_DENOISE_DIR = f"{IP_DIR}/02_Dada_Denoising"
NEG_ALIGNMENT_DIR = f"{IP_DIR}/03_Neg_Alignment"
NEG_ALIGN_FILT_DIR = f"{IP_DIR}/04_Neg_Alignment_Filter"
POS_ALIGNMENT_DIR = f"{IP_DIR}/05_Pos_Alignment"
ID_TAX_DIR = f"{IP_DIR}/06.1_IDTax_Taxonomy"
MOTHUR_TAX_DIR = f"{IP_DIR}/06.2-08_Mothur_Taxonomy"
TRIMMED_DIR = f"{IP_DIR}/-00_Trimming"
UMI_DEDUP_DIR = f"{IP_DIR}/UMIDeduplication"
NORM_COUNT_DIR = f"{IP_DIR}/CountNormalization"


TRACK_DIR = f"{OUT}/Tracking"
TAX_OUTPUT_DIR = f"{OUT}/Taxonomy"
LOG_DIR = f"{OUT}/logs"

# Sequencing metadata
R1_PRIMER = config["r1_primer"]
R2_PRIMER = config["r2_primer"]
R1_PRIMER_MOTIF_LEN = config["r1_primer_motif_len"]
R2_PRIMER_MOTIF_LEN = config["r2_primer_motif_len"]
UMI_LEN = config["umi_len"]
MAX_OFFSET = config["max_offset"]

# Reference Data
HOST_REF = config["host_ref"]
VIRAL_REF = config["viral_ref"]
BACT16S_REF = config["bact16s_ref"]

# mothur references
MOTHUR_REFERENCE = config["mothur_reference"]
MOTHUR_TAX_FILE = config["mothur_tax_file"]
# mothur params
MOTHUR_CUTOFF = config["mothur_cutoff"]
MOTHUR_METHOD = config["mothur_method"]

# Dada reference
ID_TAX_REFERENCE = config["IDTax_reference"]


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
        # require matching R2 (either gz or not)
        r2a = os.path.join(RAW, f"{base}_R2_001.fastq")
        r2b = os.path.join(RAW, f"{base}_R2_001.fastq.gz")
        if os.path.exists(r2a) or os.path.exists(r2b):
            samples.append(base)
    if not samples:
        raise ValueError(f"No samples found in {RAW} matching *_R1_001.fastq(.gz) with paired R2.")
    return samples

SAMPLES = discover_samples()
#makes samples list into an format that works better with bash for-loop in read counts rule
SAMPLES_STR = " ".join(SAMPLES)


#----------------------------
# All Rule
#----------------------------
rule all:
    input:
        abunXtype_plot = f"{TAX_OUTPUT_DIR}/abunXtype.png",
        abunXsample_plot = f"{TAX_OUTPUT_DIR}/abunXsample.png",
        abunXtypeXsample_plot = f"{TAX_OUTPUT_DIR}/abunXtypeXsample.png",
        tax_summary = f"{TAX_OUTPUT_DIR}/bacterial.ASV.ncbi20.wang.tax.summary",
        tax_count = expand(f"{TAX_OUTPUT_DIR}/bacterial.ASV.{{level}}.count", level=TAX_LEVELS),

        dada_read_counts = f"{TRACK_DIR}/dada_read_counts.tsv",
        combined_read_counts = f"{TRACK_DIR}/combined_read_counts.tsv",
        norm_seq_table = f"{NORM_COUNT_DIR}/NormSeqTable.tsv"




#----------------------------
# Record Sample Names
#----------------------------
rule sample_names:
    output:
        f"{OUT}/sample.names"
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
    shell:
        r"""
        set -euo pipefail
        bwa index "{input.ref}"
        """

#----------------------------
# Link Original Files
#----------------------------
# TODO: implement


#-------------------------------------
# 00 Normalizing and Unzipping FASTQs
#-------------------------------------
#checks if the input raw fastq is already decompressed or not for feeding input to raw normalization
def pick_raw_fastq(wc, read):
    fastq_path = os.path.join(RAW, f"{wc.s}_R{read}_001.fastq")
    gz_path = fastq_path + ".gz"
    if os.path.exists(fastq_path):
        return fastq_path
    elif os.path.exists(gz_path):
        return gz_path
    else:
        raise (ValueError(f"Missing raw FASTQ for sample {wc.s} read R{read}: {fastq_path}(.gz)"))


rule norm_fastq:
    input:
        r1_raw = lambda wc: pick_raw_fastq(wc, 1),
        r2_raw = lambda wc: pick_raw_fastq(wc, 2),
    output:
        r1_norm = f"{NORM_RAW_DIR}/{{s}}_R1_001.fastq",
        r2_norm = f"{NORM_RAW_DIR}/{{s}}_R2_001.fastq"
    log:
        f"{LOG_DIR}/00_norm/00_norm.{{s}}.log"
    threads: 2
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
        umi_tsv = f"{UMI_SELECT_DIR}/UMIMap.{{s}}.tsv",
        sel_names = f"{UMI_SELECT_DIR}/SelectedNames.{{s}}.names",
        sel_umi_r1 = f"{UMI_SELECT_DIR}/Selected.{{s}}.UMI_R1.fastq",
        sel_r1 = f"{UMI_SELECT_DIR}/Selected.{{s}}.R1.fastq",
    params:
        r2_primer_motif = R2_PRIMER[:R2_PRIMER_MOTIF_LEN]
    log:
        f"{LOG_DIR}/01_umi_select/01_umi_select.{{s}}.log"
    # TODO: address seqtk dependency (most likely using singularity)
    shell:
        r"""
        set -euo pipefail
        exec 2> "{log}"

        python3 scripts/UMISelection.py \
            --r1 "{input.r1}" \
            --r2 "{input.r2}" \
            --r2-primer-motif "{params.r2_primer_motif}" \
            --umi-len "{UMI_LEN}" \
            --max-offset "{MAX_OFFSET}" \
            --out-umi-tsv "{output.umi_tsv}" \
            --out-sel-names "{output.sel_names}" \
            --out-umi-r1 "{output.sel_umi_r1}" \
            --out-sel-r1 "{output.sel_r1}"
        """


#----------------------------
# UMI Deduplication
#----------------------------
rule umi_dedup:
    input:
        selected_umi_r1 = f"{UMI_SELECT_DIR}/Selected.{{s}}.UMI_R1.fastq",
    output:
        umi_dedup_reads = f"{UMI_DEDUP_DIR}/Deduped.{{s}}.fastq",
    params:
        AmpUMI_regex = "^" + ("I" * UMI_LEN)
    log:
        f"{LOG_DIR}/umi_dedup/umi_dedup.{{s}}.log"
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
# 02 Dada Denoising
#---------------------------- 
rule dada_denoising:
    input: 
        sample_names = f"{OUT}/sample.names",
        umi_dedup_reads = expand(f"{UMI_DEDUP_DIR}/Deduped.{{s}}.fastq", s=SAMPLES)
    output:
        filtered_reads = expand(f"{DADA_DENOISE_DIR}/filteredAndTrimmed/filtered.{{s}}.fastq", s=SAMPLES),
        seq_err_plot = f"{TRACK_DIR}/dada_error_plot.png",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        filter_stage_counts = f"{TRACK_DIR}/dada_read_counts.tsv",
        
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
    log:
        f"{LOG_DIR}/02_dada.log"
    script:
        f"{SCRIPTS}/DadaASVFilter.R"

#--------------------------
# Count Normalization
#--------------------------
rule count_normalization:
    #--------------------------
    #--------------------------
    # Normalization Reasoning
    #--------------------------
    #--------------------------
    # TSS Normalization  
    #   - for exploratory normalization starting with simplest TSS
    #   - first normalizing my row sum (seems to be the best measure of authentic reads after PCR correction)
    # TSS with Host Mapped Read Counts
    #   - normalizing whole sequence table with row sums from Human Mapped subsetted data frame
    input:
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
    output:
        norm_seq_table = f"{NORM_COUNT_DIR}/NormSeqTable.tsv"
    log:
        f"{LOG_DIR}/normalization.log"
    shell:
        r"""
        set -euo pipefail
        python scripts/SequenceTableNormalization.py \
            --in {input.seq_table} \
            --host-names {input.host_unmapped_names} \
            --out {output.norm_seq_table} \
        > "{log}" 2>&1
        """
    





#--------------------------------------------------------
# 03 Contaminant Alignment and Unmapped Selection
#--------------------------------------------------------
rule host_viral_alignment:
    input:
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        reference_fasta = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag],
        reference_bwt = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag] + ".bwt"
    output:
        sam = f"{NEG_ALIGNMENT_DIR}/{{tag}}.ASV.sam",
        unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.{{tag}}.ASV.names"
    log:
        f"{LOG_DIR}/03_alignment/03_alignment_{{tag}}.log"
    threads: 8
    # TODO: address BWA dependency (most likely using singularity)
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} "{input.reference_fasta}" "{input.rep_asv_fasta}" > "{output.sam}" 2> "{log}"
        samtools view -f 4 "{output.sam}" 2> "{log}" | cut -f1 | sort -u > "{output.unmapped_names}"
        """

#-------------------------------------
# 04 Unmapped Host and Viral Read Filter
#------------------------------------- 
rule nonhost_nonviral_filter:
    input:
        rep_asv_fasta = f"{DADA_DENOISE_DIR}/ASV.fasta",
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
        viral_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.viral.ASV.names"
    output:
        nonhost_nonviral_ASVs = f"{NEG_ALIGN_FILT_DIR}/nonhost.nonviral.ASV.fasta"
    log:
        f"{LOG_DIR}/04_nonhost_nonviral.log"
    threads: 8
    # TODO: consider consolidating with unmapped_bwa_mem and make one step??
    shell:
        r"""
        set -euo pipefail

        seqtk subseq "{input.rep_asv_fasta}" "{input.host_unmapped_names}" \
         | seqtk subseq - "{input.viral_unmapped_names}" \
         > "{output.nonhost_nonviral_ASVs}"
        """


#-------------------------------------
# 05 Positive Bacterial Alignment
#------------------------------------- 
rule bacterial_alignment:
    input:
        nonhost_nonviral_ASVs = f"{NEG_ALIGN_FILT_DIR}/nonhost.nonviral.ASV.fasta"
    output:
        bacterial_alignment = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.sam",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        bacterial_ASVs = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta"

    log:
        f"{LOG_DIR}/05_bact_map.log"
    threads: 8
    # Reverted filtering to use samtools view with positive selection for mapped reads (reverse of unmmaped fileter)
    #   was previously filtered via awk statement (this is more consistent)
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} "{BACT16S_REF}" "{input.nonhost_nonviral_ASVs}" > "{output.bacterial_alignment}" 2> "{log}"
        samtools view -F 4 "{output.bacterial_alignment}" | cut -f1 | sort -u > "{output.bacterial_names}"
        seqtk subseq "{input.nonhost_nonviral_ASVs}" "{output.bacterial_names}" > "{output.bacterial_ASVs}" 2> "{log}"
        """

#------------------------------------------------------------
# 06.1 IDTax (DECIPHER) Taxa Classification (NONFUNCTIONING)
#------------------------------------------------------------
rule IDTax_classification:
    input:
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        tax_ref = ID_TAX_REFERENCE
    output:
        taxonomy_table = f"{ID_TAX_DIR}/IDTax_Taxonomy.tsv",
    params:
        seq_table_dir = DADA_DENOISE_DIR,
        idtax_class_dir = ID_TAX_DIR,
    threads: 8
    log:
        f"{LOG_DIR}/06.1_IDTax.log"
    script:
        #would need to adjust script to deal with ealier Dada filter
        f"{SCRIPTS}/IdTaxaClassification.R"
        
#------------------------------
# 06.2 Mothur Taxa Classifcation 
#------------------------------
rule mothur_classify:
    input:
        bacterial_ASVs = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.fasta",
        mothur_ref = MOTHUR_REFERENCE,
        mothur_tax_file = MOTHUR_TAX_FILE
    output:
        taxfile = f"{MOTHUR_TAX_DIR}/bacterial.ASV.ncbi20.wang.taxonomy",
        tax_summary = f"{MOTHUR_TAX_DIR}/bacterial.ASV.ncbi20.wang.tax.summary",
        tax_summary_copy = f"{TAX_OUTPUT_DIR}/bacterial.ASV.ncbi20.wang.tax.summary"
    log:
        # TODO: determine why mothur still outputs additional .logfile
        #       piping properly into desired log but also creates its own in pwd
        #       currently just removing these after mothur run
        f"{LOG_DIR}/06.2_mothur_class.log"
    threads: 8
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
        
        # ---- Copy Taxonomy Summary to Output ---
        cp {output.tax_summary} {TAX_OUTPUT_DIR}/
        
        # ---- Remove Redundant Mothur Logs ----
        rm *.logfile
        """

#-------------------------------
# 07 Phyloseq Taxonomy Analysis
#-------------------------------
rule phyloseq_analysis:
    input:
        #seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        norm_seq_table = f"{NORM_COUNT_DIR}/NormSeqTable.tsv",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
        taxfile = f"{MOTHUR_TAX_DIR}/bacterial.ASV.ncbi20.wang.taxonomy",
    output:
        abunXtype_plot = f"{TAX_OUTPUT_DIR}/abunXtype.png",
        abunXsample_plot = f"{TAX_OUTPUT_DIR}/abunXsample.png",
        abunXtypeXsample_plot = f"{TAX_OUTPUT_DIR}/abunXtypeXsample.png"
    log:
        f"{LOG_DIR}/07_phyloseq.log"
    script:
        f"{SCRIPTS}/PhyloseqTaxAnalysis.R"


#--------------------------------
# 08 Manual ASV Taxonomy Counts
#--------------------------------
rule asv_tax_counts:
    input:
        taxfile = f"{MOTHUR_TAX_DIR}/bacterial.ASV.ncbi20.wang.taxonomy"
    output:
        tax_count = f"{TAX_OUTPUT_DIR}/bacterial.ASV.{{level}}.count"
    params:
        min_count = lambda wc: TAX_MIN[wc.level],
        cutoff = lambda wc: TAX_FIELDS[wc.level]
    log:
        f"{LOG_DIR}/08_manual_mothur_count/08_manual_mothur_count.{{level}}.log"
    shell:
        r"""
        set -euo pipefail
        exec 2> "{log}"

        # taxonomy -> (strip confidences) -> lineage col -> split by ';' -> count filtered taxa
        sed -E 's/\([^()]*\)//g' "{input.taxfile}" \
        | awk -F"\t" '{{print $2}}' \
        | cut -d ';' -f {params.cutoff} \
        | sort \
        | uniq -c \
        | awk -v m="{params.min_count}" '$1>m{{print $1, $2}}' \
        | sed -E 's/^ +//' \
        | sort -k2 \
        > "{output.tax_count}"
        """

        
# TODO: get to work with post-ASV collapse read counts
#       probably will have to do with R script (subsetting sequence table)
#-----------------------------
# -01 Read Count Calculations
#-----------------------------
rule read_counts:
    input:
        norm_raw_r1_fqs = expand(f"{NORM_RAW_DIR}/{{s}}_R1_001.fastq", s=SAMPLES),
        selected_r1_fqs = expand(f"{UMI_SELECT_DIR}/Selected.{{s}}.UMI_R1.fastq", s=SAMPLES),
        deduped_fqs = expand(f"{UMI_DEDUP_DIR}/Deduped.{{s}}.fastq", s=SAMPLES),

        dada_read_counts = f"{TRACK_DIR}/dada_read_counts.tsv",
        seq_table = f"{DADA_DENOISE_DIR}/SeqTable.tsv",
        
        host_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.host.ASV.names",
        viral_unmapped_names = f"{NEG_ALIGNMENT_DIR}/unmapped.viral.ASV.names",
        bacterial_names = f"{POS_ALIGNMENT_DIR}/bacterial.ASV.names",
    output:
        combined_read_counts = f'{TRACK_DIR}/combined_read_counts.tsv'
    params:
        samples = SAMPLES,
        normalized = NORM_RAW_DIR,
        selected = UMI_SELECT_DIR,
        deduped = UMI_DEDUP_DIR
    log:
        f"{LOG_DIR}/-01_read_count.log"
    script:
        "scripts/ReadCountCompilation.R"




#----------------------------
# -00 FASTQ to FASTA (Defunct)
#----------------------------   
# TODO: check if mothur really requires FASTA over FASTQ
rule fastq_to_fasta:
    input:
        fastq = f"{POS_ALIGNMENT_DIR}/bacterial.{{s}}.fastq"
    output:
        fasta = f"{POS_ALIGNMENT_DIR}/bacterial.{{s}}.fasta"
    log:
        f"{LOG_DIR}/-02_fq2fa/-02_fq2fa.{{s}}.log"
    # TODO: address seqkit dependency (most likely using singularity)
    shell:
        r"""
        set -euo pipefail
        seqkit fq2fa "{input.fastq}" > "{output.fasta}" 2> "{log}"
        """


#----------------------------------------
# -00 Trimming with Trimmomatic (Defunct)
#----------------------------------------   
rule trim_selected_r1:
    input:
        selected_fastq = f"{UMI_SELECT_DIR}/SelectedReads.{{s}}.R1.fastq",
        init = f"{OUT}/.init.done"
    output:
        trimmed_fastq = f"{TRIMMED_DIR}/trimmed.{{s}}.R1.fastq"
    params:
        trim_args = TRIM_ARGS_STR
    log:
        f"{LOG_DIR}/02_trim/02_trim.{{s}}.log"
    threads: 4
    # TODO: figure out trimmomatic dependency
    shell:
        r"""
        set -euo pipefail

        trimmomatic_jar="$TRIMMOJAR"
        java -jar "$trimmomatic_jar" SE -threads {threads} -phred33 "{input.selected_fastq}" "{output.trimmed_fastq}" "{params.trim_args}"
        """
