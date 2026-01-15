# TODO: shift global variable declaration to configile
import os, re, glob
from pathlib import Path

configfile: "config.yaml"

THREADS = config["threads"]

# Input and Output Directories
RAW = config["raw_dir"]
OUT = config["out_root"]

ORIG_DIR = f"{OUT}/original.fastq"
NORM_RAW_DIR = f"{OUT}/raw.norm"
CLEAN_DIR = f"{OUT}/clean.fastq"
TAX_DIR = f"{OUT}/taxonomy"
LOG_DIR = f"{OUT}/logs"

# Sequencing metadata
PRIMER_MOTIF = config["primer_motif"]
UMI_LEN = config["umi_len"]

# Reference Data
HOST_REF = config["host_ref"]
VIRAL_REF = config["viral_ref"]
BACT16S_REF = config["bact16s_ref"]
# mothur references
MOTHUR_REFERENCE = config["mothur_reference"]
MOTHUR_TAX = config["mothur_tax"]

# Trimming Parameters
TRIMMOMATIC_JAR = config["trimmomatic_jar"]
TRIM_ARGS = config["trim_args"]
TRIM_ARGS_STR =  " ".join(TRIM_ARGS)

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
        tax_count = expand(f"{TAX_DIR}/bacterial.{{s}}.counts/{{s}}.{{level}}.count", s=SAMPLES, level=TAX_LEVELS)

#----------------------------
# Construct Directories
#----------------------------
rule init_dirs:
    output:
        touch(f"{OUT}/.init.done")
    shell:
        r"""
        mkdir -p "{OUT}" "{ORIG_DIR}" "{NORM_RAW_DIR}" "{CLEAN_DIR}" "{TAX_DIR}" "{LOG_DIR}"
        touch "{output}"
        """


#----------------------------
# Record Sample Names
#----------------------------
rule sample_names:
    input:
        f"{OUT}/.init.done"
    output:
        f"{OUT}/sample.names"
    run:
        Path(output[0]).write_text("\n".join(SAMPLES) + "\n")

#----------------------------
# Link Original Files
#----------------------------
# TODO: implement


#-----------------------------------
# Normalizing and Unzipping FASTQs
#-----------------------------------
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
        init = f"{OUT}/.init.done"
    output:
        r1_norm = f"{NORM_RAW_DIR}/{{s}}_R1_001.fastq",
        r2_norm = f"{NORM_RAW_DIR}/{{s}}_R2_001.fastq"
    log:
        f"{LOG_DIR}/00_norm/00_norm.{{s}}.log"
    threads: THREADS
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

#---------------------------------------------
# UMI extraction from R2 and Select R1 Reads
#---------------------------------------------
rule umi_extract_select_r1:
    input:
        r1 = f'{NORM_RAW_DIR}/{{s}}_R1_001.fastq',
        r2 = f'{NORM_RAW_DIR}/{{s}}_R2_001.fastq',
        sample_names = f"{OUT}/sample.names",
        init = f"{OUT}/.init.done"
    output:
        umi = f"{CLEAN_DIR}/UMISelection/UMIMap.{{s}}.tsv",
        selected_names = f"{CLEAN_DIR}/UMISelection/SelectedNames.{{s}}.names",
        selected_r1 = f"{CLEAN_DIR}/UMISelection/SelectedReads.{{s}}.R1.fastq"
    log:
        f"{LOG_DIR}/01_umi_select/01_umi_select.{{s}}.log"
    threads: THREADS
    # TODO: address seqtk dependency (most likely using singularity)
    shell:
        r"""
        set -euo pipefail
        exec 2> "{log}"

        R2_IN="{input.r2}"
        R1_IN="{input.r1}"

        # extract NAME<TAB>UMI from R2; also write selected names
        awk -v motif="{PRIMER_MOTIF}" -v L={UMI_LEN} '
            NR%4==1 {{ 
                hdr=$0; 
                name=substr($0,2)
                sub(/ .*/, "", name)
                next 
            }}
            NR%4==2 {{ 
                seq=$0
                if (match(seq, "([ACGTN]{{"L"}})" motif, m)) {{
                            print name "\t" m[1]
                }}
            }}
            ' "{input.r2}" | tee "{output.umi}" | cut -f1 > "{output.selected_names}"

        # subset R1 by selected names
        seqtk subseq "{input.r1}" "{output.selected_names}" > "{output.selected_r1}"
        """


#----------------------------
# Trimming with Trimmomatic
#----------------------------   
rule trim_selected_r1:
    input:
        selected_fastq = f"{CLEAN_DIR}/UMISelection/SelectedReads.{{s}}.R1.fastq",
        init = f"{OUT}/.init.done"
    output:
        trimmed_fastq = f"{CLEAN_DIR}/Trimming/trimmed.{{s}}.R1.fastq"
    params:
        trim_args = TRIM_ARGS_STR
    log:
        f"{LOG_DIR}/02_trim/02_trim.{{s}}.log"
    threads: THREADS
    # TODO: figure out trimmomatic dependency
    shell:
        r"""
        set -euo pipefail

        trimmomatic_jar="$TRIMMOJAR"
        java -jar "$trimmomatic_jar" SE -threads {threads} -phred33 "{input.selected_fastq}" "{output.trimmed_fastq}" "{params.trim_args}"
        """


#----------------------------
# FASTQ to FASTA
#----------------------------   
rule fastq_to_fasta:
    input:
        fastq = f"{CLEAN_DIR}/Trimming/trimmed.{{s}}.R1.fastq"
    output:
        fasta = f"{CLEAN_DIR}/Trimming/trimmed.{{s}}.R1.fasta"
    log:
        f"{LOG_DIR}/03_fq2fa/03_fq2fa.{{s}}.log"
    # TODO: address seqkit dependency (most likely using singularity)
    shell:
        r"""
        set -euo pipefail
        seqkit fq2fa "{input.fastq}" > "{output.fasta}" 2> "{log}"
        """


#----------------------------
# Reference Indexing
#----------------------------  
rule index_ref_bwa:
    input:
        reference = "{ref}"
    output:
        # indexing produces many files but .bwt is the key one (will use as sentinel)
        bwt = "{ref}.bwt"
    threads: THREADS
    shell:
        r"""
        set -euo pipefail
        bwa index "{input.ref}"
        """

#--------------------------------------------------------
# Contaminant Alignment and Unmapped Selection
#--------------------------------------------------------
rule unmapped_bwa_mem:
    input:
        trimmed_fasta = f"{CLEAN_DIR}/Trimming/trimmed.{{s}}.R1.fasta",
        # What is the lambda doing???
        reference_fasta = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag],
        reference_bwt = lambda wc: {"host": HOST_REF, "viral": VIRAL_REF}[wc.tag] + ".bwt"
    output:
        sam = f"{CLEAN_DIR}/Alignment/{{tag}}.{{s}}.sam",
        unmapped_names = f"{CLEAN_DIR}/Unmapped/{{tag}}.unmapped.{{s}}.names"
    log:
        f"{LOG_DIR}/04_map/04_map_{{tag}}.{{s}}.log"
    threads: THREADS
    # TODO: address BWA dependency (most likely using singularity)
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} "{input.reference_fasta}" "{input.trimmed_fasta}" > "{output.sam}" 2> "{log}"
        samtools view -f 4 "{output.sam}" | cut -f1 | sort -u > "{output.unmapped_names}"
        """

#-------------------------------------
# Unmapped Host and Viral Read Filter
#------------------------------------- 
rule nonhost_nonviral_reads:
    input:
        trimmed_fasta = f"{CLEAN_DIR}/Trimming/trimmed.{{s}}.R1.fasta",
        host_unmapped_names = f"{CLEAN_DIR}/Unmapped/host.unmapped.{{s}}.names",
        viral_unmapped_names = f"{CLEAN_DIR}/Unmapped/viral.unmapped.{{s}}.names"
    output:
        nonhost_nonviral_reads = f"{CLEAN_DIR}/Unmapped/nonhost.nonviral.{{s}}.fasta"
    log:
        f"{LOG_DIR}/05_nonhost_nonviral/05_nonhost_nonviral{{s}}.log"
    threads: THREADS
    # TODO: consider consolidating with unmapped_bwa_mem and make one step??
    shell:
        r"""
        set -euo pipefail

        seqtk subseq "{input.trimmed_fasta}" "{input.host_unmapped_names}" \
         | seqtk subseq - "{input.viral_unmapped_names}" \
         > "{output.nonhost_nonviral_reads}"
        """


#-------------------------------------
# Positive Bacterial Alignment
#------------------------------------- 
rule bacterial_alignment:
    input:
        nonhost_nonviral_reads = f"{CLEAN_DIR}/Unmapped/nonhost.nonviral.{{s}}.fasta"
    output:
        bacterial_alignment = f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.sam",
        bacterial_names = f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.names",
        bacterial_reads = f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.fasta"

    log:
        f"{LOG_DIR}/06_bact_map/06_bact_map.{{s}}.log"
    threads: THREADS
    # Reverted filtering to use samtools view with positive selection for mapped reads (reverse of unmmaped fileter)
    #   was previously filtered via awk statement (this is more consistent)
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} "{BACT16S_REF}" "{input.nonhost_nonviral_reads}" > "{output.bacterial_alignment}" 2> "{log}"
        samtools view -F 4 "{output.bacterial_alignment}" | cut -f1 | sort -u > "{output.bacterial_names}"
        seqtk subseq "{input.nonhost_nonviral_reads}" "{output.bacterial_names}" > "{output.bacterial_reads}" 2> "{log}"
        """

#-------------------------
# Read Count Calculations
#-------------------------
rule read_counts:
    input:
        # need all samples processed to execute rule once, 
        #   consolidating data from all samples into one stats file
        raw_r1s = expand(f'{NORM_RAW_DIR}/{{s}}_R1_001.fastq', s=SAMPLES),
        host_alignment = expand(f"{CLEAN_DIR}/Alignment/host.{{s}}.sam", s=SAMPLES),
        bacterial_reads = expand(f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.fasta", s=SAMPLES)
    output:
        # TODO: consider numbering out the output subdirectories to keep track of order
        read_counts = f'{OUT}/counts/all_sample.count.stats.tsv'
    params:
        samples=SAMPLES_STR,
        norm_raw=NORM_RAW_DIR,
        selected=f'{CLEAN_DIR}/UMISelection',
        trimmed=f'{CLEAN_DIR}/Trimming',
        unmapping=f'{CLEAN_DIR}/Alignment',
        bacterial=f'{CLEAN_DIR}/Bacterial'
    log:
        f"{LOG_DIR}/-01_read_count/-01_read_count.log"
    shell:
        # TODO: take into account secondary and supplementary alignments whne counting records in sam's
        r"""
        set -euo pipefail
        exec 2> "{log}"

        {{
            set -euo pipefail
        
            echo -e "ID\tRaw_reads\tSelected_reads\tTrimmed_reads\tPos_Viral_reads\tPos_Host_reads\tPos_Bacterial_reads"
            for s in {params.samples}; do

                raw_r1_fq="{params.norm_raw}/${{s}}_R1_001.fastq"
                selected_fq="{params.selected}/SelectedReads.${{s}}.R1.fastq"
                trimmed_fq="{params.trimmed}/trimmed.${{s}}.R1.fastq"
                viral_sam="{params.unmapping}/viral.${{s}}.sam"
                host_sam="{params.unmapping}/host.${{s}}.sam"
                bacterial_fa="{params.bacterial}/bacterial.${{s}}.fasta"

                raw_r1=$(($(wc -l < "$raw_r1_fq")/4))
                selected=$(($(wc -l < "$selected_fq")/4))
                trimmed=$(($(wc -l < "$trimmed_fq")/4))
                viral=$(samtools view -F 4 "$viral_sam" | wc -l)
                host=$(samtools view -F 4 "$host_sam" | wc -l)
                bacterial=$(grep -c '^>' "$bacterial_fa" 2>/dev/null || true)

                
                echo -e "${{s}}\t${{raw_r1}}\t${{selected}}\t${{trimmed}}\t-${{viral}}\t-${{host}}\t${{bacterial}}"
            done
        }} > "{output.read_counts}"
        """

#-------------------------
# Mothur CLassifcation
#-------------------------
rule mothur_classify:
    input:
        bacterial_reads = f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.fasta",
        mothur_ref = MOTHUR_REFERENCE,
        mothur_tax = MOTHUR_TAX
    output:
        taxfile = f"{TAX_DIR}/bacterial.{{s}}.ncbi20.wang.taxonomy",
        tax_summary = f"{TAX_DIR}/bacterial.{{s}}.ncbi20.wang.tax.summary"
    log:
        # TODO: determine why mothur still outputs additional .logfile
        #       piping properly into desired log but also creats its own in pwd
        f"{LOG_DIR}/07_mothur_class/07_mothur_class.{{s}}.log"
    threads: THREADS
    shell:
        r"""
        set -euo pipefail

        # ---- Guard: empty FASTA ----
        if ! grep -q '^>' "{input.bacterial_reads}"; then
            echo "[INFO] {wildcards.s}: no bacterial reads found; skipping mothur classify.seqs" > {log}
            # creating empty but valid output to continue Snakemake
            : > "{output.taxfile}"
            : > "{output.tax_summary}"
            exit 0
        fi

        # ---- Run Mothur ----
        mothur "#set.dir(output={TAX_DIR}); 
                set.logfile(name={log});
                classify.seqs(
                    fasta={input.bacterial_reads}, 
                    reference={input.mothur_ref}, 
                    taxonomy={input.mothur_tax}, 
                    method=wang, 
                    cutoff=0, 
                    processors={threads}
                )" >> {log} 2>&1
        """


#----------------------------------
# Taxonomy Classification Counts
#----------------------------------
rule tax_counts:
    input:
        taxonomy = f"{TAX_DIR}/bacterial.{{s}}.ncbi20.wang.taxonomy"
    output:
        tax_count = f"{TAX_DIR}/bacterial.{{s}}.counts/{{s}}.{{level}}.count"
    params:
        min_count = lambda wc: TAX_MIN[wc.level],
        cutoff = lambda wc: TAX_FIELDS[wc.level]
    log:
        f"{LOG_DIR}/08_tax_count/08_tax_count.{{level}}.{{s}}.log"
    shell:
        r"""
        set -euo pipefail
        exec 2> "{log}"

        # taxonomy -> (strip confidences) -> lineage col -> split by ';' -> count filtered taxa
        sed -E 's/\([^()]*\)//g' "{input.taxonomy}" \
        | awk -F"\t" '{{print $2}}' \
        | cut -d ';' -f {params.cutoff} \
        | sort \
        | uniq -c \
        | awk -v m="{params.min_count}" '$1>m{{print $1, $2}}' \
        | sed -E 's/^ +//' \
        | sort -k2 \
        > "{output.tax_count}"
        """