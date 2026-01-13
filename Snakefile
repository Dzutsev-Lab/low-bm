# TODO: shift global variable declaration to configile
import os, re, glob
from pathlib import Path
# Input and Output Directories
RAW="Exp_Data/Single_Sample"
OUT="Exp_Output/Single_Sample"

ORIG_DIR  = f"{OUT}/original.fastq"
CLEAN_DIR = f"{OUT}/clean.fastq"
TAX_DIR   = f"{OUT}/taxonomy"
LOG_DIR   = f"{OUT}/logs"

# Sequencing metadata
PRIMER_MOTIF="GGACTAC"
UMI_LEN="16"

# Reference Data
HOST_REF="Ref_Data/Mus_musculus.GRCm38.cdna.all.fa" 
VIRAL_REF="Ref_Data/all.viral.fna"
BACT16S_REF="Ref_Data/all.rrna.bacteria"
# mothur references
MOTHUR_TEMPLATE="Ref_Data/ncbi20.fasta"
MOTHUR_TAX="Ref_Data/ncbi20.tax"

# Trimming Parameters
TRIMMOMATIC_JAR="$TRIMMOJAR"
TRIM_ARGS="AVGQUAL:30" 

# Negative Control Filtering Parameters
# TODO: figure out how to format for snakemake without config file
NEG_DB_TOP_K="10000"
NEG_AS_KEEP_LT="160"   

# Taxonomy Cutoffs
# TODO: figure out how to format for snakemake without config file
TAX_MIN = ()

THREADS = 2


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
        taxfile = expand(f"{TAX_DIR}/bacterial.{{s}}.ncbi20.wang.taxonomy", s=SAMPLES)

#----------------------------
# Construct Directories
#----------------------------
rule init_dirs:
    output:
        touch(f"{OUT}/.init.done")
    shell:
        r"""
        mkdir -p "{OUT}" "{ORIG_DIR}" "{CLEAN_DIR}" "{TAX_DIR}" "{LOG_DIR}"
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

#---------------------------------------------
# UMI extraction from R2 and Select R1 Reads
#---------------------------------------------
rule umi_extract_select_r1:
    input:
        r1 = f'{RAW}/{{s}}_R1_001.fastq',
        r2 = f'{RAW}/{{s}}_R2_001.fastq',
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
    log:
        f"{LOG_DIR}/02_trim/02_trim.{{s}}.log"
    threads: THREADS
    # TODO: figure out trimmomatic dependency
    shell:
        r"""
        set -euo pipefail

        trimmomatic_jar="$TRIMMOJAR"
        java -jar "$trimmomatic_jar" SE -threads {threads} -phred33 "{input.selected_fastq}" "{output.trimmed_fastq}" "{TRIM_ARGS}"
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
    # TODO: work to understand why we are extracting read names with awk here and did it with samtools with viral and host alignments
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} "{BACT16S_REF}" "{input.nonhost_nonviral_reads}" > "{output.bacterial_alignment}" 2> "{log}"
        awk '($1 !~ /^@/) && ($3 != "*") {{print $1}}' "{output.bacterial_alignment}" | sort -u > "{output.bacterial_names}"
        seqtk subseq "{input.nonhost_nonviral_reads}" "{output.bacterial_names}" > "{output.bacterial_reads}" 2> "{log}"
        """

#-------------------------
# Read Count Calculations
#-------------------------
rule read_counts:
    input:
        # need all samples processed to execute rule once, 
        #   consolidating data from all samples into one stats file
        raw_r1s = expand(f'{RAW}/{{s}}_R1_001.fastq', s=SAMPLES),
        host_alignment = expand(f"{CLEAN_DIR}/Alignment/host.{{s}}.sam", s=SAMPLES),
        bacterial_reads = expand(f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.fasta", s=SAMPLES)
    output:
        # TODO: consider numbering out the output subdirectories to keep track of order
        count_stats = f'{OUT}/Counts/all_sample.count.stats.tsv'
    params:
        samples=SAMPLES_STR,
        raw=RAW,
        clean=CLEAN_DIR
    shell:
        # TODO: consider shifting this to a series of python functions for easier function
        r"""
        {{
            set -euo pipefail
            echo -e "ID\tTotal_reads\tHuman_mapped_reads\tClean_reads"
            for s in {params.samples}; do
                r1="{params.raw}/${{s}}_R1_001.fastq"
                sam="{params.clean}/Alignment/host.${{s}}.sam"
                fa="{params.clean}/Bacterial/bacterial.${{s}}.fasta"

                total=$(($(wc -l < "$r1")/4))
                hum_mapped=$(samtools view -F 4 "$sam" | wc -l)
                clean=$(grep -c '^>' "$fa" 2>/dev/null || true)

                echo -e "${{s}}\t${{total}}\t${{hum_mapped}}\t${{clean}}"
            done
        }} > "{output.count_stats}"
        """

#-------------------------
# Mothur CLassifcation
#-------------------------
rule mothur_classify:
    input:
        bacterial_reads = f"{CLEAN_DIR}/Bacterial/bacterial.{{s}}.fasta",
        mothur_template = MOTHUR_TEMPLATE,
        mothur_tax = MOTHUR_TAX
    output:
        taxfile = f"{TAX_DIR}/bacterial.{{s}}.ncbi20.wang.taxonomy"
        tax_summary = f"{TAX_DIR}/bacterial.{{s}}.ncbi20.wang.tax.summary""
    log:
        f"{LOG_DIR}/07_mothur_class/07_mothur_class.{{s}}.log"
    threads: THREADS
    shell:
        r"""
        set -euo pipefail
        mothur "#set.dir(output={TAX_DIR}); classify.seqs(fasta={input.bacterial_reads}, template={input.mothur_template}, taxonomy={input.mothur_tax}, mothod=wang, cutoff=0, processors={threads})" \
            2> "{log}"
        """
