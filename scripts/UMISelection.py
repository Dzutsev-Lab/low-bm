#!/usr/bin/env python
import argparse
import sys

def fastq_iterator(handle):
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not quality:
            raise ValueError(f"Truncated FASTQ record: no qualtiy entry found for read {header}")
        # Yield helps you keep your place in an local iteration when calling a function
        # iterator placeholder when returning values
        yield header.rstrip("\n"), sequence.rstrip("\n"), plus.rstrip("\n"), quality.rstrip("\n")

def read_name(header_line: str) -> str:
    # header format: $NB####:... 1:N:0:...
    # header token to be header up to first whitespace
    if not header_line.startswith("@"):
        raise ValueError(f"Bad FASTQ header : {header_line}")
    return header_line[1:].split()[0] # taking the @ symbol out and taking first string from split

def autheticate_and_extract_umi(seq2: str, qual2: str, r2_primer_motif: str, r2_primer_skip: bool, poly_G_threshold: float, umi_len: int, max_offset: int):
    motif_len = len(r2_primer_motif)

    for offset in range(max_offset+1):
        end_umi = offset + umi_len
        end_motif = end_umi + motif_len
        if end_motif > len(seq2):
            break
        
        umi = seq2[offset:end_umi]
        motif_candidate = seq2[end_umi:end_motif]

        if motif_candidate == r2_primer_motif:
            umi_qual = qual2[offset:end_umi]
            return offset, umi, umi_qual
    
    # In case where skipping R2 primer motif check for authenticity
    #   Put after offset check to ensure we do catch all possible cases of junk beginning inserts
    #   If no authetic motif found, just return umi as first umi-length characters
    #   Offset set to -1 in this case to indicate no motif was found
    if r2_primer_skip:
        umi = seq2[:umi_len]
        umi_qual =qual2[:umi_len]
        # checking poly-G content of UMI
        if not (umi.count("G") / umi_len > poly_G_threshold):
            return -1, umi, umi_qual
    return None

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-name", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--r2-primer-motif", required=True, help="Primer/adapter motif immediately following UMI in R2 used to authenticate R1-R2 read pairings")
    ap.add_argument("--r2-primer-skip", action='store_true', help="Skipping mofit-based authentification of R1-R2 read pairs")
    ap.add_argument("--poly-G-threshold", type=float, default=1.0, help="The maximum fraction of G's in UMI before filteration, only applicable with --r2-primer-skip")
    ap.add_argument("--umi-len", type=int ,required=True)
    ap.add_argument("--max-offset", type=int, default=0, help="The maximum number of leading junk bases allowed in R2 when searching for [UMI][R2-primer-motif]")
    ap.add_argument("--out-count-summary", required=True, help="Output summary TSV file with read count columns Raw_reads (input R1), Selected_reads (ouput)")
    ap.add_argument("--out-umi-r1", required=True)
    args = ap.parse_args()

    r2_primer_motif = args.r2_primer_motif
    r2_primer_skip = args.r2_primer_skip
    poly_G_threshold = args.poly_G_threshold
    motif_len = len(r2_primer_motif)
    umi_len = args.umi_len
    max_offset = args.max_offset

    total = authentic = motif_fail = polyG_fail = short_fail = desync_fail = 0

    with \
        open(args.r1, "r") as r1_fq, open(args.r2, "r") as r2_fq, \
        open(args.out_umi_r1, "w") as out_umi_r1, \
        open(args.out_count_summary, "w") as count_summary:

        for (head1, seq1, plus1, qual1), (head2, seq2, plus2, qual2) in zip(fastq_iterator(r1_fq), fastq_iterator(r2_fq)):
            total += 1

            name1 = read_name(head1)
            name2 = read_name(head2)

            if name1 != name2:
                desync_fail += 1
                raise ValueError(f"R1/R2 Unmatched: {name1} != {name2}")
            
            if len(seq2) < umi_len + motif_len:
                short_fail += 1
                continue


            extracted = autheticate_and_extract_umi(seq2 = seq2, 
                                                    qual2 = qual2, 
                                                    r2_primer_motif = r2_primer_motif, 
                                                    r2_primer_skip = r2_primer_skip, 
                                                    poly_G_threshold = poly_G_threshold, 
                                                    umi_len = umi_len, 
                                                    max_offset = max_offset)
            if not extracted:
                if r2_primer_skip:
                    polyG_fail += 1
                else:
                    motif_fail += 1
                continue
            elif extracted[0] == -1:
                # motif search failed but we are skipping as mofit check as hard filter,
                # but must have passed poly-G content check to be here
                # will count as both a motif failure and an authentic read pair since we are keeping it despite mofit failure
                motif_fail += 1
                authentic += 1
            else:
                authentic += 1

            offset, umi, umi_qual = extracted

            out_umi_r1.write(f"{head1} UMI:{umi}\n{umi}{seq1}\n{plus1}\n{umi_qual}{qual1}\n")

        count_summary.write(f"SampleID\tRaw_reads\tSelected_reads\n")
        count_summary.write(f"{args.sample_name}\t{total}\t{authentic}\n")
            
    

    print(f"Total read pairs:           {total}", file=sys.stderr)
    print(f"Fully Authentic read pairs: {authentic}", file=sys.stderr)
    print(f"High G content UMI:         {polyG_fail}", file=sys.stderr)
    print(f"Motif mismatches:           {motif_fail}", file=sys.stderr)
    print(f"Too shorts:                 {short_fail}", file=sys.stderr)
    print(f"Desync failures:            {desync_fail}", file=sys.stderr)

if __name__ == "__main__":
    main()

