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

def autheticate_and_extract_umi(seq2: str, qual2: str, motif: str, umi_len: int, max_offset: int):
    motif_len = len(motif)

    for offset in range(max_offset+1):
        end_umi = offset + umi_len
        end_motif = end_umi + motif_len
        if end_motif > len(seq2):
            break
        
        umi = seq2[offset:end_umi]
        motif_candidate = seq2[end_umi:end_motif]

        if motif_candidate == motif:
            umi_qual = qual2[offset:end_umi]
            return offset, umi, umi_qual
    
    return None

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--r2-primer-motif", required=True, help="Primer/adapter motif immediately following UMI in R2 used to authenticate R1-R2 read pairings")
    ap.add_argument("--umi-len", type=int ,required=True)
    ap.add_argument("--max-offset", type=int, default=0, help="The maximum number of leading junk bases allowed in R2 when searching for [UMI][R2-primer-motif]")
    ap.add_argument("--out-umi-tsv", required=True)
    ap.add_argument("--out-sel-names", required=True)
    ap.add_argument("--out-umi-r1", required=True)
    ap.add_argument("--out-sel-r1", default=None, help="Optional sanity-check fastq output of selected authentic R1 reads without concatenated UMIs" )
    args = ap.parse_args()

    r2_primer_motif = args.r2_primer_motif
    motif_len = len(r2_primer_motif)
    umi_len = args.umi_len
    max_offset = args.max_offset

    total = authentic = motif_fail = short_fail = desync_fail = 0

    with \
        open(args.r1, "r") as r1_fq, open(args.r2, "r") as r2_fq, \
        open(args.out_umi_tsv, "w") as out_umi_tsv, \
        open(args.out_sel_names, "w") as out_sel_names, \
        open(args.out_umi_r1, "w") as out_umi_r1, \
        (open(args.out_sel_r1, "w") if args.out_sel_r1 else open ("/dev/null", "w")) as out_sel_r1:

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


            extracted = autheticate_and_extract_umi(seq2, qual2, r2_primer_motif, umi_len, max_offset)
            if not extracted:
                motif_fail += 1
                continue

            offset, umi, umi_qual = extracted

            out_umi_tsv.write(f"{name1}\t{umi}\t{umi_qual}\n")
            out_sel_names.write(f"{name1}\n")

            if args.out_sel_r1:
                out_sel_r1.write(f"{head1}\n{seq1}\n{plus1}\n{qual1}\n")

            out_umi_r1.write(f"{head1} UMI:{umi}\n{umi}{seq1}\n{plus1}\n{umi_qual}{qual1}\n")
            authentic += 1
    

    print(f"Total read pairs:           {total}", file=sys.stderr)
    print(f"Authentic read pairs:       {authentic}", file=sys.stderr)
    print(f"Motif mismatches:           {motif_fail}", file=sys.stderr)
    print(f"Too shorts:                 {short_fail}", file=sys.stderr)
    print(f"Desync failures:            {desync_fail}", file=sys.stderr)

if __name__ == "__main__":
    main()

