import argparse
import subprocess
import re
from Bio import SeqIO
from collections import OrderedDict

def seqID2name_mapping(tax_path):
    mapping = OrderedDict()
    with open(tax_path) as tax_file:
        for line in tax_file:
            line = line.strip()
            if not line:
                continue
            #splits line into seqence ID and lineage (ensures only split between those two)
            seqID, lineage = line.split(maxsplit=1)

            # removes any empty levels
            # why used over linearge.split(';') alone
            levels = [l for l in lineage.split(';') if l]

            species = None
            genus = None

            for l in levels:
                if l.startswith('s__') and len(l) > 3:
                    species = l[3:]
                    break
                if l.startswith('g__') and len(l) > 3:
                    genus = l[3:]
            
            lowest_level = species if species else genus

            if lowest_level:
                lowest_level = lowest_level.replace("_", " ").strip()
                lowest_level = lowest_level.strip(";,")
            
            mapping[seqID] = lowest_level

        return mapping
    
def name2taxID_mapping(names, data_dir, faulty_out):
    cmd = ['taxonkit', 'name2taxid', '--data-dir', data_dir]

    process_o = subprocess.Popen(cmd,
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 text=True)
    
    stdin = "\n".join(names) + "\n"
    out, err = process_o.communicate(stdin)

    if process_o.returncode != 0:
        raise RuntimeError(f"taxonkit failure: {err}")

    mapping = {}
    bad_lines = []

    for line in out.splitlines():
        print(line)
        line = line.strip()

        if not line:
            continue
        
        line_elements = line.split("\t")

        if len(line_elements) != 2:
            bad_lines.append(line)
            continue

        name, taxID = line_elements
        taxID = taxID.strip()

        if taxID == "0" or taxID == "":
            mapping[name] = None
        else:
            mapping[name] = taxID
    
    if bad_lines:
        with open(faulty_out, "w") as faulty_file:
            faulty_file.write("\n".join(bad_lines) + "\n")
    
    return mapping

def rewrite_fasta(fasta_in, fasta_out, seqID2taxID):
    with open(fasta_out, 'w') as fout:
        for record in SeqIO.parse(fasta_in, "fasta"):
            seqID = record.id
            taxID = seqID2taxID.get(seqID)
            if taxID:
                record.id = f"{seqID}|kraken:taxid|{taxID}"
                record.description = ""
                SeqIO.write(record, fout, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tax', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--faulty-out', required=True)
    args = parser.parse_args()

    print("Parsing taxonomy file...")
    seqID2name = seqID2name_mapping(args.tax)

    unique_names = sorted({n for n in seqID2name.values() if n})
    print(f"{len(unique_names)} unique names to resolve")

    print("Resolving names with taxonkit...")
    name2taxID = name2taxID_mapping(unique_names, args.data_dir, args.faulty_out)

    seqID2taxID = {}
    failed = []

    for seqID, name in seqID2name.items():
        taxID = name2taxID.get(name)
        if taxID:
            seqID2taxID[seqID] = taxID
        else:
            failed.append((seqID, name))
    
    print(f"{len(failed)} sequences failed taxid resolution")
    print("Writing Kraken-formatted FASTA...")
    rewrite_fasta(args.fasta, args.out, seqID2taxID)

    