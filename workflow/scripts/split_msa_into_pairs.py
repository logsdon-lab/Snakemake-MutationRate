import os
import sys

from os.path import join
from typing import Generator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    fa_msa = sys.argv[1]
    ref = sys.argv[2]
    outdir = sys.argv[3]

    os.makedirs(outdir, exist_ok=True)
    fh_msa: Generator[SeqRecord, None, None] = SeqIO.parse(fa_msa, "fasta")
    fa_records = {rec.id: rec for rec in fh_msa}

    ref_name = None
    for name in fa_records:
        if ref in name:
            ref_name = name

    rec_ref = fa_records.get(ref_name)
    if not rec_ref:
        raise ValueError(f"Missing reference fasta from {fa_msa} records. {fa_records.keys()}")

    for name, rec in fa_records.items():
        if name == ref_name:
            continue
        
        with open(join(outdir, f"{ref_name}+{name}.fa"), "w") as output_handle:
            SeqIO.write(
                [rec_ref, rec],
                output_handle,
                "fasta"
            )


if __name__ == "__main__":
    raise SystemExit(main())
