import os
import sys

import pyfaidx
from os.path import join


def main():
    fa_msa = sys.argv[1]
    ref = sys.argv[2]
    outdir = sys.argv[3]

    os.makedirs(outdir, exist_ok=True)

    if os.stat(fa_msa).st_size == 0:
        return

    fh_msa = pyfaidx.Fasta(fa_msa)

    fa_records = {rec.name: str(rec) for rec in fh_msa}

    ref_name = None
    for name in fa_records:
        if ref in name:
            ref_name = name

    rec_ref = fa_records.get(ref_name)
    if not rec_ref:
        return

    for name, rec in fa_records.items():
        if name == ref_name:
            continue

        with open(join(outdir, f"{ref_name}+{name}.fa"), "w") as output_handle:
            for rname, rseq in [(ref_name, rec_ref), (name, rec)]:
                print(f">{rname}", file=output_handle)
                print(rseq, file=output_handle)


if __name__ == "__main__":
    raise SystemExit(main())
