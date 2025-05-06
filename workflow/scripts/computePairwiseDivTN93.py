import os
import argparse
from Bio import AlignIO
from tn93 import tn93


def calculate_gc_prop(aln: str) -> int:
    return round((aln.count("G") + aln.count("C")) / (len(aln) - aln.count("-")), 3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fa_msa", nargs="*", type=str)
    parser.add_argument("--mode", nargs="?", const=1, default=2)
    args = parser.parse_args()

    files = args.fa_msa if args.fa_msa else []
    for file in files:
        pair_1, _, pair_2 = os.path.splitext(os.path.basename(file))[0].partition("+")
        alignment = AlignIO.read(file, "fasta")
        alignment1 = str(alignment[0].seq.upper())
        alignment2 = str(alignment[1].seq.upper())

        propGC_seq1 = calculate_gc_prop(alignment1)
        propGC_seq2 = calculate_gc_prop(alignment2)
        propGC = round((propGC_seq1 + propGC_seq2) / 2, 3)
        # https://en.wikipedia.org/wiki/Genetic_distance
        # https://github.com/sdwfrost/libtn93?tab=readme-ov-file#api
        distance_tn93 = tn93(alignment1, alignment2, len(alignment1), args.mode, 100)
        l_out = [pair_1, pair_2, propGC_seq1, propGC_seq2, propGC, distance_tn93]
        print("\t".join([str(x) for x in l_out]))
