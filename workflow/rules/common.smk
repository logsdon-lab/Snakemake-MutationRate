import re
import json
from collections import defaultdict
from os.path import join, basename, splitext, dirname


RGX_CHROM = re.compile(config["regex_ref_chrom"])


def get_reference_beds():
    references = {}
    for ref in config["reference"]:
        cfg_ref = {}
        cfg_ref["path"] = ref["path"]
        if isinstance(ref["bed"], str):
            beds = [ref["bed"]]
        elif isinstance(ref["bed"], list):
            beds = ref["bed"]
        else:
            raise ValueError("No correct bed type.")

        bed_map = {}
        for bed in beds:
            name, _ = splitext(basename(bed))
            bed_map[name] = bed
        cfg_ref["bed"] = bed_map

        references[ref["name"]] = cfg_ref
    return references


def get_chrom(name: str) -> str | None:
    mtch = re.search(RGX_CHROM, name)

    if mtch:
        return mtch.group(1)
    else:
        return None


def get_window_bed_record(window: str) -> str:
    # Need to perform right split in case multiple coordinates.
    # chr1:1000-100000:1-10000
    chrom, coords = window.rsplit(":", 1)
    st, end = coords.split("-")
    return f"{chrom}\\t{st}\\t{end}\\n"
