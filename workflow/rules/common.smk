import re
import json
from collections import defaultdict
from os.path import join, basename, splitext, dirname


RGX_CHROM = re.compile(r"(?:chr|hsa)([0-9XY]{1,2})")


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
