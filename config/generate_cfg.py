import os
import sys
import yaml
import copy
import glob


def main():
    with open(
        os.path.join(os.path.dirname(__file__), "config_template.yaml"), "rt"
    ) as fh:
        cfg = yaml.safe_load(fh)

    sample = cfg["samples"].pop(0)

    cfg_samples = []
    for fa_cen in glob.glob(
        "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/data/cens/*.fasta"
    ):
        cfg_sample = copy.deepcopy(sample)
        # {sm}.cen.fasta
        sm, _, _ = os.path.basename(fa_cen).split(".")
        cfg_sample["name"] = sm
        cfg_sample["path"] = fa_cen
        cfg_samples.append(cfg_sample)

    cfg_samples.extend(cfg["samples"])
    cfg["samples"] = cfg_samples

    yaml.safe_dump(cfg, sys.stdout)


if __name__ == "__main__":
    raise SystemExit(main())
