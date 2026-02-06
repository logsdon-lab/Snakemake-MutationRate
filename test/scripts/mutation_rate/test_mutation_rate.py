import os
import pytest

from ..common.integration import run_integration_test

CMD = os.path.basename(os.path.dirname(__file__))

SRC_MUT_RATE = "workflow/scripts/calculate_mutation_rate.py"
INPUT_DIR = f"test/scripts/{CMD}/input"
EXPECTED_DIR = f"test/scripts/{CMD}/expected"


@pytest.mark.parametrize(
    ["divergence", "divergence_times", "expected"],
    [
        (
            os.path.join(INPUT_DIR, f"{clade}.tsv.gz"),
            os.path.join(INPUT_DIR, f"{clade}_divergence.json"),
            os.path.join(EXPECTED_DIR, f"{clade}.bedpe.gz"),
        )
        for clade in ("chr8@9", "chrX@8", "chrY@2", "chrY@2_change_chr")
    ]
)
def test_mutation_rate(divergence, divergence_times, expected):
    cmd = ("python", SRC_MUT_RATE, "-i", divergence, "-d", divergence_times)
    run_integration_test(
        *cmd,
        expected_output=expected,
        overwrite_output=False
    )
