"""Execute the documented workflow and check directed reconstruction coverage."""

from pathlib import Path
import re
import textwrap

import numpy as np
import pytest

from gudhi.liedetect.orbit_fitter import print_hausdorff_distance


def test_user_manual_examples():
    manual = Path(__file__).resolve().parents[1] / "doc" / "liedetect_user.rst"
    examples = manual.read_text(encoding="utf-8").split("Generate synthetic data", 1)[1]
    blocks = re.findall(r"^\.\. code-block:: python\n\n((?:    [^\n]*\n|\n)+)", examples, re.MULTILINE)
    assert len(blocks) == 2
    scope = {}
    exec(compile(textwrap.dedent(blocks[0]), str(manual), "exec"), scope)
    assert scope["pts"].shape == (100, 4)
    exec(compile(textwrap.dedent(blocks[1]), str(manual), "exec"), scope)
    assert scope["rep_type"] == (3,)
    assert scope["orbit"].shape == (1000, 3)
    assert scope["orbit_fitter"].hausdorff_distances_[0] < 0.2


def test_directed_reconstruction_coverage_detects_a_missing_point():
    measured = np.array([[1.0, 0.0], [-1.0, 0.0]])
    reconstructed = measured[:1]
    assert print_hausdorff_distance(reconstructed, measured, verbose=False) == 0
    assert print_hausdorff_distance(measured, reconstructed, verbose=False) == pytest.approx(2)
