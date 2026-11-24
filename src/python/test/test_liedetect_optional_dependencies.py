"""Check that LieDetect imports only the dependencies required for each step."""

import subprocess
import sys
import textwrap


def test_optional_dependencies_load_only_when_needed():
    script = textwrap.dedent(
        """
        import importlib.abc
        import sys

        blocked = {"autograd", "pymanopt", "sklearn"}

        class BlockOptionalDependencies(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname.split(".")[0] in blocked:
                    raise ModuleNotFoundError(fullname, name=fullname)

        sys.meta_path.insert(0, BlockOptionalDependencies())

        from gudhi.datasets.linear_orbits import sample_from_lie_group_rep
        from gudhi.liedetect import OrbitFitter, Orthonormalize

        points = sample_from_lie_group_rep("torus", ((1, 2),), 24)
        fitter = OrbitFitter(Orthonormalize(points).fit_transform())
        try:
            fitter.lie_pca(5, 1)
        except ImportError as exc:
            assert "requires scikit-learn" in str(exc)
        else:
            raise AssertionError("Lie-PCA should require scikit-learn")

        blocked.remove("sklearn")
        fitter.lie_pca(5, 1)
        for missing in ("autograd", "pymanopt"):
            blocked.clear()
            blocked.add(missing)
            try:
                fitter.closest_algebra("torus", group_dim=1, frequency_max=2)
            except ImportError as exc:
                assert "requires Autograd and Pymanopt" in str(exc)
                assert exc.__cause__.name == missing
            else:
                raise AssertionError(f"Algebra fitting should require {missing}")
        """
    )
    result = subprocess.run(
        [sys.executable, "-B", "-c", script], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
