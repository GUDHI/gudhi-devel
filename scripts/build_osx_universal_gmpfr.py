import os
import subprocess
import glob
from pathlib import Path

# Oldest macOS version supported by brew bottles
OSX_VERSION = "ventura"

# Some color to distinguish outputs
GREEN = "\033[92m"
END = "\033[0m"

def command_output(cmd: str) -> str:
    print(f"{GREEN}{cmd}{END}")
    return subprocess.check_output(cmd, shell=True, text=True)

def fetch_and_extract(arch: str, pkg: str) -> None:
    """Fetch a brew bottle for given arch and package, then extract it.
    Assumes that the user has enough rights to run brew fetch.
    Parameters:
        arch (str): Architecture ('x86_64' or 'arm64').
        pkg (str): Brew package name.
    """
    if arch not in ["x86_64", "arm64"]:
        raise RuntimeError(f"Unknown architecture {arch}")
    _ = command_output(f"brew fetch --bottle-tag={arch}_{OSX_VERSION} {pkg}")
    # Get tar name from brew fetch message
    tarball = command_output(f"brew --cache --bottle-tag={arch}_{OSX_VERSION} {pkg}")
    # Here we should only find "bottle.tar.gz" in the command output
    assert "bottle.tar.gz" in tarball
    _ = command_output(f"tar xf {tarball}")


cwd = Path.cwd()

for arch in ["x86_64", "arm64"]:
    Path(f"deps-{arch}").mkdir(exist_ok=True)
    os.chdir(f"deps-{arch}")
    fetch_and_extract(arch, "gmp")
    fetch_and_extract(arch, "mpfr")
    os.chdir(cwd)

# Merge libraries in a universal directory: deps-uni
uni_lib = Path("deps-uni/lib")
uni_lib.mkdir(parents=True, exist_ok=True)

def merge_lib(pattern: str, pkg: str) -> Path:
    amd64_lib = glob.glob(f"deps-x86_64/{pkg}/*/lib/{pattern}")[0]
    arm64_lib = glob.glob(f"deps-arm64/{pkg}/*/lib/{pattern}")[0]
    out = uni_lib / os.path.basename(amd64_lib)
    _ = command_output(f"lipo -create {amd64_lib} {arm64_lib} -output {out}")
    _ = command_output(f"install_name_tool -id {cwd}/{out} {out}")
    return out


gmp = merge_lib("libgmp.*.dylib", "gmp")
gmpxx = merge_lib("libgmpxx.*.dylib", "gmp")
mpfr = merge_lib("libmpfr.*.dylib", "mpfr")


# Fix dependencies
def fix_gmp_dependency(lib: Path) -> None:
    cmd = f"otool -L {lib}"
    print(f"{GREEN}{cmd}{END}")
    cmd_output = subprocess.check_output(cmd, shell=True, text=True).split()
    # Find libgmp* but not libgmpxx*
    cmd_output = [value for value in cmd_output if (("libgmp" in value) and (str(lib) not in value))]
    # Remove duplicates from cmd_output
    bad_gmp_dependencies = list(dict.fromkeys(cmd_output))

    for bad_gmp_dependency in bad_gmp_dependencies:
        cmd = f"install_name_tool -change {bad_gmp_dependency} {cwd}/{gmp} {lib}"
        print(f"{GREEN}{cmd}{END}")
        subprocess.check_call(cmd, shell=True)


fix_gmp_dependency(mpfr)
fix_gmp_dependency(gmpxx)

# Symlinks
for lib in [gmp, gmpxx, mpfr]:
    name = str(lib).split(".")[0] + ".dylib"
    if not Path(name).exists():
        print(f"{GREEN}ln -s {str(lib)} {name}{END}")
        os.symlink(lib, name)
