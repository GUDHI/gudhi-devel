import os
import subprocess
import glob
from pathlib import Path

# Oldest macOS version supported by brew bottles
OSX_VERSION = "ventura"

# Some color to distinguish outputs
GREEN = '\033[92m'
END = '\033[0m'

def fetch_and_extract(arch: str, pkg: str):
    """Fetch a brew bottle for given arch and package, then extract it.
    Assumes that the user has enough rights to run brew fetch.
    Parameters:
        arch (str): Architecture ('x86_64' or 'arm64').
        pkg (str): Brew package name.
    """
    if not arch in ['x86_64', 'arm64']:
        raise RuntimeError(f"Unknown architecture {arch}")
    cmd = f"brew fetch --bottle-tag={arch}_{OSX_VERSION} {pkg}"
    subprocess.check_call(cmd, shell=True)
    # Get tar name from brew fetch message
    cmd = f"brew --cache --bottle-tag={arch}_{OSX_VERSION} {pkg}"
    tarball = subprocess.check_output(cmd, shell=True, text=True)
    print(f"tarball ={tarball}")
    # Here we should only find "bottle.tar.gz" in the command output
    assert "bottle.tar.gz" in tarball
    print(f"extract tarball {tarball} in {Path.cwd()}")
    subprocess.check_call(f"tar xf {tarball}", shell=True)

cwd = Path.cwd()

# amd64 deps
Path("deps-amd64").mkdir(exist_ok=True)
os.chdir("deps-amd64")
fetch_and_extract("x86_64", "gmp")
fetch_and_extract("x86_64", "mpfr")
os.chdir(cwd)

# arm64 deps
Path("deps-arm64").mkdir(exist_ok=True)
os.chdir("deps-arm64")
fetch_and_extract("arm64", "gmp")
fetch_and_extract("arm64", "mpfr")
os.chdir(cwd)

# Merge libraries
uni_lib = Path("deps-uni/lib")
uni_lib.mkdir(parents=True, exist_ok=True)

def merge_lib(pattern, pkg, outname=None):
    amd64_lib = glob.glob(f"deps-amd64/{pkg}/*/lib/{pattern}")[0]
    arm64_lib = glob.glob(f"deps-arm64/{pkg}/*/lib/{pattern}")[0]
    out = uni_lib / os.path.basename(amd64_lib)
    cmd = f"lipo -create {amd64_lib} {arm64_lib} -output {out}"
    print(GREEN + cmd + END)
    subprocess.check_call(cmd, shell=True)
    cmd = f"install_name_tool -id {cwd}/{out} {out}"
    print(GREEN + cmd + END)
    subprocess.check_call(cmd, shell=True)
    return out

gmp = merge_lib("libgmp.*.dylib", "gmp")
gmpxx = merge_lib("libgmpxx.*.dylib", "gmp")
mpfr = merge_lib("libmpfr.*.dylib", "mpfr")

# Fix dependencies
def fix_gmp_dependency(lib):
    cmd = f"otool -L {lib}"
    print(GREEN + cmd + END)
    cmd_output = subprocess.check_output(cmd, shell=True, text=True).split()
    # Find libgmp* but not libgmpxx*
    cmd_output = [value for value in cmd_output if (("libgmp" in value) and not (str(lib) in value))]
    # Remove duplicates from cmd_output
    bad_gmp_dependencies = list(dict.fromkeys(cmd_output))

    for bad_gmp_dependency in bad_gmp_dependencies:
        cmd = f"install_name_tool -change {bad_gmp_dependency} {cwd}/{gmp} {lib}"
        print(GREEN + cmd + END)
        subprocess.check_call(cmd, shell=True)

fix_gmp_dependency(mpfr)
fix_gmp_dependency(gmpxx)

# Symlinks
for lib in [gmp, gmpxx, mpfr]:
    name = str(lib).split(".")[0] + ".dylib"
    if not Path(name).exists():
        print(f"{GREEN}ln -s {str(lib)} {name}{END}")
        os.symlink(lib, name)
