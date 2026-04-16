# pyi_rth_env.py
import os
import sys
from pathlib import Path

root = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))

bin_dir = root / "bin"
lib_dir = root / "lib"
etc_dir = root / "etc" / "apptainer"
libexec_dir = root / "libexec" / "apptainer"
conf_file = etc_dir / "apptainer.conf"

# Create runtime state directories Apptainer expects.
(root / "var" / "apptainer" / "mnt" / "session").mkdir(parents=True, exist_ok=True)

# Make bundled executables discoverable by name.
old_path = os.environ.get("PATH", "")
path_parts = [str(bin_dir)]
if old_path:
    path_parts.append(old_path)
os.environ["PATH"] = os.pathsep.join(path_parts)

# Make bundled shared libraries discoverable.
# Put extraction root first for Qt/OpenGL.
# Put lib/ second for helper tools like squashfuse that may expect ../lib.
old_ld = os.environ.get("LD_LIBRARY_PATH", "")
ld_parts = [str(root), str(lib_dir)]
if old_ld:
    ld_parts.append(old_ld)
os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(ld_parts)

# Point Apptainer explicitly at bundled config and libexec trees.
if conf_file.exists():
    os.environ["APPTAINER_CONFIG_FILE"] = str(conf_file)

if etc_dir.exists():
    os.environ["APPTAINER_CONFDIR"] = str(etc_dir)

if libexec_dir.exists():
    os.environ["APPTAINER_LIBEXECDIR"] = str(libexec_dir)
