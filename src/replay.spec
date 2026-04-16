# REPLAY.spec
# Build with:
#   conda activate your_env
#   pyinstaller --noconfirm --distpath . replay.spec

import os
import subprocess
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files

APP_NAME = "REPLAY"
ENTRY_SCRIPT = "replay.py"

CONDA_PREFIX = os.environ.get("CONDA_PREFIX", "")
if not CONDA_PREFIX:
    raise RuntimeError("CONDA_PREFIX is not set. Activate the conda environment first.")

CONDA_PREFIX = Path(CONDA_PREFIX)
CONDA_BIN = CONDA_PREFIX / "bin"
CONDA_LIB = CONDA_PREFIX / "lib"
CONDA_ETC = CONDA_PREFIX / "etc"
CONDA_LIBEXEC = CONDA_PREFIX / "libexec"

APPTAINER = CONDA_BIN / "apptainer"

# Keep PySide6 data files; let standard PyInstaller hooks collect Python modules.
datas = collect_data_files("PySide6")
binaries = []


def add_binary_if_exists(src, dest_dir="."):
    src = Path(src)
    if src.exists() and src.is_file():
        binaries.append((str(src), dest_dir))
        print(f"[spec] Added binary: {src} -> {dest_dir}")
        return True
    print(f"[spec] Missing binary, skipped: {src}")
    return False


def add_data_tree_if_exists(src_dir, dest_root):
    """
    Recursively add all files under src_dir into datas as:
      (absolute_source_file, relative_destination_dir)
    """
    src_dir = Path(src_dir)
    if not src_dir.exists() or not src_dir.is_dir():
        print(f"[spec] Missing data tree, skipped: {src_dir}")
        return False

    added_any = False
    for path in src_dir.rglob("*"):
        if path.is_file():
            rel_parent = path.parent.relative_to(src_dir)
            if str(rel_parent) == ".":
                dest_dir = dest_root
            else:
                dest_dir = f"{dest_root}/{rel_parent.as_posix()}"
            datas.append((str(path), dest_dir))
            print(f"[spec] Added data file: {path} -> {dest_dir}")
            added_any = True

    if added_any:
        print(f"[spec] Added data tree: {src_dir} -> {dest_root}")
    else:
        print(f"[spec] Data tree exists but contains no files: {src_dir}")

    return added_any


def find_system_library(libname: str):
    try:
        result = subprocess.run(
            ["ldconfig", "-p"],
            capture_output=True,
            text=True,
            check=False,
        )
        for line in result.stdout.splitlines():
            if libname in line and "=>" in line:
                return line.split("=>", 1)[1].strip()
    except Exception:
        pass

    candidates = [
        f"/usr/lib/x86_64-linux-gnu/{libname}",
        f"/usr/lib64/{libname}",
        f"/usr/lib/{libname}",
        f"/lib/x86_64-linux-gnu/{libname}",
        f"/lib64/{libname}",
        f"/lib/{libname}",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return candidate
    return None


# -------------------------
# Required: apptainer binary
# -------------------------
if not APPTAINER.exists():
    raise RuntimeError(
        f"Could not find apptainer at {APPTAINER}. "
        "Install it in the active conda env or edit the spec."
    )

# Place apptainer under bin/ so the runtime hook can prepend it to PATH.
add_binary_if_exists(APPTAINER, "bin")


# --------------------------------------------------
# Extra libraries used by apptainer and helper tools
# Keep these in lib/ because helper binaries like
# squashfuse commonly expect ../lib relative to bin/.
# --------------------------------------------------
helper_conda_libs = [
    CONDA_LIB / "libseccomp.so.2",
    CONDA_LIB / "libz.so.1",
    CONDA_LIB / "liblzma.so.5",
    CONDA_LIB / "liblzo2.so.2",
    CONDA_LIB / "liblz4.so.1",
    CONDA_LIB / "libzstd.so.1",
    CONDA_LIB / "libfuse3.so.4",
    # Add more here if later ldd output shows more conda-resident deps:
    # CONDA_LIB / "libarchive.so.13",
    # CONDA_LIB / "libstdc++.so.6",
]

for lib in helper_conda_libs:
    add_binary_if_exists(lib, "lib")


# -------------------------------
# Optional: OpenGL system library
# Keep this at the extraction root, not lib/,
# because moving it broke PySide import behavior.
# -------------------------------
libopengl = find_system_library("libOpenGL.so.0")
if libopengl:
    add_binary_if_exists(libopengl, ".")
else:
    print("[spec] libOpenGL.so.0 not found; relying on target system OpenGL stack.")


# -----------------------------------------
# Required: Apptainer config tree
# buildcfg:
#   APPTAINER_CONFDIR=/.../etc/apptainer
#   APPTAINER_CONF_FILE=/.../etc/apptainer/apptainer.conf
# -----------------------------------------
possible_apptainer_etc = [
    CONDA_ETC / "apptainer",
    Path("/etc/apptainer"),
    Path("/usr/local/etc/apptainer"),
]

apptainer_etc_added = False
for etc_dir in possible_apptainer_etc:
    if add_data_tree_if_exists(etc_dir, "etc/apptainer"):
        apptainer_etc_added = True
        break

if not apptainer_etc_added:
    raise RuntimeError(
        "Could not find an Apptainer config directory. "
        "Expected one of: "
        f"{CONDA_ETC / 'apptainer'}, /etc/apptainer, /usr/local/etc/apptainer"
    )


# -----------------------------------------
# Required: Apptainer libexec tree
# buildcfg:
#   LIBEXECDIR=/.../libexec
#   PLUGIN_ROOTDIR=/.../libexec/apptainer/plugin
# -----------------------------------------
possible_apptainer_libexec = [
    CONDA_LIBEXEC / "apptainer",
    Path("/usr/libexec/apptainer"),
    Path("/usr/local/libexec/apptainer"),
]

apptainer_libexec_added = False
for libexec_dir in possible_apptainer_libexec:
    if add_data_tree_if_exists(libexec_dir, "libexec/apptainer"):
        apptainer_libexec_added = True
        break

if not apptainer_libexec_added:
    raise RuntimeError(
        "Could not find an Apptainer libexec directory. "
        "Expected one of: "
        f"{CONDA_LIBEXEC / 'apptainer'}, /usr/libexec/apptainer, /usr/local/libexec/apptainer"
    )


# -----------------------------------------
# Optional helper executables
# These were mentioned in Apptainer runtime messages.
# -----------------------------------------
helper_candidates = [
    CONDA_BIN / "squashfuse",
    CONDA_BIN / "fuse2fs",
    CONDA_BIN / "gocryptfs",
    Path("/usr/bin/squashfuse"),
    Path("/usr/bin/fuse2fs"),
    Path("/usr/bin/gocryptfs"),
]

seen_helpers = set()
for helper in helper_candidates:
    helper = Path(helper)
    if helper.exists() and helper.is_file() and helper.name not in seen_helpers:
        add_binary_if_exists(helper, "bin")
        seen_helpers.add(helper.name)


a = Analysis(
    [ENTRY_SCRIPT],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=["pyi_rth_env.py"],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name=APP_NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    onefile=True,
)
