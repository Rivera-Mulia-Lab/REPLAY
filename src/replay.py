import gzip
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml
from PySide6.QtCore import Qt, QProcess, QTimer
from PySide6.QtGui import QAction, QActionGroup, QColor, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)
from pathlib import Path
import os
import sys
import subprocess


try:
    import qdarktheme
except Exception:
    qdarktheme = None

CANONICAL_RE = re.compile(r"^(?P<desc>[^_]+)_(?P<sample>\d+)_(?P<el>[EL])_(?P<read>R[12])\.fastq\.gz$")
SAMPLE_TOKEN_RE = re.compile(r"^(sample|rep|repeat)(\d+)$", re.IGNORECASE)
PLAIN_NUMERIC_RE = re.compile(r"^(\d+)$")
TOKEN_SPLIT_RE = re.compile(r"[_\-\s\.]+")
PROGRESS_RE = re.compile(r"(?P<done>\d+) of (?P<total>\d+) steps \((?P<pct>\d+)%\) done")
NOTHING_TO_DO_TEXT = "Nothing to be done (all requested files are present and up to date)."

PAIR_FWD_DEFAULT = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
PAIR_REV_DEFAULT = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
SINGLE_FWD_DEFAULT = "AGATCGGAAGAGCACACGTCTG"

PRESET_HG38 = "Human (hg38) Default"
PRESET_MM39 = "Mouse (mm39) Default"
PRESET_CUSTOM = "Custom Configuration"

NORMALIZATION_DISPLAY_TO_VALUE = {
    "Quantile": "quantile",
    "IQR": "iqr",
    "Median": "median",
    "None": "none",
}
NORMALIZATION_VALUE_TO_DISPLAY = {v: k for k, v in NORMALIZATION_DISPLAY_TO_VALUE.items()}
SMOOTHING_DISPLAY_TO_VALUE = {
    "None": "none",
    "Gaussian": "gaussian",
    "Loess": "loess",
}
SMOOTHING_VALUE_TO_DISPLAY = {v: k for k, v in SMOOTHING_DISPLAY_TO_VALUE.items()}

DEFAULT_CONFIG = {
    "latency-wait": 120,
    "genome": "",
    "output_folder": "",
    "samples_folder": "",
    "window_sizes": [5000, 20000],
    "reads": ["R1"],
    "barcode1": SINGLE_FWD_DEFAULT,
    "barcode2": "",
    "contigs": [],
    "rt_processing": {
        "normalization": "quantile",
        "smoothing": "loess",
        "median_center": True,
        "target_file": "resources/TargetDataset_RT.txt",
        "gap_file": "resources/hg38_gaplist.tsv",
        "mask_gaps": True,
        "loess_span_bp": 500000,
        "clip_min": -8,
        "clip_max": 8,
    },
    "target": [],
}

ROW_OK_COLOR = QColor(230, 255, 230)
ROW_WARN_COLOR = QColor(255, 245, 210)
ROW_ERROR_COLOR = QColor(255, 230, 230)
ROW_EXCLUDED_COLOR = QColor(235, 235, 235)


def canonical_validate_filename(filename: str) -> Tuple[bool, str, Optional[Dict[str, str]]]:
    m = CANONICAL_RE.match(filename)
    if not m:
        return False, "Does not match canonical pattern desc_num_E/L_R1/R2.fastq.gz", None
    return True, "Valid canonical filename", m.groupdict()


def normalize_el_token(token: str) -> Optional[str]:
    token_l = token.lower()
    if token_l in {"e", "early"}:
        return "E"
    if token_l in {"l", "late"}:
        return "L"
    return None


def normalize_read_token(token: str) -> Optional[str]:
    token_l = token.lower()
    if token_l in {"r1", "read1", "r_1"}:
        return "R1"
    if token_l in {"r2", "read2", "r_2"}:
        return "R2"
    return None


def extract_sample_number(token: str) -> Optional[str]:
    m = SAMPLE_TOKEN_RE.match(token)
    if m:
        return m.group(2)
    m2 = PLAIN_NUMERIC_RE.match(token)
    if m2:
        return m2.group(1)
    return None


def strip_fastq_suffix(name: str) -> str:
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def attempted_parse_and_propose(filename: str) -> Dict[str, str]:
    result = {"status": "failed", "original": filename, "proposed": "", "reason": ""}
    valid, _, _ = canonical_validate_filename(filename)
    if valid:
        result["status"] = "valid"
        result["proposed"] = filename
        result["reason"] = "Already canonical"
        return result

    stem = strip_fastq_suffix(filename)
    tokens = [t for t in TOKEN_SPLIT_RE.split(stem) if t]
    el_hits, read_hits, sample_hits = [], [], []
    used_idx = set()

    for i, tok in enumerate(tokens):
        norm_el = normalize_el_token(tok)
        if norm_el is not None:
            el_hits.append((i, norm_el, tok))
        norm_read = normalize_read_token(tok)
        if norm_read is not None:
            read_hits.append((i, norm_read, tok))
        samp = extract_sample_number(tok)
        if samp is not None:
            sample_hits.append((i, samp, tok))

    if len(el_hits) != 1:
        result["status"] = "ambiguous" if len(el_hits) > 1 else "failed"
        result["reason"] = f"Expected exactly one E/L token, found {len(el_hits)}"
        return result
    if len(read_hits) != 1:
        result["status"] = "ambiguous" if len(read_hits) > 1 else "failed"
        result["reason"] = f"Expected exactly one read token, found {len(read_hits)}"
        return result
    if len(sample_hits) != 1:
        result["status"] = "ambiguous" if len(sample_hits) > 1 else "failed"
        result["reason"] = f"Expected exactly one sample token, found {len(sample_hits)}"
        return result

    el_idx, el_val, _ = el_hits[0]
    read_idx, read_val, _ = read_hits[0]
    sample_idx, sample_val, _ = sample_hits[0]
    used_idx.update([el_idx, read_idx, sample_idx])

    desc_tokens = [tok for i, tok in enumerate(tokens) if i not in used_idx]
    if not desc_tokens:
        result["status"] = "failed"
        result["reason"] = "No remaining tokens available for description"
        return result

    desc_raw = "-".join(desc_tokens)
    desc_clean = re.sub(r"[^A-Za-z0-9\-]+", "", desc_raw)
    if not desc_clean:
        result["status"] = "failed"
        result["reason"] = "Description became empty after sanitization"
        return result

    result["status"] = "proposed"
    result["proposed"] = f"{desc_clean}_{sample_val}_{el_val}_{read_val}.fastq.gz"
    result["reason"] = "Proposed rename from loose parse"
    return result


def dump_yaml_text(config: Dict) -> str:
    return yaml.safe_dump(config, sort_keys=False, default_flow_style=False)


def load_yaml_text(text: str) -> Dict:
    data = yaml.safe_load(text)
    return data if data is not None else {}


def ensure_parent(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)


def genome_id_from_path(genome_text: str) -> str:
    name = Path(genome_text).name
    if name.endswith(".fa.gz"):
        return name[:-6]
    if name.endswith(".fa"):
        return name[:-3]
    if name.endswith(".fasta"):
        return name[:-6]
    return Path(name).stem


def dir_with_trailing_slash(path: Path) -> str:
    return str(path) + os.sep


def get_total_system_memory_mb() -> int:
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        phys_pages = os.sysconf("SC_PHYS_PAGES")
        return max(1, int((page_size * phys_pages) / 1_000_000))
    except Exception:
        pass
    try:
        import psutil  # type: ignore
        return max(1, int(psutil.virtual_memory().total / 1_000_000))
    except Exception:
        return 0


def get_total_system_cores() -> int:
    return max(1, os.cpu_count() or 1)


class CollapsibleSection(QWidget):
    def __init__(self, title: str, content_widget: QWidget, expanded: bool = True, parent=None):
        super().__init__(parent)
        self.toggle_button = QToolButton(text=title, checkable=True, checked=expanded)
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(Qt.DownArrow if expanded else Qt.RightArrow)
        self.toggle_button.clicked.connect(self._on_toggled)
        self.content_widget = content_widget
        self.content_widget.setVisible(expanded)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toggle_button)
        layout.addWidget(self.content_widget)

    def _on_toggled(self, checked: bool):
        self.toggle_button.setArrowType(Qt.DownArrow if checked else Qt.RightArrow)
        self.content_widget.setVisible(checked)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("REPLAY")
        self.resize(1840, 1040)

        self.current_config = DEFAULT_CONFIG.copy()
        self.rename_plan: List[Dict[str, str]] = []
        self.local_run_started_at: float = 0.0
        self.run_reached_100: bool = False
        self.pending_unlock: bool = False
        self.run_process: Optional[QProcess] = None
        self.current_log_path: Optional[Path] = None
        self.current_log_pos: int = 0
        self.current_profile_mode: str = PRESET_CUSTOM
        self.dark_mode_enabled: bool = False
        self.results_index: int = -1
        self.total_system_memory_mb: int = get_total_system_memory_mb()
        self.total_system_cores: int = get_total_system_cores()
        self.base_font_point_size: float = max(8.0, QApplication.instance().font().pointSizeF() if QApplication.instance() else 9.0)
        self.current_font_size: str = "Small"

        self.log_tail_timer = QTimer(self)
        self.log_tail_timer.setInterval(1000)
        self.log_tail_timer.timeout.connect(self.poll_latest_snakemake_log)

        self._build_ui()
        self.populate_form_from_config(DEFAULT_CONFIG)
        self.set_default_workflow_dir()
        self.refresh_preset_availability()
        if self.hg38_available():
            self.combo_preset.setCurrentText(PRESET_HG38)
            self.apply_default_hg38()
        else:
            self.combo_preset.setCurrentText(PRESET_CUSTOM)
        self.refresh_yaml_preview()
        self.update_barcode_defaults(force_if_blank=True)
        self.update_file_check_status()
        self.update_file_action_buttons()
        self.apply_theme(False)
        self.apply_font_size("Small")
        self.update_resource_limits_status()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self):
        self._build_menu()
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        self.header_widget = QWidget()
        header_layout = QHBoxLayout(self.header_widget)
        header_layout.setContentsMargins(4, 4, 4, 4)
        header_layout.setSpacing(10)
        self.logo_label = QLabel()
        self.logo_label.setMinimumWidth(180)
        self.logo_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        header_layout.addWidget(self.logo_label)
        header_text = QVBoxLayout()
        header_text.setContentsMargins(0, 0, 0, 0)
        header_text.setSpacing(0)
        self.header_text_label = QLabel(
            'A reproducible and user-friendly application for DNA replication timing analysis<br>'
            '<a href="https://github.com/Rivera-Mulia-Lab/REPLAY">https://github.com/Rivera-Mulia-Lab/REPLAY</a>'
        )
        self.header_text_label.setTextFormat(Qt.RichText)
        self.header_text_label.setOpenExternalLinks(True)
        self.header_text_label.setWordWrap(True)

        header_text.addWidget(self.header_text_label)
        header_layout.addLayout(header_text, 1)
        root.addWidget(self.header_widget)

        self.main_splitter = QSplitter(Qt.Horizontal)
        root.addWidget(self.main_splitter)
        root.setStretch(0, 15)
        root.setStretch(1, 85)

        self.tabs = QTabWidget()
        self.tabs.currentChanged.connect(self.on_tab_changed)
        self.main_splitter.addWidget(self.tabs)

        self.right_panel = QWidget()
        right_layout = QVBoxLayout(self.right_panel)
        self.main_splitter.addWidget(self.right_panel)
        self.main_splitter.setSizes([1200, 620])

        self._build_files_tab()
        self._build_config_tab()
        self._build_adv_tab()
        self._build_results_tab()

        self.progress_container = QWidget()
        progress_layout = QVBoxLayout(self.progress_container)
        progress_layout.setContentsMargins(0, 0, 0, 0)
        self.progress_label = QLabel("No Snakemake progress detected yet")
        progress_layout.addWidget(self.progress_label)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)
        right_layout.addWidget(self.progress_container)

        yaml_container = QWidget()
        yaml_layout = QVBoxLayout(yaml_container)
        yaml_layout.setContentsMargins(0, 0, 0, 0)
        self.yaml_preview = QPlainTextEdit()
        self.yaml_preview.setReadOnly(False)
        yaml_layout.addWidget(self.yaml_preview, stretch=1)
        yaml_btn_row = QHBoxLayout()
        self.btn_refresh_yaml = QPushButton("Refresh YAML Preview")
        self.btn_refresh_yaml.clicked.connect(self.refresh_yaml_preview)
        yaml_btn_row.addWidget(self.btn_refresh_yaml)
        self.btn_import_preview = QPushButton("Import Preview From File")
        self.btn_import_preview.clicked.connect(self.load_yaml_file)
        yaml_btn_row.addWidget(self.btn_import_preview)
        self.btn_export_preview = QPushButton("Export Preview To File")
        self.btn_export_preview.clicked.connect(self.save_yaml_file)
        yaml_btn_row.addWidget(self.btn_export_preview)
        yaml_layout.addLayout(yaml_btn_row)
        self.yaml_section = CollapsibleSection("YAML Preview", yaml_container, expanded=False)
        right_layout.addWidget(self.yaml_section, stretch=3)

        log_container = QWidget()
        log_layout = QVBoxLayout(log_container)
        log_layout.setContentsMargins(0, 0, 0, 0)
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        log_layout.addWidget(self.log_box, stretch=1)
        self.log_section = CollapsibleSection("Application Log", log_container, expanded=False)
        right_layout.addWidget(self.log_section, stretch=2)

        self.update_header_height()

    def _build_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        act_load = QAction("Load YAML...", self)
        act_load.triggered.connect(self.load_yaml_file)
        file_menu.addAction(act_load)
        act_save = QAction("Save YAML...", self)
        act_save.triggered.connect(self.save_yaml_file)
        file_menu.addAction(act_save)
        act_exit = QAction("Exit", self)
        act_exit.triggered.connect(self.close)
        file_menu.addAction(act_exit)

        options_menu = menubar.addMenu("Options")
        self.act_dark_mode = QAction("Dark Mode", self)
        self.act_dark_mode.setCheckable(True)
        self.act_dark_mode.toggled.connect(self.apply_theme)
        options_menu.addAction(self.act_dark_mode)

        font_menu = options_menu.addMenu("Font Size")
        self.font_action_group = QActionGroup(self)
        self.font_action_group.setExclusive(True)
        self.font_actions = {}
        for label in ["Small", "Medium", "Large"]:
            act = QAction(label, self)
            act.setCheckable(True)
            if label == "Small":
                act.setChecked(True)
            act.triggered.connect(lambda checked, size_label=label: self.apply_font_size(size_label))
            self.font_action_group.addAction(act)
            self.font_actions[label] = act
            font_menu.addAction(act)

    def _make_prefixed_path_row(self, prefix_text: str, line_edit: QLineEdit, browse_callback=None):
        row = QHBoxLayout()
        row.addWidget(QLabel(prefix_text))
        row.addWidget(line_edit)
        if browse_callback is not None:
            btn = QPushButton("Browse")
            btn.clicked.connect(browse_callback)
            row.addWidget(btn)
        wrap = QWidget()
        wrap.setLayout(row)
        return wrap

    def _make_browse_row(self, line_edit: QLineEdit, browse_callback):
        row = QHBoxLayout()
        row.addWidget(line_edit)
        btn = QPushButton("Browse")
        btn.clicked.connect(browse_callback)
        row.addWidget(btn)
        wrap = QWidget()
        wrap.setLayout(row)
        return wrap

    def _build_files_tab(self):
        tab = QWidget()
        self.tabs.addTab(tab, "Files")
        outer = QVBoxLayout(tab)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll)
        host = QWidget()
        scroll.setWidget(host)
        layout = QVBoxLayout(host)

        input_box = QGroupBox("Input folder and sample staging")
        input_form = QFormLayout(input_box)
        self.edit_input_dir = QLineEdit()
        self.edit_input_dir.textChanged.connect(self.on_input_dir_changed)
        input_form.addRow("Raw Reads Folder", self._make_browse_row(self.edit_input_dir, self.choose_input_dir))

        self.edit_files_samples_folder = QLineEdit()
        self.edit_files_samples_folder.textChanged.connect(self.sync_samples_folder_from_files_tab)
        input_form.addRow("Pipeline Reads Folder", self._make_browse_row(self.edit_files_samples_folder, self.choose_pipeline_reads_dir))

        self.edit_output_folder = QLineEdit()
        self.edit_output_folder.textChanged.connect(self.update_file_check_status)
        input_form.addRow("Output Folder", self._make_browse_row(self.edit_output_folder, self.choose_output_folder))
        layout.addWidget(CollapsibleSection("Input folder and sample staging", input_box, expanded=True))

        read_box = QGroupBox("Read Options")
        read_form = QFormLayout(read_box)
        mode_row = QHBoxLayout()
        self.check_single = QCheckBox("Single-Read")
        self.check_paired = QCheckBox("Paired-End")
        self.check_single.stateChanged.connect(self.on_read_mode_changed)
        self.check_paired.stateChanged.connect(self.on_read_mode_changed)
        mode_row.addWidget(self.check_single)
        mode_row.addWidget(self.check_paired)
        mode_wrap = QWidget()
        mode_wrap.setLayout(mode_row)
        read_form.addRow("Sequencing Format", mode_wrap)
        self.edit_barcode1 = QLineEdit()
        self.edit_barcode1.textChanged.connect(self.refresh_yaml_preview)
        read_form.addRow("Forward Barcode", self.edit_barcode1)
        self.edit_barcode2 = QLineEdit()
        self.edit_barcode2.textChanged.connect(self.refresh_yaml_preview)
        read_form.addRow("Reverse Barcode", self.edit_barcode2)
        layout.addWidget(CollapsibleSection("Read Options", read_box, expanded=True))

        config_box = QGroupBox("Configuration")
        config_form = QFormLayout(config_box)
        self.combo_preset = QComboBox()
        self.combo_preset.addItems([PRESET_HG38, PRESET_MM39, PRESET_CUSTOM])
        self.combo_preset.currentTextChanged.connect(self.on_preset_changed)
        config_form.addRow("Preset", self.combo_preset)
        dl_row = QHBoxLayout()
        self.btn_download_hg38 = QPushButton("Download hg38")
        self.btn_download_hg38.clicked.connect(lambda: self.download_genome("hg38"))
        dl_row.addWidget(self.btn_download_hg38)
        self.btn_download_mm39 = QPushButton("Download mm39")
        self.btn_download_mm39.clicked.connect(lambda: self.download_genome("mm39"))
        dl_row.addWidget(self.btn_download_mm39)
        dl_wrap = QWidget()
        dl_wrap.setLayout(dl_row)
        config_form.addRow("Genome Downloads", dl_wrap)
        self.combo_normalization = QComboBox()
        self.combo_normalization.addItems(list(NORMALIZATION_DISPLAY_TO_VALUE.keys()))
        self.combo_normalization.currentIndexChanged.connect(self.update_file_check_status)
        config_form.addRow("Normalization", self.combo_normalization)
        self.combo_smoothing = QComboBox()
        self.combo_smoothing.addItems(list(SMOOTHING_DISPLAY_TO_VALUE.keys()))
        self.combo_smoothing.setCurrentText("Loess")
        self.combo_smoothing.currentIndexChanged.connect(self.refresh_yaml_preview)
        config_form.addRow("Smoothing", self.combo_smoothing)

        mem_row = QHBoxLayout()
        self.spin_runtime_memory = QSpinBox()
        self.spin_runtime_memory.setRange(1, 1_000_000)
        self.spin_runtime_memory.setValue(12)
        self.spin_runtime_memory.valueChanged.connect(self.update_resource_limits_status)
        mem_row.addWidget(self.spin_runtime_memory)
        self.combo_runtime_memory_unit = QComboBox()
        self.combo_runtime_memory_unit.addItems(["MB", "GB"])
        self.combo_runtime_memory_unit.setCurrentText("GB")
        self.combo_runtime_memory_unit.currentIndexChanged.connect(self.update_resource_limits_status)
        mem_row.addWidget(self.combo_runtime_memory_unit)
        mem_wrap = QWidget()
        mem_wrap.setLayout(mem_row)
        config_form.addRow("Memory", mem_wrap)

        core_row = QHBoxLayout()
        self.spin_runtime_cores = QSpinBox()
        self.spin_runtime_cores.setRange(1, max(1, self.total_system_cores if self.total_system_cores else 1024))
        self.spin_runtime_cores.setValue(2)
        self.spin_runtime_cores.valueChanged.connect(self.update_resource_limits_status)
        core_row.addWidget(self.spin_runtime_cores)
        core_wrap = QWidget()
        core_wrap.setLayout(core_row)
        config_form.addRow("Cores", core_wrap)

        self.btn_runtime_max = QPushButton("Max")
        self.btn_runtime_max.clicked.connect(self.apply_max_runtime_resources)
        config_form.addRow("System Max", self.btn_runtime_max)

        layout.addWidget(CollapsibleSection("Configuration", config_box, expanded=True))

        action_box = QGroupBox("File actions")
        action_layout = QVBoxLayout(action_box)
        top_row = QHBoxLayout()
        self.btn_scan = QPushButton("Scan + Auto-populate")
        self.btn_scan.clicked.connect(self.scan_input_folder)
        top_row.addWidget(self.btn_scan)
        self.btn_copy_or_rename = QPushButton("Copy")
        self.btn_copy_or_rename.clicked.connect(self.apply_primary_file_action)
        top_row.addWidget(self.btn_copy_or_rename)
        self.btn_symlink = QPushButton("Symlink")
        self.btn_symlink.clicked.connect(lambda: self.apply_rename_plan(mode="symlink"))
        top_row.addWidget(self.btn_symlink)
        action_layout.addLayout(top_row)
        overwrite_row = QHBoxLayout()
        self.check_overwrite = QCheckBox("Allow overwrite on copy/rename/symlink")
        self.check_overwrite.stateChanged.connect(self.update_file_action_buttons)
        overwrite_row.addWidget(self.check_overwrite)
        self.label_overwrite_warning = QLabel("")
        self.label_overwrite_warning.setStyleSheet("color: red;")
        overwrite_row.addWidget(self.label_overwrite_warning)
        overwrite_row.addStretch(1)
        action_layout.addLayout(overwrite_row)
        layout.addWidget(CollapsibleSection("File actions", action_box, expanded=True))

        table_box = QGroupBox("Files Table")
        table_layout = QVBoxLayout(table_box)
        self.table_files = QTableWidget(0, 10)
        self.table_files.setHorizontalHeaderLabels([
            "Original", "Include", "Description", "Replicate", "E/L", "Reads", "Proposed", "Successful Copy/Link", "Row Check", "Group Check"
        ])
        self.table_files.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table_layout.addWidget(self.table_files)
        layout.addWidget(CollapsibleSection("Files Table", table_box, expanded=True))

        run_box = QGroupBox("Local run")
        run_layout = QVBoxLayout(run_box)
        self.label_file_check_status = QLabel("File checks not yet run")
        run_layout.addWidget(self.label_file_check_status)
        self.label_action_status = QLabel("Validation state unknown")
        run_layout.addWidget(self.label_action_status)
        self.label_resource_error = QLabel("")
        self.label_resource_error.setWordWrap(True)
        run_layout.addWidget(self.label_resource_error)
        self.btn_check_files = QPushButton("Check Required Files")
        self.btn_check_files.clicked.connect(self.update_file_check_status)
        run_layout.addWidget(self.btn_check_files)
        self.btn_run_local = QPushButton("Run Local")
        self.btn_run_local.clicked.connect(self.run_local)
        run_layout.addWidget(self.btn_run_local)
        layout.addWidget(CollapsibleSection("Local run", run_box, expanded=True))
        layout.addStretch(1)

    def _build_config_tab(self):
        tab = QWidget()
        self.tabs.addTab(tab, "Configuration")
        self.config_tab_index = self.tabs.indexOf(tab)
        outer = QVBoxLayout(tab)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll)
        form_host = QWidget()
        scroll.setWidget(form_host)
        layout = QVBoxLayout(form_host)

        general_box = QGroupBox("Custom Configuration")
        general_form = QFormLayout(general_box)
        self.edit_workflow_dir = QLineEdit()
        self.edit_workflow_dir.textChanged.connect(self.on_workflow_dir_changed)
        workflow_row = QHBoxLayout()
        workflow_row.addWidget(self.edit_workflow_dir)
        self.btn_browse_workflow_dir = QPushButton("Browse")
        self.btn_browse_workflow_dir.clicked.connect(self.choose_workflow_dir)
        workflow_row.addWidget(self.btn_browse_workflow_dir)
        workflow_wrap = QWidget()
        workflow_wrap.setLayout(workflow_row)
        general_form.addRow("Workflow directory", workflow_wrap)

        self.edit_genome = QLineEdit()
        self.edit_genome.textChanged.connect(self.on_genome_changed)
        self.genome_prefix = self._make_prefixed_path_row("(workflow_dir)/", self.edit_genome, self.choose_reference_genome)
        general_form.addRow("Reference Genome", self.genome_prefix)

        self.edit_target_file = QLineEdit("resources/TargetDataset_RT.txt")
        self.edit_target_file.textChanged.connect(self.update_file_check_status)
        self.target_prefix = self._make_prefixed_path_row("(workflow_dir)/", self.edit_target_file, self.choose_target_file)
        general_form.addRow("target_file", self.target_prefix)

        self.edit_gap_file = QLineEdit("resources/hg38_gaplist.tsv")
        self.edit_gap_file.textChanged.connect(self.update_file_check_status)
        self.gap_prefix = self._make_prefixed_path_row("(workflow_dir)/", self.edit_gap_file, self.choose_gap_file)
        general_form.addRow("gap_file", self.gap_prefix)

        self.check_median_center = QCheckBox()
        general_form.addRow("median_center", self.check_median_center)
        self.check_mask_gaps = QCheckBox()
        self.check_mask_gaps.stateChanged.connect(self.update_file_check_status)
        general_form.addRow("mask_gaps", self.check_mask_gaps)
        layout.addWidget(general_box)

        window_box = QGroupBox("Window sizes")
        window_layout = QGridLayout(window_box)
        self.list_window_sizes = QListWidget()
        self.list_window_sizes.setSelectionMode(QAbstractItemView.MultiSelection)
        for val in [5000, 10000, 20000, 50000, 100000]:
            item = QListWidgetItem(str(val))
            self.list_window_sizes.addItem(item)
            if val in {5000, 20000}:
                item.setSelected(True)
        window_layout.addWidget(QLabel("window_sizes"), 0, 0)
        window_layout.addWidget(self.list_window_sizes, 1, 0)
        win_btns = QVBoxLayout()
        self.edit_new_window_size = QLineEdit()
        self.edit_new_window_size.setPlaceholderText("Add custom size")
        win_btns.addWidget(self.edit_new_window_size)
        self.btn_add_window = QPushButton("Add window size")
        self.btn_add_window.clicked.connect(self.add_window_size)
        win_btns.addWidget(self.btn_add_window)
        self.btn_remove_window = QPushButton("Remove selected")
        self.btn_remove_window.clicked.connect(self.remove_selected_window_sizes)
        win_btns.addWidget(self.btn_remove_window)
        win_wrap = QWidget()
        win_wrap.setLayout(win_btns)
        window_layout.addWidget(win_wrap, 1, 1)
        layout.addWidget(window_box)

        contig_box = QGroupBox("Contigs to Align")
        contig_layout = QVBoxLayout(contig_box)
        self.list_contigs = QListWidget()
        self.list_contigs.setSelectionMode(QAbstractItemView.MultiSelection)
        self.list_contigs.itemSelectionChanged.connect(self.refresh_yaml_preview)
        contig_layout.addWidget(self.list_contigs)
        row1 = QHBoxLayout()
        self.btn_make_contigs = QPushButton("Generate Contigs from Genome")
        self.btn_make_contigs.clicked.connect(self.generate_contigs_from_genome)
        row1.addWidget(self.btn_make_contigs)
        self.btn_reload_contigs = QPushButton("Reload Contigs File")
        self.btn_reload_contigs.clicked.connect(lambda: self.populate_contigs_from_expected_file(auto_select_all=True))
        row1.addWidget(self.btn_reload_contigs)
        contig_layout.addLayout(row1)
        row2 = QHBoxLayout()
        self.btn_select_all_contigs = QPushButton("Select All")
        self.btn_select_all_contigs.clicked.connect(self.select_all_contigs)
        row2.addWidget(self.btn_select_all_contigs)
        self.btn_deselect_all_contigs = QPushButton("Deselect All")
        self.btn_deselect_all_contigs.clicked.connect(self.deselect_all_contigs)
        row2.addWidget(self.btn_deselect_all_contigs)
        self.btn_write_selected_contigs = QPushButton("Overwrite contigs file with selection")
        self.btn_write_selected_contigs.clicked.connect(self.write_selected_contigs_to_file)
        row2.addWidget(self.btn_write_selected_contigs)
        contig_layout.addLayout(row2)
        layout.addWidget(contig_box)
        layout.addStretch(1)

    def _build_adv_tab(self):
        tab = QWidget()
        self.tabs.addTab(tab, "Advanced Settings")
        outer = QVBoxLayout(tab)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll)
        form_host = QWidget()
        scroll.setWidget(form_host)
        layout = QVBoxLayout(form_host)
        adv_box = QGroupBox("Advanced settings")
        adv_form = QFormLayout(adv_box)
        self.spin_loess_span = QSpinBox()
        self.spin_loess_span.setRange(0, 10_000_000)
        adv_form.addRow("Smoothing Span bp", self.spin_loess_span)
        self.spin_clip_min = QSpinBox()
        self.spin_clip_min.setRange(-999999, 999999)
        adv_form.addRow("clip_min", self.spin_clip_min)
        self.spin_clip_max = QSpinBox()
        self.spin_clip_max.setRange(-999999, 999999)
        adv_form.addRow("clip_max", self.spin_clip_max)
        self.spin_latency = QSpinBox()
        self.spin_latency.setRange(0, 999999)
        adv_form.addRow("latency-wait", self.spin_latency)
        layout.addWidget(adv_box)
        self.btn_unlock = QPushButton("Unlock Snakemake")
        self.btn_unlock.clicked.connect(self.unlock_snakemake)
        layout.addWidget(self.btn_unlock)
        layout.addStretch(1)

    def _build_results_tab(self):
        self.results_tab = QWidget()
        self.tabs.addTab(self.results_tab, "Results")
        self.results_index = self.tabs.indexOf(self.results_tab)
        self.tabs.setTabEnabled(self.results_index, False)
        layout = QHBoxLayout(self.results_tab)
        self.results_sample_list = QListWidget()
        self.results_sample_list.currentRowChanged.connect(self.update_results_view)
        layout.addWidget(self.results_sample_list, 1)
        right = QWidget()
        right_layout = QVBoxLayout(right)
        layout.addWidget(right, 5)

        self.label_norm_results_dir = QLabel("Normalized Results: ")
        right_layout.addWidget(self.label_norm_results_dir)
        self.label_raw_results_dir = QLabel("Raw Results: ")
        right_layout.addWidget(self.label_raw_results_dir)

        top_controls = QHBoxLayout()
        top_controls.addWidget(QLabel("Bin size"))
        self.combo_results_bin_size = QComboBox()
        self.combo_results_bin_size.currentIndexChanged.connect(self.update_results_view)
        top_controls.addWidget(self.combo_results_bin_size)
        self.toggle_reads = QCheckBox("Reads")
        self.toggle_reads.setChecked(True)
        self.toggle_reads.stateChanged.connect(self.update_results_visibility)
        top_controls.addWidget(self.toggle_reads)
        self.toggle_windows = QCheckBox("Windows")
        self.toggle_windows.setChecked(True)
        self.toggle_windows.stateChanged.connect(self.update_results_visibility)
        top_controls.addWidget(self.toggle_windows)
        self.toggle_acf = QCheckBox("ACF")
        self.toggle_acf.setChecked(True)
        self.toggle_acf.stateChanged.connect(self.update_results_visibility)
        top_controls.addWidget(self.toggle_acf)
        self.toggle_raw = QCheckBox("Raw Log2")
        self.toggle_raw.setChecked(True)
        self.toggle_raw.stateChanged.connect(self.update_results_visibility)
        top_controls.addWidget(self.toggle_raw)
        self.toggle_rawkde = QCheckBox("Raw KDE")
        self.toggle_rawkde.setChecked(True)
        self.toggle_rawkde.stateChanged.connect(self.update_results_visibility)
        top_controls.addWidget(self.toggle_rawkde)
        self.toggle_kde = QCheckBox("KDE")
        self.toggle_kde.setChecked(True)
        self.toggle_kde.stateChanged.connect(self.update_results_visibility)
        top_controls.addWidget(self.toggle_kde)
        top_controls.addStretch(1)
        right_layout.addLayout(top_controls)

        self.results_header = QLabel("No sample selected")
        right_layout.addWidget(self.results_header)
        self.results_scroll = QScrollArea()
        self.results_scroll.setWidgetResizable(True)
        right_layout.addWidget(self.results_scroll, 1)
        self.results_content = QWidget()
        self.results_scroll.setWidget(self.results_content)
        self.results_content_layout = QHBoxLayout(self.results_content)

        self.col_reads = self._make_result_column("Reads")
        self.col_windows = self._make_result_column("Windows")
        self.col_acf = self._make_result_column("ACF")
        self.col_raw = self._make_result_column("Raw Log2")
        self.col_rawkde = self._make_result_column("Raw KDE")
        self.col_kde = self._make_result_column("KDE")
        for col in [self.col_reads, self.col_windows, self.col_acf, self.col_raw, self.col_rawkde, self.col_kde]:
            self.results_content_layout.addWidget(col["box"])
        self.results_content_layout.addStretch(1)

    def _make_result_column(self, title: str):
        box = QGroupBox(title)
        layout = QVBoxLayout(box)
        path_label = QLabel("")
        path_label.setWordWrap(True)
        layout.addWidget(path_label)
        if title in {"Reads", "Windows"}:
            label_e = QLabel("E")
            img_e = QLabel("Missing")
            img_e.setAlignment(Qt.AlignCenter)
            label_l = QLabel("L")
            img_l = QLabel("Missing")
            img_l.setAlignment(Qt.AlignCenter)
            layout.addWidget(label_e)
            layout.addWidget(img_e)
            layout.addWidget(label_l)
            layout.addWidget(img_l)
            return {"box": box, "type": "pair", "e": img_e, "l": img_l, "path": path_label}
        img = QLabel("Missing")
        img.setAlignment(Qt.AlignCenter)
        layout.addWidget(img)
        return {"box": box, "type": "single", "img": img, "path": path_label}

    # ------------------------------------------------------------------
    # Theme / header
    # ------------------------------------------------------------------
    def fallback_dark_stylesheet(self) -> str:
        return """
        QWidget { background-color: #202124; color: #e8eaed; }
        QLineEdit, QTextEdit, QPlainTextEdit, QListWidget, QTableWidget, QComboBox, QSpinBox, QScrollArea {
            background-color: #2b2c2f; color: #e8eaed; border: 1px solid #5f6368;
        }
        QPushButton { background-color: #3c4043; color: #e8eaed; border: 1px solid #5f6368; padding: 4px; }
        QPushButton:disabled { color: #9aa0a6; }
        QGroupBox { border: 1px solid #5f6368; margin-top: 8px; }
        QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 3px 0 3px; }
        QHeaderView::section { background-color: #3c4043; color: #e8eaed; }
        """

    def apply_theme(self, enabled: bool):
        self.dark_mode_enabled = bool(enabled)
        app = QApplication.instance()
        applied = False
        if self.dark_mode_enabled:
            if qdarktheme is not None:
                try:
                    if hasattr(qdarktheme, "setup_theme"):
                        qdarktheme.setup_theme("dark")
                    elif hasattr(qdarktheme, "load_stylesheet"):
                        try:
                            app.setStyleSheet(qdarktheme.load_stylesheet("dark"))
                        except TypeError:
                            app.setStyleSheet(qdarktheme.load_stylesheet(theme="dark"))
                    applied = True
                except Exception:
                    applied = False
            if not applied:
                app.setStyleSheet(self.fallback_dark_stylesheet())
        else:
            if qdarktheme is not None:
                try:
                    if hasattr(qdarktheme, "setup_theme"):
                        qdarktheme.setup_theme("light")
                        applied = True
                    elif hasattr(qdarktheme, "load_stylesheet"):
                        app.setStyleSheet("")
                        applied = True
                except Exception:
                    applied = False
            if not applied:
                app.setStyleSheet("")
        self.update_header_logo()
        self.apply_font_size(getattr(self, "current_font_size", "Small"))

    def update_header_logo(self):
        workflow_dir = self.edit_workflow_dir.text().strip() if hasattr(self, "edit_workflow_dir") else os.getcwd()
        resources = Path(workflow_dir) / "resources"
        dark_logo = resources / "darktheme_logo.png"
        light_logo = resources / "logo.png"
        chosen = dark_logo if (self.dark_mode_enabled and dark_logo.exists()) else light_logo
        if chosen.exists():
            pix = QPixmap(str(chosen))
            if not pix.isNull():
                self.logo_label.setPixmap(pix.scaled(220, 90, Qt.KeepAspectRatio, Qt.SmoothTransformation))
                self.logo_label.setText("")
                return
        self.logo_label.setPixmap(QPixmap())
        self.logo_label.setText("REPLAY")

    def update_header_height(self):
        self.header_widget.setFixedHeight(max(110, int(self.height() * 0.15)))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.update_header_height()

    def apply_font_size(self, size_label: str):
        scale_map = {"Small": 1.0, "Medium": 1.2, "Large": 1.4}
        self.current_font_size = size_label
        app = QApplication.instance()
        font = app.font()
        font.setPointSizeF(self.base_font_point_size * scale_map.get(size_label, 1.0))
        app.setFont(font)

        # Force the updated font through the existing widget tree as well.
        # This makes font changes apply immediately even when a dark theme
        # stylesheet is already active.
        for widget in app.allWidgets():
            widget.setFont(font)
            widget.updateGeometry()
            widget.update()

        if hasattr(self, "font_actions") and size_label in self.font_actions:
            self.font_actions[size_label].setChecked(True)
        self.table_files.setStyleSheet(
            "QTableWidget, QTableWidget::item, "
            "QTableWidget QLineEdit, QTableWidget QComboBox { color: black; }"
        )
        app.processEvents()
        self.update_header_height()

    def get_requested_memory_mb(self) -> int:
        value = int(self.spin_runtime_memory.value())
        unit = self.combo_runtime_memory_unit.currentText()
        return value if unit == "MB" else value * 1000

    def get_runtime_limits_mb(self) -> Tuple[int, int]:
        requested_mb = self.get_requested_memory_mb()
        memory_limit_mb = max(1, int(requested_mb * 0.99))
        memory_reservation_mb = max(1, int(requested_mb * 0.90))
        return memory_limit_mb, memory_reservation_mb

    def apply_max_runtime_resources(self):
        if self.total_system_memory_mb > 0:
            self.combo_runtime_memory_unit.setCurrentText("MB")
            self.spin_runtime_memory.setValue(self.total_system_memory_mb)
        if self.total_system_cores > 0:
            self.spin_runtime_cores.setValue(self.total_system_cores)
        self.update_resource_limits_status()

    def update_resource_limits_status(self) -> bool:
        errors = []
        mem_ok = True
        cores_ok = True
        requested_mem_mb = self.get_requested_memory_mb()
        requested_cores = int(self.spin_runtime_cores.value())

        if self.total_system_memory_mb > 0 and requested_mem_mb > self.total_system_memory_mb:
            mem_ok = False
            errors.append(f"Requested memory ({requested_mem_mb} MB) exceeds total system memory ({self.total_system_memory_mb} MB).")

        if self.total_system_cores > 0 and requested_cores > self.total_system_cores:
            cores_ok = False
            errors.append(f"Requested cores ({requested_cores}) exceeds total system cores ({self.total_system_cores}).")

        self.spin_runtime_memory.setStyleSheet("background-color: #ffdddd;" if not mem_ok else "")
        self.spin_runtime_cores.setStyleSheet("background-color: #ffdddd;" if not cores_ok else "")

        if errors:
            self.label_resource_error.setText(" ".join(errors))
            self.label_resource_error.setStyleSheet("color: red; font-weight: bold;")
            return False

        self.label_resource_error.setText(
            f"System resources detected: {self.total_system_memory_mb} MB RAM, {self.total_system_cores} cores."
            if self.total_system_memory_mb and self.total_system_cores else ""
        )
        self.label_resource_error.setStyleSheet("color: green;")
        return True

    # ------------------------------------------------------------------
    # Presets and defaults
    # ------------------------------------------------------------------
    def hg38_available(self) -> bool:
        wf = Path(self.edit_workflow_dir.text().strip())
        return (wf / "genome" / "hg38" / "hg38.fa").exists() or (wf / "genome" / "hg38" / "hg38.fa.gz").exists()

    def mm39_available(self) -> bool:
        wf = Path(self.edit_workflow_dir.text().strip())
        return (wf / "genome" / "mm39" / "mm39.fa").exists() or (wf / "genome" / "mm39" / "mm39.fa.gz").exists()

    def refresh_preset_availability(self):
        model = self.combo_preset.model()
        hg38_item = model.item(0)
        mm39_item = model.item(1)
        if hg38_item is not None:
            hg38_item.setEnabled(self.hg38_available())
        if mm39_item is not None:
            mm39_item.setEnabled(self.mm39_available())
        self.btn_download_hg38.setEnabled(not self.hg38_available())
        self.btn_download_mm39.setEnabled(not self.mm39_available())
        if self.combo_preset.currentText() == PRESET_HG38 and not self.hg38_available():
            self.combo_preset.setCurrentText(PRESET_CUSTOM)
        if self.combo_preset.currentText() == PRESET_MM39 and not self.mm39_available():
            self.combo_preset.setCurrentText(PRESET_CUSTOM)

    def apply_default_hg38(self):
        self.edit_genome.setText("genome/hg38/hg38.fa")
        self.edit_gap_file.setText("resources/hg38_gaplist.tsv")
        self.edit_target_file.setText("resources/TargetDataset_RT.txt")
        self.combo_normalization.setCurrentText("Quantile")
        self.combo_smoothing.setCurrentText("Loess")
        self.check_median_center.setChecked(True)
        self.check_mask_gaps.setChecked(True)
        self.select_window_sizes([5000, 20000])
        self.select_named_contigs([f"chr{i}" for i in range(1, 23)] + ["chrX"])

    def apply_default_mm39(self):
        self.edit_genome.setText("genome/mm39/mm39.fa")
        self.edit_gap_file.setText("")
        self.edit_target_file.setText("resources/TargetDataset_RT.txt")
        self.combo_normalization.setCurrentText("Quantile")
        self.combo_smoothing.setCurrentText("Loess")
        self.check_median_center.setChecked(True)
        self.check_mask_gaps.setChecked(False)
        self.select_window_sizes([5000, 20000])
        self.select_named_contigs([f"chr{i}" for i in range(1, 20)] + ["chrX"])

    def on_preset_changed(self, new_mode: str):
        self.current_profile_mode = new_mode
        if new_mode == PRESET_HG38:
            self.apply_default_hg38()
        elif new_mode == PRESET_MM39:
            self.apply_default_mm39()
        self.refresh_yaml_preview()
        self.update_file_check_status()

    def download_genome(self, genome_name: str):
        workflow_dir = Path(self.edit_workflow_dir.text().strip())
        if not workflow_dir:
            QMessageBox.warning(self, "Missing workflow directory", "Please set the workflow directory first.")
            return
        dest = workflow_dir / "genome" / genome_name
        ensure_parent(dest / "placeholder")
        url = f"http://hgdownload.soe.ucsc.edu/goldenPath/{genome_name}/bigZips/{genome_name}.fa.gz"
        cmd = ["wget", url, "-P", str(dest)]
        self.log("Downloading genome: " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=False)
            gz_path = dest / f"{genome_name}.fa.gz"
            fa_path = dest / f"{genome_name}.fa"
            if gz_path.exists():
                if shutil.which("gunzip"):
                    subprocess.run(["gunzip", "-f", str(gz_path)], check=False)
                else:
                    with gzip.open(gz_path, "rb") as src, open(fa_path, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                    gz_path.unlink(missing_ok=True)
                self.log(f"Unzipped genome to: {fa_path}")
            if fa_path.exists():
                self.edit_genome.setText(str(Path("genome") / genome_name / f"{genome_name}.fa"))
                self.generate_contigs_from_genome()
            self.refresh_preset_availability()
            if genome_name == "hg38" and self.hg38_available():
                self.combo_preset.setCurrentText(PRESET_HG38)
            elif genome_name == "mm39" and self.mm39_available():
                self.combo_preset.setCurrentText(PRESET_MM39)
            self.update_file_check_status()
            self.log(f"Finished genome download and setup for {genome_name}")
        except Exception as e:
            QMessageBox.critical(self, "Download failed", str(e))

    # ------------------------------------------------------------------
    # Generic UI helpers
    # ------------------------------------------------------------------
    def set_default_workflow_dir(self):
        cwd = Path(os.getcwd()).resolve()
        self.edit_workflow_dir.setText(str(cwd))
        if not self.edit_output_folder.text().strip():
            self.edit_output_folder.setText(str((cwd / "output" / "test").resolve()))
        if not self.edit_files_samples_folder.text().strip():
            self.edit_files_samples_folder.setText(str((cwd / "samples" / "run1").resolve()))
        self.update_header_logo()

    def on_workflow_dir_changed(self):
        self.refresh_preset_availability()
        self.populate_contigs_from_expected_file(auto_select_all=True)
        self.update_file_check_status()
        self.update_results_overview_labels()
        self.update_header_logo()

    def on_input_dir_changed(self):
        self.update_file_check_status()
        self.update_file_action_buttons()

    def on_genome_changed(self):
        self.populate_contigs_from_expected_file(auto_select_all=True)
        self.update_file_check_status()

    def on_tab_changed(self, idx: int):
        if self.results_index < 0:
            return
        show_right = idx != self.results_index
        self.right_panel.setVisible(show_right)
        if show_right:
            self.main_splitter.setSizes([1200, 620])
        else:
            self.main_splitter.setSizes([1820, 0])

    def log(self, msg: str):
        self.log_box.append(msg)

    def choose_workflow_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Select workflow directory")
        if path:
            self.edit_workflow_dir.setText(path)

    def choose_input_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Select raw reads folder")
        if path:
            self.edit_input_dir.setText(path)

    def choose_pipeline_reads_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Select pipeline reads folder")
        if path:
            self.edit_files_samples_folder.setText(path)

    def choose_output_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Select output folder")
        if path:
            self.edit_output_folder.setText(path)

    def choose_reference_genome(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select reference genome")
        if path:
            self.edit_genome.setText(self.to_relative_if_possible(path))
            self.update_file_check_status()

    def choose_target_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select target file")
        if path:
            self.edit_target_file.setText(self.to_relative_if_possible(path))
            self.update_file_check_status()

    def choose_gap_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select gap file")
        if path:
            self.edit_gap_file.setText(self.to_relative_if_possible(path))
            self.update_file_check_status()

    def to_relative_if_possible(self, path: str) -> str:
        workflow_dir = self.edit_workflow_dir.text().strip()
        try:
            return str(Path(path).resolve().relative_to(Path(workflow_dir).resolve()))
        except Exception:
            return path

    def resolve_workflow_path(self, text: str) -> Path:
        p = Path(text)
        if p.is_absolute():
            return p
        return Path(self.edit_workflow_dir.text().strip()) / p

    def resolved_samples_folder_path(self) -> Optional[Path]:
        text = self.edit_files_samples_folder.text().strip()
        if not text:
            return None
        return Path(text).resolve()

    def sync_samples_folder_from_files_tab(self):
        self.update_file_check_status()
        self.refresh_yaml_preview()
        self.update_file_action_buttons()

    def current_read_mode(self) -> str:
        return "paired" if self.check_paired.isChecked() else "single"

    def on_read_mode_changed(self):
        sender = self.sender()
        if sender == self.check_single and self.check_single.isChecked():
            self.check_paired.blockSignals(True)
            self.check_paired.setChecked(False)
            self.check_paired.blockSignals(False)
        elif sender == self.check_paired and self.check_paired.isChecked():
            self.check_single.blockSignals(True)
            self.check_single.setChecked(False)
            self.check_single.blockSignals(False)
        if not self.check_single.isChecked() and not self.check_paired.isChecked():
            self.check_single.blockSignals(True)
            self.check_single.setChecked(True)
            self.check_single.blockSignals(False)
        self.update_barcode_defaults(force_if_blank=False)
        self.validate_file_table()
        self.refresh_yaml_preview()

    def update_barcode_defaults(self, force_if_blank: bool = False):
        mode = self.current_read_mode()
        current_b1 = self.edit_barcode1.text().strip()
        current_b2 = self.edit_barcode2.text().strip()
        if mode == "paired":
            self.edit_barcode2.setEnabled(True)
            if force_if_blank or current_b1 in {"", SINGLE_FWD_DEFAULT, PAIR_FWD_DEFAULT}:
                self.edit_barcode1.setText(PAIR_FWD_DEFAULT)
            if force_if_blank or current_b2 in {"", PAIR_REV_DEFAULT}:
                self.edit_barcode2.setText(PAIR_REV_DEFAULT)
        else:
            self.edit_barcode2.setEnabled(False)
            if force_if_blank or current_b1 in {"", SINGLE_FWD_DEFAULT, PAIR_FWD_DEFAULT}:
                self.edit_barcode1.setText(SINGLE_FWD_DEFAULT)
            self.edit_barcode2.setText("")

    def select_window_sizes(self, sizes: List[int]):
        wanted = {str(s) for s in sizes}
        existing = {self.list_window_sizes.item(i).text() for i in range(self.list_window_sizes.count())}
        for size in wanted - existing:
            self.list_window_sizes.addItem(QListWidgetItem(size))
        for i in range(self.list_window_sizes.count()):
            self.list_window_sizes.item(i).setSelected(self.list_window_sizes.item(i).text() in wanted)

    def select_named_contigs(self, contigs: List[str]):
        wanted = set(contigs)
        for i in range(self.list_contigs.count()):
            self.list_contigs.item(i).setSelected(self.list_contigs.item(i).text() in wanted)
        self.refresh_yaml_preview()

    # ------------------------------------------------------------------
    # YAML / config save
    # ------------------------------------------------------------------
    def load_yaml_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load YAML", filter="YAML Files (*.yaml *.yml)")
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = yaml.safe_load(fh) or {}
            self.populate_form_from_config(data)
            self.refresh_yaml_preview()
            self.log(f"Loaded YAML: {path}")
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e))

    def save_yaml_file(self):
        self.refresh_yaml_preview()
        path, _ = QFileDialog.getSaveFileName(self, "Save YAML", filter="YAML Files (*.yaml *.yml)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(self.yaml_preview.toPlainText())
            self.log(f"Saved YAML: {path}")
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

    def populate_form_from_config(self, cfg: Dict):
        self.spin_latency.setValue(int(cfg.get("latency-wait", 120)))
        out_folder = str(cfg.get("output_folder", ""))
        self.edit_output_folder.setText(out_folder if out_folder else str((Path(os.getcwd()) / "output" / "test").resolve()))
        samples_folder = str(cfg.get("samples_folder", ""))
        self.edit_files_samples_folder.setText(samples_folder if samples_folder else str((Path(os.getcwd()) / "samples" / "run1").resolve()))
        genome_path = str(cfg.get("genome", ""))
        self.edit_genome.setText(genome_path if not Path(genome_path).is_absolute() else self.to_relative_if_possible(genome_path))
        self.edit_barcode1.setText(str(cfg.get("barcode1", SINGLE_FWD_DEFAULT)))
        self.edit_barcode2.setText(str(cfg.get("barcode2", "")))

        selected_ws = {str(v) for v in cfg.get("window_sizes", [])}
        for i in range(self.list_window_sizes.count()):
            self.list_window_sizes.item(i).setSelected(self.list_window_sizes.item(i).text() in selected_ws)

        reads = set(cfg.get("reads", []))
        self.check_single.blockSignals(True)
        self.check_paired.blockSignals(True)
        self.check_single.setChecked(reads == {"R1"})
        self.check_paired.setChecked(reads == {"R1", "R2"})
        if not self.check_single.isChecked() and not self.check_paired.isChecked():
            self.check_single.setChecked(True)
        self.check_single.blockSignals(False)
        self.check_paired.blockSignals(False)

        rt = cfg.get("rt_processing", {})
        self.combo_normalization.setCurrentText(NORMALIZATION_VALUE_TO_DISPLAY.get(str(rt.get("normalization", "quantile")).lower(), "Quantile"))
        self.combo_smoothing.setCurrentText(SMOOTHING_VALUE_TO_DISPLAY.get(str(rt.get("smoothing", "loess")).lower(), "Loess"))
        self.check_median_center.setChecked(bool(rt.get("median_center", True)))
        self.edit_target_file.setText(str(rt.get("target_file", "resources/TargetDataset_RT.txt")))
        self.edit_gap_file.setText(str(rt.get("gap_file", "resources/hg38_gaplist.tsv")))
        self.check_mask_gaps.setChecked(bool(rt.get("mask_gaps", True)))
        self.spin_loess_span.setValue(int(rt.get("loess_span_bp", 500000)))
        self.spin_clip_min.setValue(int(rt.get("clip_min", -8)))
        self.spin_clip_max.setValue(int(rt.get("clip_max", 8)))

        existing_contigs = cfg.get("contigs", []) or []
        self.populate_contigs_from_expected_file(auto_select_all=(len(existing_contigs) == 0))
        if existing_contigs:
            self.select_named_contigs(existing_contigs)

        self.update_barcode_defaults(force_if_blank=False)
        self.update_file_check_status()
        self.update_file_action_buttons()

    def inferred_target_files(self) -> List[str]:
        valid_targets = []
        for row in range(self.table_files.rowCount()):
            include_widget = self.table_files.cellWidget(row, 1)
            include = include_widget.currentText() if include_widget else "Yes"
            proposed_item = self.table_files.item(row, 6)
            row_check_item = self.table_files.item(row, 8)
            group_check_item = self.table_files.item(row, 9)
            proposed = proposed_item.text() if proposed_item else ""
            row_check = row_check_item.text() if row_check_item else ""
            group_check = group_check_item.text() if group_check_item else ""
            if include == "Yes" and proposed and row_check == "OK" and group_check == "OK":
                valid_targets.append(proposed)
        return sorted(valid_targets)

    def get_selected_contigs(self) -> List[str]:
        return [item.text() for item in self.list_contigs.selectedItems()]

    def form_to_config(self) -> Dict:
        reads = ["R1"] if self.current_read_mode() == "single" else ["R1", "R2"]
        window_sizes = [int(self.list_window_sizes.item(i).text()) for i in range(self.list_window_sizes.count()) if self.list_window_sizes.item(i).isSelected()]
        if not window_sizes:
            raise ValueError("At least one window size must be selected.")
        genome_full = str(self.resolve_workflow_path(self.edit_genome.text().strip()).resolve())
        output_full = str(Path(self.edit_output_folder.text().strip()).resolve()) if self.edit_output_folder.text().strip() else ""
        samples_full = str(self.resolved_samples_folder_path()) if self.resolved_samples_folder_path() is not None else ""
        return {
            "latency-wait": self.spin_latency.value(),
            "genome": genome_full,
            "output_folder": output_full,
            "samples_folder": samples_full,
            "window_sizes": window_sizes,
            "reads": reads,
            "barcode1": self.edit_barcode1.text().strip(),
            "barcode2": self.edit_barcode2.text().strip() if self.current_read_mode() == "paired" else "",
            "contigs": self.get_selected_contigs(),
            "rt_processing": {
                "normalization": NORMALIZATION_DISPLAY_TO_VALUE[self.combo_normalization.currentText()],
                "smoothing": SMOOTHING_DISPLAY_TO_VALUE[self.combo_smoothing.currentText()],
                "median_center": self.check_median_center.isChecked(),
                "target_file": self.edit_target_file.text().strip(),
                "gap_file": self.edit_gap_file.text().strip(),
                "mask_gaps": self.check_mask_gaps.isChecked(),
                "loess_span_bp": self.spin_loess_span.value(),
                "clip_min": self.spin_clip_min.value(),
                "clip_max": self.spin_clip_max.value(),
            },
            "target": self.inferred_target_files(),
        }

    def refresh_yaml_preview(self):
        try:
            cfg = self.form_to_config()
            self.current_config = cfg
            self.yaml_preview.setPlainText(dump_yaml_text(cfg))
            self.populate_results_bin_sizes()
            self.update_results_overview_labels()
        except Exception as e:
            self.log(f"YAML preview refresh failed: {e}")

    def save_generated_config(self):
        workflow_dir = self.edit_workflow_dir.text().strip()
        if not workflow_dir:
            raise ValueError("Missing workflow directory")
        self.refresh_yaml_preview()
        config_data = load_yaml_text(self.yaml_preview.toPlainText())
        cfg_path = Path(workflow_dir) / "config" / "generated_config.yaml"
        ensure_parent(cfg_path)
        with open(cfg_path, "w", encoding="utf-8") as fh:
            fh.write(dump_yaml_text(config_data))
        self.log(f"Saved generated config: {cfg_path}")

    # ------------------------------------------------------------------
    # Window sizes / contigs
    # ------------------------------------------------------------------
    def add_window_size(self):
        text = self.edit_new_window_size.text().strip()
        if not text.isdigit():
            QMessageBox.warning(self, "Invalid size", "Window size must be a positive integer.")
            return
        existing = {self.list_window_sizes.item(i).text() for i in range(self.list_window_sizes.count())}
        if text not in existing:
            item = QListWidgetItem(text)
            item.setSelected(True)
            self.list_window_sizes.addItem(item)
        self.edit_new_window_size.clear()
        self.refresh_yaml_preview()

    def remove_selected_window_sizes(self):
        for item in self.list_window_sizes.selectedItems():
            self.list_window_sizes.takeItem(self.list_window_sizes.row(item))
        self.refresh_yaml_preview()

    def expected_contigs_file_path(self) -> Optional[Path]:
        workflow_dir = self.edit_workflow_dir.text().strip()
        genome_text = self.edit_genome.text().strip()
        if not workflow_dir or not genome_text:
            return None
        return Path(workflow_dir) / "resources" / "contigs" / f"{genome_id_from_path(genome_text)}_contigs.txt"

    def generate_contigs_file_for_genome(self, genome_path: Path):
        out_path = self.expected_contigs_file_path()
        if out_path is None:
            return
        ensure_parent(out_path)
        contigs = []
        with open(genome_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    if contig:
                        contigs.append(contig)
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("\n".join(contigs) + ("\n" if contigs else ""))
        self.log(f"Wrote contigs file: {out_path}")

    def generate_contigs_from_genome(self):
        genome_path = self.resolve_workflow_path(self.edit_genome.text().strip())
        if not genome_path.exists():
            QMessageBox.warning(self, "Genome missing", f"Reference genome not found: {genome_path}")
            return
        try:
            self.generate_contigs_file_for_genome(genome_path)
            self.populate_contigs_from_expected_file(auto_select_all=True)
            self.update_file_check_status()
        except Exception as e:
            QMessageBox.critical(self, "Contig generation failed", str(e))

    def populate_contigs_from_expected_file(self, auto_select_all: bool = True):
        contigs_path = self.expected_contigs_file_path()
        self.list_contigs.clear()
        if contigs_path is None or not contigs_path.exists():
            return
        try:
            with open(contigs_path, "r", encoding="utf-8", errors="replace") as fh:
                contigs = [line.strip() for line in fh if line.strip()]
            for contig in contigs:
                item = QListWidgetItem(contig)
                self.list_contigs.addItem(item)
                if auto_select_all:
                    item.setSelected(True)
            self.refresh_yaml_preview()
        except Exception as e:
            self.log(f"Failed to load contigs file: {e}")

    def select_all_contigs(self):
        for i in range(self.list_contigs.count()):
            self.list_contigs.item(i).setSelected(True)
        self.refresh_yaml_preview()

    def deselect_all_contigs(self):
        for i in range(self.list_contigs.count()):
            self.list_contigs.item(i).setSelected(False)
        self.refresh_yaml_preview()

    def write_selected_contigs_to_file(self):
        out_path = self.expected_contigs_file_path()
        if out_path is None:
            QMessageBox.warning(self, "Missing contigs path", "Set workflow directory and genome first.")
            return
        ensure_parent(out_path)
        selected = self.get_selected_contigs()
        try:
            with open(out_path, "w", encoding="utf-8") as fh:
                fh.write("\n".join(selected) + ("\n" if selected else ""))
            self.log(f"Overwrote contigs file with {len(selected)} selected contigs: {out_path}")
            self.populate_contigs_from_expected_file(auto_select_all=False)
            self.select_named_contigs(selected)
            self.update_file_check_status()
        except Exception as e:
            QMessageBox.critical(self, "Contig write failed", str(e))

    # ------------------------------------------------------------------
    # File table / validation
    # ------------------------------------------------------------------
    def _make_line_edit(self, text: str) -> QLineEdit:
        widget = QLineEdit(text)
        widget.editingFinished.connect(self.validate_file_table)
        return widget

    def _make_combo(self, values: List[str], current: str) -> QComboBox:
        widget = QComboBox()
        widget.addItems(values)
        idx = widget.findText(current)
        if idx >= 0:
            widget.setCurrentIndex(idx)
        widget.currentIndexChanged.connect(self.validate_file_table)
        return widget

    def _make_include_combo(self, current: str = "Yes") -> QComboBox:
        widget = QComboBox()
        widget.addItems(["Yes", "No"])
        widget.setCurrentText(current)
        widget.currentIndexChanged.connect(self.validate_file_table)
        return widget

    def set_table_item(self, row: int, col: int, text: str):
        item = self.table_files.item(row, col)
        if item is None:
            item = QTableWidgetItem(text)
            self.table_files.setItem(row, col, item)
        else:
            item.setText(text)

    def get_row_values(self, row: int):
        original = self.table_files.item(row, 0).text()
        include_widget = self.table_files.cellWidget(row, 1)
        include = include_widget.currentText().strip() if include_widget else "Yes"
        desc_widget = self.table_files.cellWidget(row, 2)
        rep_widget = self.table_files.cellWidget(row, 3)
        el_widget = self.table_files.cellWidget(row, 4)
        read_widget = self.table_files.cellWidget(row, 5)
        desc = desc_widget.text().strip() if desc_widget else ""
        replicate = rep_widget.text().strip() if rep_widget else ""
        el = el_widget.currentText().strip() if el_widget else ""
        read = read_widget.currentText().strip() if read_widget else ""
        return include, original, desc, replicate, el, read

    def color_row(self, row: int, color: Optional[QColor]):
        for col in [0, 6, 7, 8, 9]:
            item = self.table_files.item(row, col)
            if item is not None:
                if color is None:
                    item.setBackground(QColor())
                else:
                    item.setBackground(color)
        for col in [1, 2, 3, 4, 5]:
            widget = self.table_files.cellWidget(row, col)
            if widget is not None:
                widget.setStyleSheet("") if color is None else widget.setStyleSheet(f"background-color: {color.name()};")

    def scan_input_folder(self):
        input_text = self.edit_input_dir.text().strip()
        input_dir = Path(input_text) if input_text else Path("")
        if not input_text:
            QMessageBox.warning(self, "Missing raw reads folder", "Please select a raw reads folder.")
            return
        if not input_dir.exists() or not input_dir.is_dir():
            QMessageBox.warning(self, "Invalid raw reads folder", f"Not a valid directory: {input_dir}")
            return
        fastqs = sorted([x.name for x in input_dir.iterdir() if x.is_file() and x.name.lower().endswith(".fastq.gz")])
        self.table_files.setRowCount(0)
        self.rename_plan = []
        for fn in fastqs:
            row = self.table_files.rowCount()
            self.table_files.insertRow(row)
            self.set_table_item(row, 0, fn)
            self.table_files.setCellWidget(row, 1, self._make_include_combo("Yes"))
            valid, _, parsed = canonical_validate_filename(fn)
            if valid and parsed:
                desc, replicate, el, read = parsed["desc"], parsed["sample"], parsed["el"], parsed["read"]
                reason = "Canonical filename"
            else:
                plan = attempted_parse_and_propose(fn)
                reason = plan.get("reason", "")
                proposed = plan.get("proposed", "")
                match = CANONICAL_RE.match(proposed) if proposed else None
                if match:
                    desc, replicate, el, read = match.group("desc"), match.group("sample"), match.group("el"), match.group("read")
                else:
                    desc, replicate, el, read = strip_fastq_suffix(fn), "", "E", "R1"
            self.table_files.setCellWidget(row, 2, self._make_line_edit(desc))
            self.table_files.setCellWidget(row, 3, self._make_line_edit(replicate))
            self.table_files.setCellWidget(row, 4, self._make_combo(["E", "L"], el))
            self.table_files.setCellWidget(row, 5, self._make_combo(["R1", "R2"], read))
            self.set_table_item(row, 6, "")
            self.set_table_item(row, 7, "")
            self.set_table_item(row, 8, reason)
            self.set_table_item(row, 9, "")
        self.log(f"Scanned, auto-populated, and validating {len(fastqs)} FASTQ files from {input_dir}")
        self.validate_file_table()
        self.update_file_check_status()

    def build_proposed_name(self, desc: str, replicate: str, el: str, read: str) -> str:
        return f"{desc}_{replicate}_{el}_{read}.fastq.gz"

    def validate_file_table(self):
        selected_reads = ["R1"] if self.current_read_mode() == "single" else ["R1", "R2"]
        combo_counts = {}
        group_presence = {}
        row_cache = []
        for row in range(self.table_files.rowCount()):
            include, original, desc, replicate, el, read = self.get_row_values(row)
            row_messages = []
            proposed = ""
            if include == "Yes":
                if not desc:
                    row_messages.append("Missing description")
                if "_" in desc:
                    row_messages.append("Description cannot contain underscores")
                if not replicate.isdigit():
                    row_messages.append("Replicate must be numeric")
                if el not in {"E", "L"}:
                    row_messages.append("E/L must be E or L")
                if read not in {"R1", "R2"}:
                    row_messages.append("Read must be R1 or R2")
                if read not in selected_reads:
                    row_messages.append("Read not enabled by selected pipeline reads")
                if not row_messages:
                    proposed = self.build_proposed_name(desc, replicate, el, read)
                    combo_key = (desc, replicate, el, read)
                    combo_counts[combo_key] = combo_counts.get(combo_key, 0) + 1
                    group_presence.setdefault((desc, replicate), set()).add((el, read))
            else:
                if desc and replicate and el and read:
                    proposed = self.build_proposed_name(desc, replicate, el, read)
            row_cache.append({
                "row": row,
                "include": include,
                "original": original,
                "desc": desc,
                "replicate": replicate,
                "el": el,
                "read": read,
                "proposed": proposed,
                "row_messages": row_messages,
            })

        self.rename_plan = []
        pipeline_dir = self.resolved_samples_folder_path()
        for info in row_cache:
            row_messages = list(info["row_messages"])
            group_messages = []
            proposed = info["proposed"]
            copy_success = False
            if proposed and pipeline_dir is not None:
                copy_success = (pipeline_dir / proposed).exists()

            if info["include"] == "Yes" and proposed:
                combo_key = (info["desc"], info["replicate"], info["el"], info["read"])
                if combo_counts.get(combo_key, 0) > 1:
                    group_messages.append("Duplicate desc/replicate/E-L/read combination")
                required = {(phase, r) for phase in ["E", "L"] for r in selected_reads}
                present = group_presence.get((info["desc"], info["replicate"]), set())
                missing = sorted(required - present)
                if missing:
                    group_messages.append("Missing pair(s): " + ", ".join([f"{p}_{r}" for p, r in missing]))

            if info["include"] == "No":
                row_check = "Excluded"
                group_check = "Excluded"
                row_ok = True
                group_ok = True
            else:
                row_check = "OK" if not row_messages else "; ".join(row_messages)
                group_check = "OK" if proposed and not group_messages and not row_messages else "; ".join(group_messages)
                row_ok = not row_messages
                group_ok = proposed and not group_messages and row_ok

            self.set_table_item(info["row"], 6, proposed)
            self.set_table_item(info["row"], 7, "Yes" if copy_success else "")
            self.set_table_item(info["row"], 8, row_check)
            self.set_table_item(info["row"], 9, group_check)

            status = "valid" if (row_ok and group_ok and info["original"] == proposed) else ("proposed" if (row_ok and proposed) else "failed")
            self.rename_plan.append({
                "status": status,
                "include": info["include"],
                "original": info["original"],
                "proposed": proposed,
                "reason": row_check if row_check != "OK" else group_check,
                "row_ok": row_ok,
                "group_ok": group_ok,
            })

            if info["include"] == "No":
                self.color_row(info["row"], ROW_EXCLUDED_COLOR)
            elif copy_success:
                self.color_row(info["row"], ROW_OK_COLOR)
            elif row_ok and group_ok and info["original"] == proposed:
                self.color_row(info["row"], ROW_OK_COLOR)
            elif row_ok and proposed:
                self.color_row(info["row"], ROW_WARN_COLOR)
            else:
                self.color_row(info["row"], ROW_ERROR_COLOR)

        self.log(f"Files table validated; inferred {len(self.inferred_target_files())} valid target files")
        self.refresh_yaml_preview()
        self.update_file_check_status()
        self.update_file_action_buttons()

    def destinations_would_overwrite(self) -> bool:
        samples_base = self.resolved_samples_folder_path()
        return bool(samples_base and any((samples_base / target).exists() for target in self.inferred_target_files()))

    def table_fully_validated(self) -> bool:
        included_plans = [plan for plan in self.rename_plan if plan.get("include") == "Yes"]
        targets = self.inferred_target_files()
        return len(targets) > 0 and all(plan.get("row_ok") and plan.get("group_ok") for plan in included_plans)

    def input_matches_samples_folder(self) -> bool:
        input_text = self.edit_input_dir.text().strip()
        samples_base = self.resolved_samples_folder_path()
        if not input_text or samples_base is None:
            return False
        try:
            return Path(input_text).resolve() == samples_base.resolve()
        except Exception:
            return False

    def update_file_action_buttons(self):
        validated = self.table_fully_validated()
        overwrite_needed = self.destinations_would_overwrite()
        overwrite_allowed = self.check_overwrite.isChecked()
        same_dir = self.input_matches_samples_folder()
        self.btn_copy_or_rename.setText("Rename" if same_dir else "Copy")
        can_copy = validated and (overwrite_allowed or not overwrite_needed)
        self.btn_copy_or_rename.setEnabled(can_copy)
        self.btn_symlink.setEnabled(can_copy and not same_dir)
        self.label_overwrite_warning.setText("would overwrite" if overwrite_needed and not overwrite_allowed else "")
        if not validated:
            self.label_action_status.setText("Validation errors present")
            self.label_action_status.setStyleSheet("color: red;")
        else:
            self.label_action_status.setText("Ready")
            self.label_action_status.setStyleSheet("color: green;")

    def apply_primary_file_action(self):
        self.apply_rename_plan(mode="rename" if self.input_matches_samples_folder() else "copy")

    def apply_rename_plan(self, mode: str = "copy"):
        input_text = self.edit_input_dir.text().strip()
        if not input_text:
            QMessageBox.warning(self, "Missing raw reads folder", "Please select a raw reads folder.")
            return
        self.validate_file_table()
        if not self.rename_plan:
            QMessageBox.information(self, "No plan", "Please scan files first.")
            return
        overwrite = self.check_overwrite.isChecked()
        input_base = Path(input_text).resolve()
        samples_base = self.resolved_samples_folder_path()
        if samples_base is None:
            QMessageBox.warning(self, "Missing pipeline reads folder", "Please set a pipeline reads folder.")
            return
        samples_base.mkdir(parents=True, exist_ok=True)
        changed = skipped = errors = 0
        for plan in self.rename_plan:
            if plan.get("include") != "Yes":
                skipped += 1
                continue
            original = plan["original"]
            proposed = plan.get("proposed", "")
            if not proposed or not plan.get("row_ok") or not plan.get("group_ok"):
                skipped += 1
                continue
            src = input_base / original
            dst = samples_base / proposed
            self.log(f"{mode.upper()}: {original} -> {dst}")
            if not src.exists():
                self.log(f"Missing source file: {src}")
                errors += 1
                continue
            if dst.exists() and not overwrite:
                self.log(f"Destination exists, skipping: {dst}")
                errors += 1
                continue
            try:
                if mode == "copy":
                    if dst.exists() and overwrite and (dst.is_symlink() or dst.is_file()):
                        dst.unlink()
                    shutil.copy2(src, dst)
                elif mode == "symlink":
                    if dst.exists() and overwrite:
                        dst.unlink()
                    dst.symlink_to(src)
                elif mode == "rename":
                    if dst.exists() and overwrite:
                        dst.unlink()
                    src.rename(dst)
                changed += 1
            except Exception as e:
                self.log(f"Failed {mode} for {original}: {e}")
                errors += 1
        self.log(f"File action finished: mode={mode}, changed={changed}, skipped={skipped}, errors={errors}")
        self.validate_file_table()
        self.update_file_check_status()
        self.update_file_action_buttons()

    # ------------------------------------------------------------------
    # File checks
    # ------------------------------------------------------------------
    def collect_required_file_checks(self) -> List[Tuple[str, Path, bool]]:
        checks = []
        workflow_dir = Path(self.edit_workflow_dir.text().strip())
        samples_base = self.resolved_samples_folder_path()
        genome = self.resolve_workflow_path(self.edit_genome.text().strip()) if self.edit_genome.text().strip() else Path("")
        target_file = self.resolve_workflow_path(self.edit_target_file.text().strip()) if self.edit_target_file.text().strip() else Path("")
        gap_file = self.resolve_workflow_path(self.edit_gap_file.text().strip()) if self.edit_gap_file.text().strip() else Path("")
        contigs_file = self.expected_contigs_file_path()
        if samples_base is not None:
            for target in self.inferred_target_files():
                checks.append((f"sample target: {target}", samples_base / target, True))
        checks.append(("reference genome", genome, True))
        if NORMALIZATION_DISPLAY_TO_VALUE[self.combo_normalization.currentText()] == "quantile":
            checks.append(("target_file", target_file, True))
        if self.check_mask_gaps.isChecked():
            checks.append(("gap_file", gap_file, True))
        checks.append(("scripts/rpkm.sh", workflow_dir / "scripts" / "rpkm.sh", True))
        checks.append(("scripts/normalize_smooth_rt.py", workflow_dir / "scripts" / "normalize_smooth_rt.py", True))
        if contigs_file is not None:
            checks.append(("contigs file", contigs_file, True))
        checks.append(("Snakefile", workflow_dir / "Snakefile", True))
        checks.append(("snakemake_repliseq.sif", workflow_dir / "snakemake_repliseq.sif", True))
        return checks

    def update_file_check_status(self):
        self.update_resource_limits_status()
        missing = []
        for label, path, required in self.collect_required_file_checks():
            if required and (not str(path) or not path.exists()):
                missing.append(f"{label}: {path}")
        if missing:
            self.label_file_check_status.setText("Missing required files")
            self.label_file_check_status.setStyleSheet("color: red;")
            self.log("File check found missing items:")
            for item in missing:
                self.log(f"  {item}")
            return False
        self.label_file_check_status.setText("All required files found")
        self.label_file_check_status.setStyleSheet("color: green;")
        return True

    # ------------------------------------------------------------------
    # Snakemake execution
    # ------------------------------------------------------------------
    def reset_progress(self):
        self.progress_bar.setValue(0)
        self.progress_bar.setStyleSheet("")
        self.progress_label.setStyleSheet("")
        self.progress_label.setText("No Snakemake progress detected yet")
        self.run_reached_100 = False

    def mark_progress_error(self):
        self.progress_bar.setStyleSheet("QProgressBar::chunk { background-color: red; }")
        self.progress_label.setStyleSheet("color: red; font-weight: bold;")
        text = self.progress_label.text()
        if ":ERROR" not in text:
            if text:
                self.progress_label.setText(text + ":ERROR")
            else:
                self.progress_label.setText("ERROR")

    def process_runtime_line(self, line: str):
        self.log(line)
        self.update_progress_from_line(line)
        if NOTHING_TO_DO_TEXT in line:
            self.populate_results_tab()
        if "Local run finished with exit_code=1" in line:
            self.mark_progress_error()

    def update_progress_from_line(self, line: str):
        match = PROGRESS_RE.search(line)
        if not match:
            return
        done = int(match.group("done"))
        total = int(match.group("total"))
        pct = int(match.group("pct"))
        self.progress_bar.setValue(max(0, min(100, pct)))
        self.progress_label.setText(f"{done} of {total} steps ({pct}%) done")
        if done == total and pct == 100:
            self.run_reached_100 = True

    def start_log_tailing(self):
        self.current_log_path = None
        self.current_log_pos = 0
        self.log_tail_timer.start()

    def stop_log_tailing(self):
        self.log_tail_timer.stop()
        self.current_log_path = None
        self.current_log_pos = 0

    def find_latest_snakemake_log(self, workflow_dir: str) -> Optional[Path]:
        log_dir = Path(workflow_dir) / ".snakemake" / "log"
        if not log_dir.exists() or not log_dir.is_dir():
            return None
        eligible = [p for p in log_dir.iterdir() if p.is_file() and p.stat().st_mtime >= self.local_run_started_at]
        if not eligible:
            return None
        return max(eligible, key=lambda p: p.stat().st_mtime)

    def poll_latest_snakemake_log(self):
        workflow_dir = self.edit_workflow_dir.text().strip()
        if not workflow_dir or self.local_run_started_at == 0.0 or self.pending_unlock:
            return
        latest = self.find_latest_snakemake_log(workflow_dir)
        if latest is None:
            return
        if self.current_log_path != latest:
            self.current_log_path = latest
            self.current_log_pos = 0
            self.log(f"Switched to latest Snakemake log: {latest.name}")
        try:
            with open(self.current_log_path, "r", encoding="utf-8", errors="replace") as fh:
                fh.seek(self.current_log_pos)
                new_text = fh.read()
                self.current_log_pos = fh.tell()
        except Exception as e:
            self.log(f"Failed to read Snakemake log: {e}")
            return
        if not new_text:
            return
        for line in new_text.splitlines():
            if line.strip():
                self.process_runtime_line(line)

    def on_run_process_stdout(self):
        if self.run_process is None:
            return
        text = bytes(self.run_process.readAllStandardOutput()).decode("utf-8", errors="replace")
        for line in text.splitlines():
            if line.strip():
                self.process_runtime_line(line)

    def on_run_process_stderr(self):
        if self.run_process is None:
            return
        text = bytes(self.run_process.readAllStandardError()).decode("utf-8", errors="replace")
        for line in text.splitlines():
            if line.strip():
                self.process_runtime_line(line)

    def on_run_finished(self, exit_code: int, exit_status):
        self.poll_latest_snakemake_log()
        self.stop_log_tailing()
        if self.pending_unlock:
            self.log(f"Snakemake unlock finished with exit_code={exit_code}")
        else:
            self.log(f"Local run finished with exit_code={exit_code}")
            if exit_code == 0:
                self.progress_bar.setValue(100)
                if "No Snakemake progress" in self.progress_label.text():
                    self.progress_label.setText("Run completed")
                if self.run_reached_100:
                    self.populate_results_tab()
            else:
                self.mark_progress_error()
        self.run_process = None
        self.pending_unlock = False
        self.local_run_started_at = 0.0
        self.update_file_check_status()

    def build_snakemake_args(self, unlock: bool = False) -> List[str]:
        workflow_dir = self.edit_workflow_dir.text().strip()
        workdir = os.path.join(workflow_dir, "apptainer")
        os.makedirs(workdir, exist_ok=True)
    
        _, memory_reservation_mb = self.get_runtime_limits_mb()
        requested_cores = int(self.spin_runtime_cores.value())
    
        args = [
            "exec",
            "--bind",
            workflow_dir,
            "--pwd",
            workflow_dir,
            "-c",
            "--workdir",
            workdir + "/",
            "--unsquash",
            "snakemake_repliseq.sif",
            "snakemake",
            "--snakefile",
            "Snakefile",
            "--cores",
            str(requested_cores),
            "--resources",
            f"mem_mb={memory_reservation_mb}",
            "--configfile",
            "config/generated_config.yaml",
            "--latency-wait",
            str(self.spin_latency.value()),
        ]
    
        if unlock:
            args.append("--unlock")
    
        return args



    def run_local(self):
        workflow_dir = self.edit_workflow_dir.text().strip()
        if not workflow_dir:
            QMessageBox.warning(self, "Missing workflow directory", "Please set the workflow directory.")
            return
        singularity_cmd = shutil.which("singularity") or shutil.which("apptainer")
        if singularity_cmd is None:
            QMessageBox.warning(self, "Container runtime not found", "Neither singularity nor apptainer is available in PATH.")
            return
        if not self.update_file_check_status():
            QMessageBox.warning(self, "Missing required files", "Resolve file check errors before running locally.")
            return
        if not self.update_resource_limits_status():
            QMessageBox.warning(self, "Invalid runtime resources", "Requested memory or cores exceed system resources.")
            return
        if self.run_process is not None and self.run_process.state() != QProcess.NotRunning:
            QMessageBox.information(self, "Run already active", "A process is already in progress.")
            return
        try:
            self.save_generated_config()
            self.log("Running local Snakemake:")
            self.reset_progress()
            self.local_run_started_at = time.time()
            self.pending_unlock = False
            self.run_process = QProcess(self)
            self.run_process.setWorkingDirectory(workflow_dir)
            self.run_process.readyReadStandardOutput.connect(self.on_run_process_stdout)
            self.run_process.readyReadStandardError.connect(self.on_run_process_stderr)
            self.run_process.finished.connect(self.on_run_finished)
            program = singularity_cmd
            args = self.build_snakemake_args(unlock=False)
            self.log(" ".join([program] + args))
            self.start_log_tailing()
            self.run_process.start(program, args)
            if not self.run_process.waitForStarted(3000):
                self.stop_log_tailing()
                raise RuntimeError("Failed to start local Snakemake process.")
        except Exception as e:
            self.stop_log_tailing()
            self.local_run_started_at = 0.0
            QMessageBox.critical(self, "Run failed", str(e))

    def unlock_snakemake(self):
        workflow_dir = self.edit_workflow_dir.text().strip()
        if not workflow_dir:
            QMessageBox.warning(self, "Missing workflow directory", "Please set the workflow directory.")
            return
        singularity_cmd = shutil.which("singularity") or shutil.which("apptainer")
        if singularity_cmd is None:
            QMessageBox.warning(self, "Container runtime not found", "Neither singularity nor apptainer is available in PATH.")
            return
        if self.run_process is not None and self.run_process.state() != QProcess.NotRunning:
            QMessageBox.information(self, "Process already active", "Another process is already running.")
            return
        if not self.update_resource_limits_status():
            QMessageBox.warning(self, "Invalid runtime resources", "Requested memory or cores exceed system resources.")
            return
        try:
            self.save_generated_config()
            self.log("Running Snakemake unlock:")
            self.pending_unlock = True
            self.local_run_started_at = time.time()
            self.run_process = QProcess(self)
            self.run_process.setWorkingDirectory(workflow_dir)
            self.run_process.readyReadStandardOutput.connect(self.on_run_process_stdout)
            self.run_process.readyReadStandardError.connect(self.on_run_process_stderr)
            self.run_process.finished.connect(self.on_run_finished)
            program = singularity_cmd
            args = self.build_snakemake_args(unlock=True)
            self.log(" ".join([program] + args))
            self.run_process.start(program, args)
            if not self.run_process.waitForStarted(3000):
                raise RuntimeError("Failed to start unlock process.")
        except Exception as e:
            self.pending_unlock = False
            self.local_run_started_at = 0.0
            QMessageBox.critical(self, "Unlock failed", str(e))

    # ------------------------------------------------------------------
    # Results tab
    # ------------------------------------------------------------------
    def populate_results_bin_sizes(self):
        current = self.combo_results_bin_size.currentText()
        sizes = [str(size) for size in sorted([int(self.list_window_sizes.item(i).text()) for i in range(self.list_window_sizes.count()) if self.list_window_sizes.item(i).isSelected()])]
        self.combo_results_bin_size.blockSignals(True)
        self.combo_results_bin_size.clear()
        self.combo_results_bin_size.addItems(sizes)
        if current in sizes:
            self.combo_results_bin_size.setCurrentText(current)
        self.combo_results_bin_size.blockSignals(False)

    def collect_result_samples(self) -> List[Tuple[str, str]]:
        pairs = set()
        for target in self.inferred_target_files():
            m = CANONICAL_RE.match(target)
            if m:
                pairs.add((m.group("desc"), m.group("sample")))
        return sorted(pairs)

    def result_base_dir(self) -> Path:
        text = self.edit_output_folder.text().strip()
        return Path(text).resolve() if text else Path(os.getcwd()).resolve()

    def update_results_overview_labels(self):
        outdir = self.result_base_dir()
        self.label_norm_results_dir.setText(f"Normalized Results: {dir_with_trailing_slash(outdir / 'bg' / 'norm')}")
        self.label_raw_results_dir.setText(f"Raw Results: {dir_with_trailing_slash(outdir / 'bg' / 'log2')}")

    def populate_results_tab(self):
        samples = self.collect_result_samples()
        self.results_sample_list.clear()
        for desc, rep in samples:
            self.results_sample_list.addItem(f"{desc}_{rep}")
        self.populate_results_bin_sizes()
        self.update_results_overview_labels()
        self.tabs.setTabEnabled(self.results_index, True)
        if samples:
            self.results_sample_list.setCurrentRow(0)
            self.update_results_view()
        self.log(f"Results tab populated with {len(samples)} sample pair(s)")

    def set_image_label(self, label: QLabel, path: Path):
        if path.exists():
            pixmap = QPixmap(str(path))
            if not pixmap.isNull():
                label.setPixmap(pixmap.scaled(280, 220, Qt.KeepAspectRatio, Qt.SmoothTransformation))
                label.setText("")
                return
        label.setPixmap(QPixmap())
        label.setText(f"Missing\n{path.name}")

    def update_results_visibility(self):
        self.col_reads["box"].setVisible(self.toggle_reads.isChecked())
        self.col_windows["box"].setVisible(self.toggle_windows.isChecked())
        self.col_acf["box"].setVisible(self.toggle_acf.isChecked())
        self.col_raw["box"].setVisible(self.toggle_raw.isChecked())
        self.col_rawkde["box"].setVisible(self.toggle_rawkde.isChecked())
        self.col_kde["box"].setVisible(self.toggle_kde.isChecked())

    def update_results_view(self):
        row = self.results_sample_list.currentRow()
        if row < 0:
            self.results_header.setText("No sample selected")
            return
        item_text = self.results_sample_list.item(row).text()
        desc, rep = item_text.rsplit("_", 1)
        bin_size = self.combo_results_bin_size.currentText() or ""
        genome_id = genome_id_from_path(self.edit_genome.text().strip())
        outdir = self.result_base_dir()
        self.results_header.setText(f"Sample: {item_text}")

        reads_dir = outdir / "figures" / "read_count" / "all"
        windows_dir = outdir / "figures" / "windows"
        raw_dir = outdir / "figures" / "raw" / "2Mb"
        acf_dir = outdir / "figures" / "acf"
        raw_kde_dir = outdir / "figures" / "raw"
        kde_dir = outdir / "figures" / "log2"

        self.col_reads["path"].setText(dir_with_trailing_slash(reads_dir))
        self.col_windows["path"].setText(dir_with_trailing_slash(windows_dir))
        self.col_acf["path"].setText(dir_with_trailing_slash(acf_dir))
        self.col_raw["path"].setText(dir_with_trailing_slash(raw_dir))
        self.col_rawkde["path"].setText(dir_with_trailing_slash(raw_kde_dir))
        self.col_kde["path"].setText(dir_with_trailing_slash(kde_dir))

        self.set_image_label(self.col_reads["e"], reads_dir / f"{desc}_{rep}_E.png")
        self.set_image_label(self.col_reads["l"], reads_dir / f"{desc}_{rep}_L.png")
        self.set_image_label(self.col_windows["e"], windows_dir / f"{desc}_{rep}_E.png")
        self.set_image_label(self.col_windows["l"], windows_dir / f"{desc}_{rep}_L.png")
        self.set_image_label(self.col_raw["img"], raw_dir / f"{desc}_{rep}_{bin_size}_{genome_id}_log2RT_2Mb_raw.png")
        self.set_image_label(self.col_acf["img"], acf_dir / f"{desc}_{rep}_{bin_size}_{genome_id}_log2RT.png")
        self.set_image_label(self.col_rawkde["img"], raw_kde_dir / f"{desc}_{rep}_{bin_size}_{genome_id}_log2RT_rawkdeplot.png")
        self.set_image_label(self.col_kde["img"], kde_dir / f"{desc}_{rep}_{bin_size}_{genome_id}_log2RT.png")
        self.update_results_visibility()


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
