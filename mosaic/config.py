"""
mosaic/config.py
Config schema, validation, and YAML I/O.
All fields and defaults match MOSAIC Design Document v0.7.
"""
from __future__ import annotations

import dataclasses as dc
import os
import sys
from pathlib import Path
from typing import Any, Optional

import yaml


def _expand_env(obj: Any) -> Any:
    """Recursively expand ``${VAR}``/``$VAR`` environment variables and a leading
    ``~`` in every string of a parsed-YAML structure. This lets config files use
    portable paths such as ``${MOSAIC_REF}/star_index`` so the same config works
    on any machine once the reviewer exports MOSAIC_REF. Strings without ``$`` or
    ``~`` are returned unchanged."""
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    if isinstance(obj, str):
        return os.path.expanduser(os.path.expandvars(obj))
    return obj


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

class ConfigError(ValueError):
    """Raised for config validation failures."""
    pass


def _require(value: Any, field: str, allowed: Optional[list] = None) -> None:
    """Assert a required field is present and optionally within allowed values."""
    if value is None:
        raise ConfigError(f"Required field '{field}' is missing.")
    if allowed and value not in allowed:
        raise ConfigError(
            f"Field '{field}' = {value!r} is not valid. "
            f"Allowed: {allowed}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Organism normalisation
# ─────────────────────────────────────────────────────────────────────────────

_ORGANISM_ALIASES: dict[str, str] = {
    "human": "human", "homo sapiens": "human", "9606": "human",
    "mouse": "mouse", "mus musculus": "mouse", "10090": "mouse",
    "rat": "rat", "rattus norvegicus": "rat", "10116": "rat",
    "zebrafish": "zebrafish", "danio rerio": "zebrafish", "7955": "zebrafish",
}

_GENOME_BUILD_DEFAULT: dict[str, str] = {
    "human": "hg38",
    "mouse": "mm39",
    "rat":   "rn7",
    "zebrafish": "GRCz11",
}

_GO_DB: dict[str, dict] = {
    "human":     {"r_db": "org.Hs.eg.db", "kegg": "hsa"},
    "mouse":     {"r_db": "org.Mm.eg.db", "kegg": "mmu"},
    "rat":       {"r_db": "org.Rn.eg.db", "kegg": "rno"},
    "zebrafish": {"r_db": "org.Dr.eg.db", "kegg": "dre"},
}


def normalise_organism(raw: str) -> str:
    key = str(raw).strip().lower()
    if key in _ORGANISM_ALIASES:
        return _ORGANISM_ALIASES[key]
    # unknown — return as-is, warn downstream
    return key


# ─────────────────────────────────────────────────────────────────────────────
# Sub-configs (dataclasses)
# ─────────────────────────────────────────────────────────────────────────────

@dc.dataclass
class ProjectConfig:
    name: str
    wd: Path
    organism: str                        # normalised
    genome_build: str
    random_seed: int = 42
    biological_context: str = "other"   # tumor|aging|congenital|other

    def validate(self) -> None:
        _require(self.name, "project.name")
        _require(self.organism, "project.organism")
        _require(self.biological_context, "project.biological_context",
                 ["tumor", "aging", "congenital", "other"])
        if not self.wd:
            raise ConfigError("project.wd must be specified.")


@dc.dataclass
class LibraryConfig:
    prep: str = "rRNA_depletion"
    layout: str = "PE"                  # PE|SE
    read_length: int = 150

    def validate(self) -> None:
        if self.prep.lower() != "rrna_depletion":
            raise ConfigError(
                f"library.prep = {self.prep!r} is not supported. "
                "MOSAIC v1.0 requires rRNA_depletion. "
                "polyA selection produces incorrect results for TE, satellite, "
                "circRNA and lncRNA — the primary targets of this pipeline."
            )
        _require(self.layout, "library.layout", ["PE", "SE"])


@dc.dataclass
class ReferenceDbConfig:
    path: Path = Path("/shared/reference_db")
    auto_update: bool = False


@dc.dataclass
class ModulesConfig:
    te_locus_specific: bool = True
    te_locus_tool: str = "telescope"    # telescope|squire
    circrna: bool = True
    alt_splicing: bool = True
    de_novo: bool = False               # controls --dta in HISAT2
    microbial: bool = False             # OPTIONAL module


@dc.dataclass
class TrimmingConfig:
    min_phred: int = 20
    min_length: int = 36
    min_length_warning: int = 75
    low_complexity_filter: bool = False  # MUST be False — TE/satellite safety
    pe_correction: bool = True
    trim_poly_x: bool = False            # MUST be False — satellite safety
    cut_right: bool = True               # sliding window 3' trimming
    cut_window_size: int = 4
    cut_mean_quality: int = 20
    disable_adapter_trimming: bool = False  # True for fusion-suspected tumors

    def validate(self) -> None:
        if self.low_complexity_filter:
            raise ConfigError(
                "trimming.low_complexity_filter must be False. "
                "LINE, SINE and satellite families have low complexity by definition — "
                "enabling this filter removes the target signal."
            )
        if self.trim_poly_x:
            raise ConfigError(
                "trimming.trim_poly_x must be False. "
                "Satellite sequences (TTAGGG, CATTC/GAATG) would be trimmed as poly-X."
            )
        if self.min_length < 20:
            raise ConfigError("trimming.min_length must be >= 20.")


@dc.dataclass
class MappingConfig:
    bowtie2_score_min: str = "L,0,1.6"
    hisat2_dta: str = "auto"          # auto|true|false
    hisat2_allow_discordant: bool = True  # False for normal samples

    def validate(self) -> None:
        _require(self.hisat2_dta, "mapping.hisat2_dta", ["auto", "true", "false"])


@dc.dataclass
class QualityConfig:
    strandedness_concordance_threshold: float = 0.90
    min_reads_warning: int = 5_000_000
    min_reads_error: int = 1_000_000
    bam_mapping_rate_warning: float = 30.0
    min_length_kraken: int = 75        # separate threshold for Kraken2

    def validate(self) -> None:
        if not 0.0 < self.strandedness_concordance_threshold <= 1.0:
            raise ConfigError(
                "quality.strandedness_concordance_threshold must be in (0, 1]."
            )
        if self.min_reads_error >= self.min_reads_warning:
            raise ConfigError(
                "quality.min_reads_error must be < quality.min_reads_warning."
            )


@dc.dataclass
class CircRNAConfig:
    min_bsj_reads: int = 5
    min_samples_per_group: str = "all"  # 'all' or integer as string
    de_test: str = "auto"               # auto|deseq2|wilcoxon


@dc.dataclass
class DEExtendedConfig:
    interactions: bool = False
    all_pairwise: bool = True
    multiple_testing: str = "bonferroni"  # bonferroni|IHW
    timeseries: bool = False
    timeseries_mode: str = "trend"         # trend|pairwise
    timeseries_splines_df: str = "auto"    # auto or integer
    timeseries_clustering: str = "mfuzz"   # mfuzz|kmeans|hierarchical


@dc.dataclass
class DEConfig:
    mode: str = "simple"                # simple|extended
    lfc_threshold: float = 1.0
    padj_threshold: float = 0.05
    lfc_shrink_method: str = "ashr"     # ashr|apeglm
    sizefactors_from: str = "canonical_stable"
    sizefactors_min_stable_genes: int = 500  # fallback threshold
    independent_filtering_te: bool = False   # OFF for TE/satellite
    extended: DEExtendedConfig = dc.field(default_factory=DEExtendedConfig)

    def validate(self) -> None:
        _require(self.mode, "de.mode", ["simple", "extended"])
        _require(self.lfc_shrink_method, "de.lfc_shrink_method", ["ashr", "apeglm"])
        if self.lfc_threshold < 0:
            raise ConfigError("de.lfc_threshold must be >= 0.")
        if not 0.0 < self.padj_threshold < 1.0:
            raise ConfigError("de.padj_threshold must be in (0, 1).")


@dc.dataclass
class TEConfig:
    go_universe: str = "te_host_expressed"
    batch_readlen_covariate: str = "auto"   # auto|true|false
    promoter_window_upstream: int = 2000



# ─────────────────────────────────────────────────────────────────────────────
# New v2 config sections
# ─────────────────────────────────────────────────────────────────────────────

@dc.dataclass
class DatabasesConfig:
    """Paths to all reference databases."""
    db_root:       str = ""
    genome_fa:     str = ""
    gtf:           str = ""
    star_index:    str = ""
    te_consensus:  str = ""
    te_satellite:  str = ""
    rmsk_bed:      str = ""
    te_gtf:        str = ""
    kraken2_db:    str = ""
    metaphlan_db:  str = ""
    chocophlan_db: str = ""
    uniref_db:     str = ""
    arriba_bin:    str = ""
    arriba_db:     str = ""
    ciri2_path:    str = "CIRI2.pl"
    auto_download: bool = True

    def get(self, key: str, default: str = "") -> str:
        return getattr(self, key, default) or default


@dc.dataclass
class ThresholdsConfig:
    """Analysis thresholds — production values."""
    padj:               float = 0.05
    lfc:                float = 1.0
    min_cov_editing:    int   = 5
    min_score_editing:  float = 1.5
    min_ef_editing:     float = 0.05
    editing_model:      str   = "DirMult:PHRED"
    min_reads_splicing: int   = 10
    min_clip_te:        int   = 30
    min_support_te:     int   = 2
    min_reads_meta:     int   = 10
    min_samples_meta:   int   = 2


@dc.dataclass
class EnvironmentsConfig:
    """Conda environment names for each tool group."""
    main:   str          = "mosaic"
    humann: str          = "humann3"
    arriba: Optional[str] = None   # None = use arriba_bin directly


@dc.dataclass
class SamplesConfig:
    """Sample metadata and FASTQ locations."""
    fastq_dir: str = ""
    coldata:   str = ""
    download:  list = dc.field(default_factory=list)

    # Populated at load time from coldata TSV
    sample_ids:  list = dc.field(default_factory=list, init=False)
    coldata_df:  Any  = dc.field(default=None, init=False)

    def load_coldata(self) -> None:
        """Load coldata TSV and extract sample_ids."""
        import csv
        p = Path(self.coldata)
        if not p.exists():
            return
        with open(p) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        self.sample_ids = [r["sample_id"] for r in rows]
        # Store as simple list of dicts (avoid pandas dependency)
        self.coldata_df = rows

    def get_condition_samples(self, condition: str) -> list:
        """Return sample_ids for a given condition."""
        if not self.coldata_df:
            return []
        return [r["sample_id"] for r in self.coldata_df
                if r.get("condition") == condition]


@dc.dataclass
class ModulesV2Config:
    """Per-module enable/disable flags (v2 schema)."""
    m00_download:      bool = False
    m02_trim:          bool = True
    m04_star:          bool = True
    m05_te_count:      bool = True
    m06_circrna:       bool = True
    m07_strandedness:  bool = True
    m08_deseq2:        bool = True
    m08b_qc:           bool = True
    m09_satellite:     bool = True
    m10_go:            bool = True
    m11_te_host:       bool = True
    m11b_promoter_te:  bool = True
    m11c_te_age:       bool = True
    m11d_te_de:        bool = True
    m12_tetranscripts: bool = True
    m15_fusion:        bool = True
    m16_stringtie:     bool = True
    m17_rmats:         bool = True
    m18_network:       bool = True
    m19_editing:       bool = True
    m20_meta:          bool = True
    m20b_humann:       bool = True
    m23_te_insertions: bool = True
    m14_report:        bool = True

    def enabled(self) -> list:
        """Return list of enabled module names."""
        return [k for k, v in vars(self).items() if v and k.startswith("m")]


# ─────────────────────────────────────────────────────────────────────────────
# Root config
# ─────────────────────────────────────────────────────────────────────────────

@dc.dataclass
class MosaicConfig:
    project:      ProjectConfig
    library:      LibraryConfig      = dc.field(default_factory=LibraryConfig)
    reference_db: ReferenceDbConfig  = dc.field(default_factory=ReferenceDbConfig)
    modules:      ModulesConfig      = dc.field(default_factory=ModulesConfig)
    trimming:     TrimmingConfig     = dc.field(default_factory=TrimmingConfig)
    mapping:      MappingConfig      = dc.field(default_factory=MappingConfig)
    quality:      QualityConfig      = dc.field(default_factory=QualityConfig)
    circrna:      CircRNAConfig      = dc.field(default_factory=CircRNAConfig)
    de:           DEConfig           = dc.field(default_factory=DEConfig)
    te:           TEConfig           = dc.field(default_factory=TEConfig)
    # v2 fields
    databases:    DatabasesConfig    = dc.field(default_factory=DatabasesConfig)
    thresholds:   ThresholdsConfig   = dc.field(default_factory=ThresholdsConfig)
    environments: EnvironmentsConfig = dc.field(default_factory=EnvironmentsConfig)
    samples:      SamplesConfig      = dc.field(default_factory=SamplesConfig)
    modules_v2:   ModulesV2Config    = dc.field(default_factory=ModulesV2Config)
    comparisons:  list               = dc.field(default_factory=list)
    threads:      int                = 20
    memory_gb:    int                = 30
    # Robustness: if True, a non-critical module failure is logged and the run
    # continues (backbone modules m02/m04/m05/m08 are always fatal on error).
    continue_on_error: bool          = True

    # Derived (set during validate)
    genome_build_resolved: str = dc.field(default="", init=False)
    go_db: dict = dc.field(default_factory=dict, init=False)

    def validate(self) -> None:
        """Run full validation. Raises ConfigError on first problem."""
        self.project.validate()
        self.library.validate()
        self.trimming.validate()
        self.mapping.validate()
        self.quality.validate()
        self.de.validate()

        # Resolve genome build
        if self.project.genome_build in ("auto", "", None):
            build = _GENOME_BUILD_DEFAULT.get(self.project.organism)
            if build is None:
                raise ConfigError(
                    f"No default genome build for organism {self.project.organism!r}. "
                    "Please set project.genome_build explicitly."
                )
            self.genome_build_resolved = build
        else:
            self.genome_build_resolved = self.project.genome_build

        # Resolve GO database
        self.go_db = _GO_DB.get(self.project.organism, {})
        if not self.go_db:
            # Not a known organism — warn, don't error
            self.go_db = {
                "r_db": "UNKNOWN — specify manually",
                "kegg": "UNKNOWN"
            }

        # apeglm only works cleanly with simple 2-group contrasts
        if (self.de.lfc_shrink_method == "apeglm"
                and self.de.mode == "extended"
                and self.de.extended.all_pairwise):
            raise ConfigError(
                "de.lfc_shrink_method='apeglm' is incompatible with "
                "de.extended.all_pairwise=True (apeglm requires a single "
                "contrast as the first model coefficient). "
                "Use lfc_shrink_method='ashr' for multi-contrast designs."
            )

        # hisat2_dta resolution
        if self.mapping.hisat2_dta == "auto":
            self._hisat2_dta_resolved = self.modules.de_novo
        else:
            self._hisat2_dta_resolved = self.mapping.hisat2_dta == "true"

    @property
    def hisat2_dta(self) -> bool:
        return getattr(self, "_hisat2_dta_resolved", False)

    @property
    def is_tumor(self) -> bool:
        return self.project.biological_context == "tumor"

    @property
    def go_r_db(self) -> str:
        return self.go_db.get("r_db", "UNKNOWN")

    @property
    def go_kegg(self) -> str:
        return self.go_db.get("kegg", "UNKNOWN")


# ─────────────────────────────────────────────────────────────────────────────
# YAML I/O
# ─────────────────────────────────────────────────────────────────────────────

def _deep_get(d: dict, *keys, default=None) -> Any:
    """Safe nested dict access."""
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k, default)
    return d


def load_config(path: str | Path) -> MosaicConfig:
    """Load, parse and validate a MOSAIC config YAML file."""
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"Config file not found: {p}")

    with open(p) as fh:
        raw = yaml.safe_load(fh)

    if not isinstance(raw, dict):
        raise ConfigError("Config YAML must be a mapping at the top level.")

    # Expand ${ENV} / ~ in all string values for portable, machine-independent paths.
    raw = _expand_env(raw)

    proj_raw = raw.get("project", {})
    organism = normalise_organism(proj_raw.get("organism", ""))
    genome_build = proj_raw.get("genome_build", "auto")
    if genome_build == "auto":
        genome_build = _GENOME_BUILD_DEFAULT.get(organism, "auto")

    project = ProjectConfig(
        name=proj_raw.get("name", ""),
        wd=Path(proj_raw.get("wd", ".")),
        organism=organism,
        genome_build=genome_build,
        random_seed=proj_raw.get("random_seed", 42),
        biological_context=proj_raw.get("biological_context", "other"),
    )

    lib_raw = raw.get("library", {})
    library = LibraryConfig(
        prep=lib_raw.get("prep", "rRNA_depletion"),
        layout=lib_raw.get("layout", "PE"),
    )

    ref_raw = raw.get("reference_db", {})
    reference_db = ReferenceDbConfig(
        path=Path(ref_raw.get("path", "/shared/reference_db")),
        auto_update=ref_raw.get("auto_update", False),
    )

    mod_raw = raw.get("modules", {})
    modules = ModulesConfig(
        te_locus_specific=mod_raw.get("te_locus_specific", True),
        te_locus_tool=mod_raw.get("te_locus_tool", "telescope"),
        circrna=mod_raw.get("circrna", True),
        alt_splicing=mod_raw.get("alt_splicing", True),
        de_novo=mod_raw.get("de_novo", False),
        microbial=mod_raw.get("microbial", False),
    )

    tr_raw = raw.get("trimming", {})
    trimming = TrimmingConfig(
        min_phred=tr_raw.get("min_phred", 20),
        min_length=tr_raw.get("min_length", 36),
        min_length_warning=tr_raw.get("min_length_warning", 75),
        low_complexity_filter=tr_raw.get("low_complexity_filter", False),
        pe_correction=tr_raw.get("pe_correction", True),
        trim_poly_x=tr_raw.get("trim_poly_x", False),
        cut_right=tr_raw.get("cut_right", True),
        cut_window_size=tr_raw.get("cut_window_size", 4),
        cut_mean_quality=tr_raw.get("cut_mean_quality", 20),
        disable_adapter_trimming=tr_raw.get("disable_adapter_trimming", False),
    )

    map_raw = raw.get("mapping", {})
    mapping = MappingConfig(
        bowtie2_score_min=map_raw.get("bowtie2_score_min", "L,0,1.6"),
        hisat2_dta=str(map_raw.get("hisat2_dta", "auto")),
        hisat2_allow_discordant=map_raw.get("hisat2_allow_discordant", True),
    )

    q_raw = raw.get("quality", {})
    quality = QualityConfig(
        strandedness_concordance_threshold=q_raw.get("strandedness_concordance_threshold", 0.90),
        min_reads_warning=q_raw.get("min_reads_warning", 5_000_000),
        min_reads_error=q_raw.get("min_reads_error", 1_000_000),
        bam_mapping_rate_warning=q_raw.get("bam_mapping_rate_warning", 30.0),
        min_length_kraken=q_raw.get("min_length_kraken", 75),
    )

    circ_raw = raw.get("circrna", {})
    circrna = CircRNAConfig(
        min_bsj_reads=circ_raw.get("min_bsj_reads", 5),
        min_samples_per_group=str(circ_raw.get("min_samples_per_group", "all")),
        de_test=circ_raw.get("de_test", "auto"),
    )

    de_raw = raw.get("de", {})
    ext_raw = de_raw.get("extended", {})
    de = DEConfig(
        mode=de_raw.get("mode", "simple"),
        lfc_threshold=float(de_raw.get("lfc_threshold", 1.0)),
        padj_threshold=float(de_raw.get("padj_threshold", 0.05)),
        lfc_shrink_method=de_raw.get("lfc_shrink_method", "ashr"),
        sizefactors_from=de_raw.get("sizefactors_from", "canonical_stable"),
        sizefactors_min_stable_genes=int(de_raw.get("sizefactors_min_stable_genes", 500)),
        independent_filtering_te=de_raw.get("independent_filtering_te", False),
        extended=DEExtendedConfig(
            interactions=ext_raw.get("interactions", False),
            all_pairwise=ext_raw.get("all_pairwise", True),
            multiple_testing=ext_raw.get("multiple_testing", "bonferroni"),
            timeseries=ext_raw.get("timeseries", False),
            timeseries_mode=ext_raw.get("timeseries_mode", "trend"),
            timeseries_splines_df=str(ext_raw.get("timeseries_splines_df", "auto")),
            timeseries_clustering=ext_raw.get("timeseries_clustering", "mfuzz"),
        ),
    )

    te_raw = raw.get("te", {})
    te = TEConfig(
        go_universe=te_raw.get("go_universe", "te_host_expressed"),
        batch_readlen_covariate=str(te_raw.get("batch_readlen_covariate", "auto")),
        promoter_window_upstream=int(te_raw.get("promoter_window_upstream", 2000)),
    )

    # ── v2 fields ────────────────────────────────────────────────────────────
    db_raw = raw.get("databases", {})
    databases = DatabasesConfig(
        db_root=db_raw.get("db_root", ""),
        genome_fa=db_raw.get("genome_fa", ""),
        gtf=db_raw.get("gtf", ""),
        star_index=db_raw.get("star_index", ""),
        te_consensus=db_raw.get("te_consensus", ""),
        te_satellite=db_raw.get("te_satellite", ""),
        rmsk_bed=db_raw.get("rmsk_bed", ""),
        te_gtf=db_raw.get("te_gtf", ""),
        kraken2_db=db_raw.get("kraken2_db", ""),
        metaphlan_db=db_raw.get("metaphlan_db", ""),
        chocophlan_db=db_raw.get("chocophlan_db", ""),
        uniref_db=db_raw.get("uniref_db", ""),
        arriba_bin=db_raw.get("arriba_bin", ""),
        arriba_db=db_raw.get("arriba_db", ""),
        ciri2_path=db_raw.get("ciri2_path", "CIRI2.pl"),
        auto_download=db_raw.get("auto_download", True),
    )

    thr_raw = raw.get("thresholds", {})
    thresholds = ThresholdsConfig(
        padj=float(thr_raw.get("padj", 0.05)),
        lfc=float(thr_raw.get("lfc", 1.0)),
        min_cov_editing=int(thr_raw.get("min_cov_editing", 5)),
        min_score_editing=float(thr_raw.get("min_score_editing", 1.5)),
        min_ef_editing=float(thr_raw.get("min_ef_editing", 0.05)),
        editing_model=thr_raw.get("editing_model", "DirMult:PHRED"),
        min_reads_splicing=int(thr_raw.get("min_reads_splicing", 10)),
        min_clip_te=int(thr_raw.get("min_clip_te", 30)),
        min_support_te=int(thr_raw.get("min_support_te", 2)),
        min_reads_meta=int(thr_raw.get("min_reads_meta", 10)),
        min_samples_meta=int(thr_raw.get("min_samples_meta", 2)),
    )

    env_raw = raw.get("environments", {})
    environments = EnvironmentsConfig(
        main=env_raw.get("main", "mosaic"),
        humann=env_raw.get("humann", "humann3"),
        arriba=env_raw.get("arriba", None),
    )

    smp_raw = raw.get("samples", {})
    samples = SamplesConfig(
        fastq_dir=smp_raw.get("fastq_dir", ""),
        coldata=smp_raw.get("coldata", ""),
        download=smp_raw.get("download", []),
    )
    samples.load_coldata()

    mod2_raw = raw.get("modules", {})
    # Support both old (te_locus_specific etc.) and new (m08_deseq2 etc.) keys
    modules_v2 = ModulesV2Config(
        m00_download=mod2_raw.get("m00_download", False),
        m02_trim=mod2_raw.get("m02_trim", True),
        m04_star=mod2_raw.get("m04_star", True),
        m05_te_count=mod2_raw.get("m05_te_count", True),
        m06_circrna=mod2_raw.get("m06_circrna", mod2_raw.get("circrna", True)),
        m07_strandedness=mod2_raw.get("m07_strandedness", True),
        m08_deseq2=mod2_raw.get("m08_deseq2", True),
        m08b_qc=mod2_raw.get("m08b_qc", True),
        m09_satellite=mod2_raw.get("m09_satellite", True),
        m10_go=mod2_raw.get("m10_go", True),
        m11_te_host=mod2_raw.get("m11_te_host", True),
        m11b_promoter_te=mod2_raw.get("m11b_promoter_te", True),
        m11c_te_age=mod2_raw.get("m11c_te_age", True),
        m11d_te_de=mod2_raw.get("m11d_te_de", True),
        m12_tetranscripts=mod2_raw.get("m12_tetranscripts",
                                        mod2_raw.get("te_locus_specific", True)),
        m15_fusion=mod2_raw.get("m15_fusion", True),
        m16_stringtie=mod2_raw.get("m16_stringtie",
                                    mod2_raw.get("de_novo", True)),
        m17_rmats=mod2_raw.get("m17_rmats",
                                mod2_raw.get("alt_splicing", True)),
        m18_network=mod2_raw.get("m18_network", True),
        m19_editing=mod2_raw.get("m19_editing", True),
        m20_meta=mod2_raw.get("m20_meta",
                               mod2_raw.get("microbial", True)),
        m20b_humann=mod2_raw.get("m20b_humann", True),
        m23_te_insertions=mod2_raw.get("m23_te_insertions", True),
        m14_report=mod2_raw.get("m14_report", True),
    )

    comparisons = raw.get("comparisons", [])
    # Fallback: build from project biological_context
    if not comparisons and project.biological_context == "tumor":
        comparisons = [{"treatment": "tumor", "control": "normal"}]

    threads   = int(raw.get("threads", 20))
    memory_gb = int(raw.get("memory_gb", 30))

    # Also read sequencing section
    seq_raw = raw.get("sequencing", {})
    if seq_raw.get("type") and not lib_raw.get("layout"):
        library.layout = seq_raw["type"]
    if seq_raw.get("strandedness"):
        library._strandedness_config = seq_raw["strandedness"]
    if seq_raw.get("read_length"):
        library.read_length = int(seq_raw["read_length"])

    cfg = MosaicConfig(
        project=project,
        library=library,
        reference_db=reference_db,
        modules=modules,
        trimming=trimming,
        mapping=mapping,
        quality=quality,
        circrna=circrna,
        de=de,
        te=te,
        databases=databases,
        thresholds=thresholds,
        environments=environments,
        samples=samples,
        modules_v2=modules_v2,
        comparisons=comparisons,
        threads=threads,
        memory_gb=memory_gb,
        continue_on_error=bool(raw.get("continue_on_error", True)),
    )
    cfg.validate()
    return cfg


def generate_default_config(organism: str = "human",
                            context: str = "tumor",
                            layout: str = "PE") -> str:
    """Return a fully documented default YAML config string."""
    org = normalise_organism(organism)
    build = _GENOME_BUILD_DEFAULT.get(org, "hg38")
    return f"""\
# MOSAIC config — generated by: mosaic init
# Edit paths and parameters, then run: mosaic check config.yaml

project:
  name: my_project
  wd: /path/to/workdir          # output directory
  organism: {org}               # human | mouse | rat | zebrafish
  genome_build: {build}         # hg38 | hg19 | mm39 | mm10 | rn7
  annotation_version: gencode_v49
  biological_context: {context} # tumor | aging | congenital | other
  random_seed: 42

sequencing:
  type: {layout}                # PE | SE
  read_length: 150
  strandedness: auto            # auto | RF | FR | unstranded

library:
  prep: rRNA_depletion          # rRNA_depletion | polyA (polyA blocks TE)
  layout: {layout}              # PE | SE

samples:
  fastq_dir: /path/to/fastq
  coldata: /path/to/coldata.tsv # columns: sample_id, condition [, batch]
  # download:                   # optional: auto-download FASTQs
  #   - sample_id: SRR19388325
  #     source: sra             # ena | sra | url | ssh | local
  #     accession: SRR19388325

comparisons:
  - treatment: tumor
    control: normal
  # - treatment: treated
  #   control: untreated

databases:
  db_root: /path/to/references
  genome_fa:     /path/to/references/genome/genome.fa
  gtf:           /path/to/references/annotation.gtf
  star_index:    /path/to/references/star_index
  te_consensus:  /path/to/references/dfam/pseudochromosomes/consensus_TE.fa
  te_satellite:  /path/to/references/dfam/pseudochromosomes/consensus_satellite.fa
  rmsk_bed:      /path/to/references/rmsk/rmsk.txt.gz
  te_gtf:        /path/to/references/rmsk/TE.gtf
  kraken2_db:    /path/to/references/kraken2_db
  metaphlan_db:  /path/to/references/humann3_db/metaphlan_db
  chocophlan_db: /path/to/references/humann3_db/chocophlan
  uniref_db:     /path/to/references/humann3_db/uniref
  arriba_bin:    /path/to/references/arriba_v2.4.0/arriba
  arriba_db:     /path/to/references/arriba_v2.4.0/database
  auto_download: true

environments:
  main:   mosaic
  humann: humann3
  arriba: null                  # null = use arriba_bin directly

thresholds:
  padj:               0.05
  lfc:                1.0
  min_cov_editing:    5
  min_score_editing:  1.5
  min_ef_editing:     0.05
  editing_model:      "DirMult:PHRED"
  min_reads_splicing: 10
  min_clip_te:        30
  min_support_te:     2
  min_reads_meta:     10
  min_samples_meta:   2

modules:
  m00_download:      false
  m02_trim:          true
  m04_star:          true
  m05_te_count:      true
  m06_circrna:       true
  m07_strandedness:  true
  m08_deseq2:        true
  m09_satellite:     true
  m10_go:            true
  m11_te_host:       true
  m11b_promoter_te:  true
  m11c_te_age:       true
  m11d_te_de:        true
  m12_tetranscripts: true
  m15_fusion:        true
  m16_stringtie:     true
  m17_rmats:         true
  m18_network:       true
  m19_editing:       true
  m20_meta:          true
  m20b_humann:       true
  m23_te_insertions: true
  m14_report:        true

threads: 20
memory_gb: 30

trimming:
  min_phred: 20
  min_length: 36
  low_complexity_filter: false  # NEVER true — removes TE reads
  trim_poly_x: false            # NEVER true — removes satellite reads

report:
  title: "MOSAIC Analysis Report"
  remove_test_notes: true
  produce_pdf: false
"""
