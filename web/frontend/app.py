"""
ROMARINOVA Scientific Platform — Streamlit app.

Multi-page interface for exploring the ROMARINOVA phytochemical screening
pipeline: compound atlas, α-glucosidase shortlist, biological target viewer,
methodology, evidence explorer, and a biodiversity map of Tunisian rosemary.

Requires Streamlit >= 1.40 (uses st.navigation, st.Page, st.fragment).
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable
import os
import re
import base64
import html as _html
import logging
import traceback
import uuid as _uuid
from urllib.parse import quote
from io import BytesIO
from datetime import datetime
import plotly
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import json
import plotly.express as px
import plotly.graph_objects as go
import pydeck as pdk
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import py3Dmol

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("romarinova")

# ==========================================================
# 1. CONFIGURATION & DATA LOADING
# ==========================================================

st.set_page_config(
    page_title="ROMARINOVA | Scientific Platform",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

BASE_DIR = Path(__file__).resolve().parent
ASSETS_DIR = BASE_DIR / "assets"
PROJECT_ROOT = Path(os.environ.get("ROMARINOVA_ROOT", BASE_DIR.parent.parent)).resolve()
RESULTS_DIR = PROJECT_ROOT / "results"
MODELS_DIR = PROJECT_ROOT / "models"

CONTACT_URL = os.environ.get(
    "ROMARINOVA_CONTACT", "https://github.com/GreenDevilz/Romarinova-AI"
)

try:
    BUILD_DATE = datetime.fromtimestamp(Path(__file__).stat().st_mtime).strftime("%Y-%m-%d")
except OSError:
    BUILD_DATE = "unknown"

# ---- Domain constants (validated at load time by _validate_enums) -----------
AD_WITHIN = "within_AD"
AD_EXTRAPOLATION = "extrapolation"

DIAG_SENSITIVITY = "sensitivity"
DIAG_FILTERED = "filtered_but_interesting"
DIAG_DOCKING = "docking_failures"
DIAG_MISSING = "missing_structures"
DIAG_DISAGREEMENT = "disagreement_patterns"

DATA_CONTRACT = {
    "screened": [
        "compound_id", "compound_name", "canonical_smiles", "inchikey",
        "ML_rank_score", "similarity_support",
        "docking_affinity_kcal_per_mol", "docking_missing",
        "ad_flag", "structure_status", "known_active_flag",
        "developability_score",
    ],
    "pair_sl": [
        "compound_id", "compound_name",
        "final_rank", "ligand_only_rank",
        "final_score", "ligand_only_score",
        "ML_rank_score", "similarity_support", "docking_score",
    ],
    "comp_sl": [
        "compound_id", "compound_name",
        "final_rank", "ligand_only_rank",
        "final_score", "ligand_only_score",
        "ML_rank_score", "similarity_support", "docking_score",
    ],
    "eval_table": ["model", "split", "pr_auc"],
    "diag": ["section"],
}

def _parse_st_version(v: str) -> tuple[int, int, int]:
    m = re.match(r"(\d+)\.(\d+)\.?(\d+)?", str(v))
    if not m:
        return (0, 0, 0)
    return tuple(int(g) if g else 0 for g in m.groups()[:3])  # type: ignore[return-value]

_USE_WIDTH_STRETCH = _parse_st_version(st.__version__) >= (1, 49, 0)

_FULL_WIDTH_KW = (
    {"width": "stretch"} if _USE_WIDTH_STRETCH else {"use_container_width": True}
)
_DATAFRAME_WIDTH_KW = _FULL_WIDTH_KW
_CHART_WIDTH_KW = _FULL_WIDTH_KW
_BUTTON_WIDTH_KW = _FULL_WIDTH_KW

_BOOL_COLS: dict[str, bool] = {
    "docking_missing": True,
    "known_active_flag": False,
}

_DEFAULT_WEIGHTS = {
    "ML": 0.30, "similarity": 0.20, "developability": 0.20,
    "novelty": 0.10, "docking": 0.20,
}

_FORCE_OVERRIDE_KEYS = {
    "N_star", "U_star", "ml_rank_top_band_count", "ml_rank_top_band_value",
}


def to_bool(series: pd.Series, missing: bool = False) -> pd.Series:
    if series.dtype == bool:
        return series
    if pd.api.types.is_numeric_dtype(series):
        def _num(v):
            if pd.isna(v):
                return missing
            try:
                fv = float(v)
            except (TypeError, ValueError):
                return missing
            if fv == 1.0:
                return True
            if fv == 0.0:
                return False
            return missing
        return series.map(_num).astype(bool)

    _TRUE_STRINGS = {"true", "1", "1.0", "yes"}
    _FALSE_STRINGS = {"false", "0", "0.0", "no"}
    s = series.astype(str).str.strip().str.lower()
    confirmed = s.isin(_TRUE_STRINGS | _FALSE_STRINGS)
    values = s.isin(_TRUE_STRINGS)
    return values.where(confirmed, missing).astype(bool)


def _to_int_or_none(v) -> int | None:
    if v is None:
        return None
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _is_missing(v) -> bool:
    if v is None or v is pd.NA:
        return True
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def _unique_screened() -> pd.DataFrame:
    """Return `screened` deduplicated to one row per compound."""
    if "compound_id" not in screened.columns:
        return screened
    return screened.dropna(subset=["compound_id"]).drop_duplicates("compound_id")


def _count_extrapolated_unique_compounds() -> int:
    """Number of *unique* compounds flagged as extrapolation."""
    if "ad_flag" not in screened.columns:
        return 0
    _uniq = _unique_screened()
    return int((_uniq["ad_flag"].astype(str).str.strip() == AD_EXTRAPOLATION).sum())


def _derive_stats(screened, pair_sl, comp_sl, diag, meta) -> tuple[dict, list[str]]:
    out: dict = {}
    consistency: list[str] = []

    _uniq_screened = (
        screened.dropna(subset=["compound_id"]).drop_duplicates("compound_id")
        if "compound_id" in screened.columns else screened
    )

    if "compound_id" in screened.columns:
        out["screened_compound_n"] = int(screened["compound_id"].nunique())

    n_star_eff: int | None = None
    if "compound_id" in pair_sl.columns:
        n_star_eff = int(len(pair_sl))
        _meta_n = _to_int_or_none(meta.get("N_star"))
        if _meta_n is None or _meta_n != n_star_eff:
            out["N_star"] = n_star_eff
            if _meta_n is not None:
                consistency.append(
                    f"metadata N_star={_meta_n} but primary_shortlist.csv has "
                    f"{n_star_eff} rows — using {n_star_eff}."
                )
    else:
        n_star_eff = _to_int_or_none(meta.get("N_star"))

    if "compound_id" in comp_sl.columns:
        _actual_u = int(comp_sl["compound_id"].nunique())
        _meta_u = _to_int_or_none(meta.get("U_star"))
        if _meta_u is None or _meta_u != _actual_u:
            out["U_star"] = _actual_u
            if _meta_u is not None:
                consistency.append(
                    f"metadata U_star={_meta_u} but primary_shortlist_compound_level.csv "
                    f"has {_actual_u} unique compounds — using {_actual_u}."
                )

    if "ML_rank_score" in _uniq_screened.columns and "compound_id" in _uniq_screened.columns:
        vals = pd.to_numeric(_uniq_screened["ML_rank_score"], errors="coerce").dropna()
        if not vals.empty:
            max_val = round(float(vals.max()), 6)
            top_count = int((vals.round(6) == max_val).sum())
            out["ml_rank_top_band_value"] = max_val
            out["ml_rank_top_band_count"] = top_count

    if "section" in diag.columns:
        sens = diag[diag["section"] == DIAG_SENSITIVITY]
        if not sens.empty:
            row = sens.iloc[0]
            if "pair_overlap" in diag.columns:
                po_i = _to_int_or_none(row.get("pair_overlap"))
                if po_i is not None and n_star_eff is not None and n_star_eff > 0:
                    out["pair_overlap"] = f"{po_i}/{n_star_eff}"
            sl = row.get("sensitivity_label")
            if not _is_missing(sl) and str(sl).strip():
                out["sensitivity_label"] = str(sl).strip()

    ad_missing = None
    if "ad_flag" in _uniq_screened.columns:
        ad_missing = int(
            (_uniq_screened["ad_flag"].astype(str).str.strip() != AD_WITHIN).sum()
        )
    dock_missing = None
    if "docking_missing" in _uniq_screened.columns:
        dock_missing = int(
            _uniq_screened["docking_missing"].fillna(True).astype(bool).sum()
        )
    if ad_missing is not None and dock_missing is not None:
        if ad_missing == 0 and dock_missing == 0:
            out["coverage_gap"] = "None"
        else:
            parts = []
            if ad_missing:
                parts.append(f"{ad_missing} AD")
            if dock_missing:
                parts.append(f"{dock_missing} docking")
            out["coverage_gap"] = ", ".join(parts) + " missing"

    return out, consistency


@st.cache_data(ttl=None, show_spinner="Loading pipeline artefacts…", max_entries=1)
def load_data():
    screened = pd.read_csv(RESULTS_DIR / "screened_compounds_annotated.csv", encoding="utf-8")
    pair_sl = pd.read_csv(RESULTS_DIR / "primary_shortlist.csv", encoding="utf-8")
    comp_sl = pd.read_csv(RESULTS_DIR / "primary_shortlist_compound_level.csv", encoding="utf-8")
    diag = pd.read_csv(RESULTS_DIR / "diagnostics.csv", encoding="utf-8")
    meta = json.loads((RESULTS_DIR / "shortlist_metadata.json").read_text(encoding="utf-8"))
    eval_table = pd.read_csv(RESULTS_DIR / "evaluation_table.csv", encoding="utf-8")
    dock_cfg = json.loads((MODELS_DIR / "docking_config.json").read_text(encoding="utf-8"))

    for df in (screened, pair_sl, comp_sl):
        for f, default in _BOOL_COLS.items():
            if f in df.columns and df[f].dtype != bool:
                df[f] = to_bool(df[f], missing=default)

    derived, consistency = _derive_stats(screened, pair_sl, comp_sl, diag, meta)
    filled_from_derivation: list[str] = []
    for k, v in derived.items():
        old = meta.get(k)
        was_missing = k not in meta or old in (None, "", "—")
        if k in _FORCE_OVERRIDE_KEYS:
            if old != v:
                filled_from_derivation.append(k)
            meta[k] = v
        elif was_missing:
            meta[k] = v
            filled_from_derivation.append(k)

    return (
        screened, pair_sl, comp_sl, diag, meta, eval_table, dock_cfg,
        filled_from_derivation, consistency,
    )


def _load_or_stop():
    try:
        return load_data()
    except FileNotFoundError as e:
        st.error(f"Missing data file: `{e}`. Ensure the pipeline has been run.")
        st.stop()
    except Exception as e:
        st.error(f"Error loading data: `{e}`")
        log.exception("load_data failed")
        with st.expander("Traceback"):
            st.code(traceback.format_exc())
        st.stop()


(
    screened, pair_sl, comp_sl, diag, meta, eval_table, dock_cfg,
    _filled_from_derivation, _CONSISTENCY_WARNINGS,
) = _load_or_stop()


def validate_columns(name: str, df: pd.DataFrame, required: list[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        st.error(f"Data contract violated in `{name}` — missing columns: {missing}")
        st.stop()


def select_existing(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [c for c in cols if c in df.columns]


for _name, _df in [
    ("screened", screened), ("pair_sl", pair_sl), ("comp_sl", comp_sl),
    ("diag", diag), ("eval_table", eval_table),
]:
    validate_columns(_name, _df, DATA_CONTRACT[_name])


def _validate_enums() -> list[str]:
    issues: list[str] = []
    if "ad_flag" in screened.columns:
        vals = set(screened["ad_flag"].dropna().astype(str).str.strip().unique())
        bad = vals - {AD_WITHIN, AD_EXTRAPOLATION}
        if bad:
            issues.append(f"Unexpected `ad_flag` values: {sorted(bad)}")
    return issues


_ENUM_ISSUES = _validate_enums()


def _resolve_test_pr_auc(eval_table: pd.DataFrame, meta: dict) -> tuple[float | None, str | None]:
    preferred = meta.get("selected_model") or "Selected RF"
    row = eval_table[(eval_table["model"] == preferred) & (eval_table["split"] == "test")]
    if not row.empty:
        return float(row.iloc[0]["pr_auc"]), preferred
    test_rows = eval_table[eval_table["split"] == "test"]
    if not test_rows.empty:
        best = test_rows.loc[test_rows["pr_auc"].idxmax()]
        return float(best["pr_auc"]), str(best["model"])
    return None, None


_RF_PR_AUC, _RF_MODEL = _resolve_test_pr_auc(eval_table, meta)


@st.cache_data(show_spinner=False)
def _artefact_mtime() -> str | None:
    paths = [
        RESULTS_DIR / "screened_compounds_annotated.csv",
        RESULTS_DIR / "primary_shortlist.csv",
        RESULTS_DIR / "shortlist_metadata.json",
    ]
    mtimes = []
    for p in paths:
        try:
            mtimes.append(p.stat().st_mtime)
        except OSError:
            continue
    if not mtimes:
        return None
    return datetime.fromtimestamp(max(mtimes)).strftime("%Y-%m-%d %H:%M")


# ==========================================================
# 2. HELPERS
# ==========================================================

CHART_SMALL = 300
CHART_MEDIUM = 400
CHART_LARGE = 480

LEGEND_TOP_RIGHT = dict(
    orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
)


def esc(value) -> str:
    if value is None:
        return ""
    return _html.escape(str(value))


_BOLD_STAR_RE = re.compile(r"\*\*(?!\s)(.+?)(?<!\s)\*\*")
_CODE_RE = re.compile(r"`([^`]+?)`")


def safe_md(text) -> str:
    s = _html.escape("" if text is None else str(text))
    stash: list[str] = []
    token = f"\x00ROM_{_uuid.uuid4().hex}_"
    token_re = re.compile(re.escape(token) + r"(\d+)\x00")

    def _stash(m):
        stash.append(m.group(1))
        return f"{token}{len(stash) - 1}\x00"

    s = _CODE_RE.sub(_stash, s)
    s = _BOLD_STAR_RE.sub(r"<strong>\1</strong>", s)
    s = token_re.sub(lambda m: f"<code>{stash[int(m.group(1))]}</code>", s)
    return s


def meta_val(key, default="—"):
    v = meta.get(key)
    if _is_missing(v):
        return default
    if isinstance(v, str) and (not v.strip() or v == "—"):
        return default
    return v


def row_get(row, key, default="—"):
    try:
        if key not in row.index:
            return default
    except AttributeError:
        return default
    v = row[key]
    if _is_missing(v):
        return default
    if isinstance(v, str) and not v.strip():
        return default
    return v


def fmt(value, digits: int = 3, dash: str = "—") -> str:
    if value is None:
        return dash
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if pd.isna(f):
        return dash
    return f"{f:.{digits}f}"


@st.cache_data(show_spinner=False, max_entries=1024)
def mol_image(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    img = Image.open(BytesIO(drawer.GetDrawingText()))
    img.load()
    return img


@st.cache_data(show_spinner=False, max_entries=8)
def _compute_descriptors(smiles_tuple: tuple[str, ...]) -> list[dict]:
    from rdkit.Chem import Descriptors, Lipinski
    out: list[dict] = []
    for smi in smiles_tuple:
        if not isinstance(smi, str) or not smi.strip():
            out.append({"MW": None, "LogP": None, "HBD": None, "HBA": None})
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            out.append({"MW": None, "LogP": None, "HBD": None, "HBA": None})
            continue
        try:
            out.append({
                "MW": float(Descriptors.MolWt(m)),
                "LogP": float(Descriptors.MolLogP(m)),
                "HBD": int(Lipinski.NumHDonors(m)),
                "HBA": int(Lipinski.NumHAcceptors(m)),
            })
        except Exception:
            log.exception("descriptor computation failed for %s", smi)
            out.append({"MW": None, "LogP": None, "HBD": None, "HBA": None})
    return out


@st.cache_data(show_spinner=False, max_entries=32)
def render_pdb_html(
    pdb_path: str, cx: float, cy: float, cz: float,
    width: int = 600, height: int = 450,
    show_surface: bool = False, show_site: bool = True,
    pdb_mtime: float = 0.0,
) -> str | None:
    p = Path(pdb_path)
    if not p.exists():
        return None
    pdb_str = p.read_text(encoding="utf-8")
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_str, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})

    if show_surface:
        try:
            if show_site:
                view.addSurface(
                    py3Dmol.VDW,
                    {"opacity": 0.20, "color": "white", "depthWrite": False},
                )
            else:
                view.addSurface(
                    py3Dmol.VDW,
                    {"opacity": 0.55, "color": "white"},
                )
        except Exception:
            log.exception("py3Dmol surface failed")

    if show_site:
        view.addSphere({
            "center": {"x": cx, "y": cy, "z": cz},
            "radius": 6,
            "color": "yellow",
            "opacity": 1.0,
        })
    view.zoomTo()

    for meth in ("_make_html", "to_html", "_repr_html_"):
        fn = getattr(view, meth, None)
        if fn is None:
            continue
        try:
            result = fn()
        except Exception:
            log.exception("py3Dmol.%s failed", meth)
            continue
        if isinstance(result, str) and result:
            return result

    raise RuntimeError("py3Dmol API changed.")


def parse_grid_center(center) -> tuple[float, float, float]:
    if isinstance(center, dict):
        try:
            return (float(center.get("x", 0.0)), float(center.get("y", 0.0)),
                    float(center.get("z", 0.0)))
        except (TypeError, ValueError):
            return (0.0, 0.0, 0.0)
    if isinstance(center, (list, tuple)) and len(center) >= 3:
        try:
            return (float(center[0]), float(center[1]), float(center[2]))
        except (TypeError, ValueError):
            return (0.0, 0.0, 0.0)
    return (0.0, 0.0, 0.0)


def best_docking_row(rows: pd.DataFrame) -> tuple[float | None, str | None]:
    if not {"docking_missing", "docking_affinity_kcal_per_mol"}.issubset(rows.columns):
        return (None, None)
    valid = rows.loc[~rows["docking_missing"].fillna(True).astype(bool)].dropna(
        subset=["docking_affinity_kcal_per_mol"])
    if valid.empty:
        return (None, None)
    idx = valid["docking_affinity_kcal_per_mol"].idxmin()
    affinity = float(valid.loc[idx, "docking_affinity_kcal_per_mol"])
    variety = None
    if "variety" in valid.columns:
        v = valid.loc[idx, "variety"]
        if not _is_missing(v) and isinstance(v, str) and v.strip():
            variety = v.strip()
    return (affinity, variety)


_LOGO_NAMES = (
    "logo.png", "logo.jpg", "logo.jpeg", "logo.webp",
    "romarinova.png", "romarinova.jpg", "romarinova.jpeg", "romarinova.webp",
)


def find_logo():
    if not ASSETS_DIR.exists():
        return None
    for name in _LOGO_NAMES:
        candidate = ASSETS_DIR / name
        if candidate.exists():
            return candidate
    for extension in ("*.png", "*.jpg", "*.jpeg", "*.webp"):
        files = sorted(ASSETS_DIR.glob(extension))
        if files:
            return files[0]
    return None


@st.cache_data(show_spinner=False, max_entries=4)
def _logo_data_uri(path_str: str) -> str | None:
    p = Path(path_str)
    if not p.exists():
        return None
    try:
        raw = p.read_bytes()
    except OSError:
        log.exception("logo read failed: %s", path_str)
        return None
    ext = p.suffix.lstrip(".").lower().replace("jpg", "jpeg") or "png"
    return f"data:image/{ext};base64,{base64.b64encode(raw).decode('ascii')}"


_PLOTLY_BASE = dict(
    template="plotly_white",
    plot_bgcolor="white",
    paper_bgcolor="white",
    font=dict(color="#24332C", family="DM Sans, sans-serif", size=13),
    title_font=dict(color="#173D32", size=17, family="DM Sans, sans-serif"),
    legend=dict(font=dict(color="#24332C", size=13), title_font=dict(color="#173D32", size=13)),
    hoverlabel=dict(font=dict(color="#24332C", size=13), bgcolor="white", bordercolor="#DCE3D9"),
    margin=dict(l=40, r=40, t=60, b=40),
)


def style_plotly(fig, *, title=None, height=None, legend=None, **extra):
    layout = dict(_PLOTLY_BASE)
    if legend is not None:
        layout["legend"] = {**_PLOTLY_BASE.get("legend", {}), **legend}
    if title is not None:
        layout["title"] = title
    if height is not None:
        layout["height"] = height
    layout.update(extra)
    fig.update_layout(**layout)
    return fig


def metric(label, value, help_text: str | None = None):
    kwargs = {}
    if help_text:
        kwargs["help"] = help_text
    st.metric(label, value, **kwargs)


# ==========================================================
# 3. IDENTITE VISUELLE
# ==========================================================

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Playfair+Display:wght@600;700&display=swap');

:root {
    --forest: #173D32; --green: #28634E; --cream: #F7F6EF; --ink: #24332C; --muted: #4A5550;
    --text-color: #24332C;
    --background-color: #F7F6EF;
    --secondary-background-color: #FFFFFF;
    --primary-color: #28634E;
}

html, body, .stApp, .stMarkdown, .stText, [data-testid="stSidebar"] {
    font-family: 'DM Sans', sans-serif;
}
.stApp { background: var(--cream); color: var(--ink); font-size: 16.5px; line-height: 1.6; }
.stMarkdown, .stText { font-size: 16.5px; line-height: 1.6; }

/* ---------- SIDEBAR ---------- */
[data-testid="stSidebar"] { background: var(--forest); font-size: 15px; line-height: 1.55; }

[data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3,
[data-testid="stSidebar"] h4, [data-testid="stSidebar"] p, [data-testid="stSidebar"] span,
[data-testid="stSidebar"] label, [data-testid="stSidebar"] li, [data-testid="stSidebar"] small,
[data-testid="stSidebar"] strong, [data-testid="stSidebar"] [data-testid="stCaptionContainer"],
[data-testid="stSidebar"] [data-testid="stCaptionContainer"] * { color: #F7F6EF !important; }

[data-testid="stSidebar"] button,
[data-testid="stSidebar"] [data-testid^="stBaseButton"] {
    background-color: #F7F6EF !important;
    border: 1px solid #DCE3D9 !important;
    font-size: 14.5px !important;
}
[data-testid="stSidebar"] button:hover,
[data-testid="stSidebar"] [data-testid^="stBaseButton"]:hover {
    background-color: #EDEBDD !important;
    border-color: #6B947A !important;
}
[data-testid="stSidebar"] button:focus,
[data-testid="stSidebar"] [data-testid^="stBaseButton"]:focus {
    background-color: #FFFFFF !important;
    border-color: #6B947A !important;
}
[data-testid="stSidebar"] button, [data-testid="stSidebar"] button p,
[data-testid="stSidebar"] button span, [data-testid="stSidebar"] button div,
[data-testid="stSidebar"] [data-testid^="stBaseButton"],
[data-testid="stSidebar"] [data-testid^="stBaseButton"] p,
[data-testid="stSidebar"] [data-testid^="stBaseButton"] span,
[data-testid="stSidebar"] [data-testid^="stBaseButton"] div { color: var(--ink) !important; }

[data-testid="stSidebar"] .romarinova-logo {
    width: 100%;
    max-width: 380px;
    height: auto;
    display: block;
    margin: 0 0 16px 0;
    border-radius: 8px;
}

/* ---------- MAIN AREA ---------- */
[data-testid="stAppViewContainer"] :where(p, span, label, li, h1, h2, h3, h4, h5, h6, small, strong, em) {
    color: var(--ink);
}
[data-testid="stAppViewContainer"] :where([data-testid="stCaptionContainer"], [data-testid="stCaptionContainer"] *),
[data-testid="stAppViewContainer"] :where(.stCaption, .stCaption *) {
    color: var(--muted) !important;
    font-size: 13.5px !important;
    line-height: 1.5;
}

/* Markdown h4 used as in-page sub-headers (Atlas, Shortlist) */
.stMarkdown h4 {
    font-size: 18px;
    line-height: 1.4;
    letter-spacing: -0.1px;
}

[data-testid="stAppViewContainer"] [data-testid="stAlert"] {
    border-radius: 12px !important;
    font-size: 15px;
}
[data-testid="stAppViewContainer"] [data-testid="stAlert"] p,
[data-testid="stAppViewContainer"] [data-testid="stAlert"] * {
    color: var(--ink) !important;
    font-size: 15px !important;
}
[data-testid="stSidebar"] [data-testid="stAlert"] {
    background-color: #FFFFFF !important;
    border-left: 4px solid #C7A66A !important;
    border-radius: 12px !important;
}
[data-testid="stSidebar"] [data-testid="stAlert"] p,
[data-testid="stSidebar"] [data-testid="stAlert"] * {
    color: var(--ink) !important;
    font-size: 14px !important;
}

/* ---------- PILLS ---------- */
.health-pill {
    display: block;
    padding: 9px 13px;
    border-radius: 10px;
    font-size: 13px;
    font-weight: 600;
    letter-spacing: 0.2px;
    text-align: center;
    margin-bottom: 9px;
}
.health-pill.ok { background: #28634E; color: #FFFFFF !important; }
.health-pill.warn { background: #C7A66A; color: #24332C !important; }

/* ---------- LINKS / DOWNLOADS / PRIMARY ---------- */
[data-testid="stLinkButton"], [data-testid="stLinkButton"] *,
[data-testid^="stBaseLinkButton"], [data-testid^="stBaseLinkButton"] * {
    color: #FFFFFF !important;
    font-size: 14.5px !important;
}
[data-testid="stDownloadButton"], [data-testid="stDownloadButton"] *,
[data-testid^="stBaseDownloadButton"], [data-testid^="stBaseDownloadButton"] * {
    color: #FFFFFF !important;
    font-size: 14.5px !important;
}
[data-testid="stAppViewContainer"] [data-testid^="stBaseButton-primary"],
[data-testid="stAppViewContainer"] [data-testid^="stBaseButton-primary"] * {
    color: #FFFFFF !important;
    font-size: 14.5px !important;
}

/* ---------- CODE ---------- */
[data-testid="stCode"] {
    background-color: #173D32 !important;
    border-radius: 12px !important;
}
[data-testid="stCode"] pre,
[data-testid="stCode"] code,
[data-testid="stCode"] code span,
[data-testid="stCode"] pre code,
[data-testid="stCode"] pre code span,
[data-testid="stCode"] .token {
    color: #FFFFFF !important;
    background: transparent !important;
    -webkit-text-fill-color: #FFFFFF !important;
    font-size: 14px !important;
}

/* ---------- METRICS ---------- */
[data-testid="stMetric"] {
    background: #FFFFFF !important;
    border: 1px solid #E5E8DF !important;
    border-radius: 14px !important;
    padding: 18px 22px !important;
    min-height: auto !important;
}
[data-testid="stMetricLabel"] > div, [data-testid="stMetricLabel"] {
    color: var(--muted) !important;
    font-size: 14px !important;
    letter-spacing: 0.15px !important;
}
[data-testid="stMetricValue"] > div, [data-testid="stMetricValue"] {
    color: var(--forest) !important;
    font-size: 34px !important;
    font-weight: 700 !important;
    line-height: 1.2 !important;
}

/* ---------- PLOTLY ---------- */
.js-plotly-plot .plotly text,
.js-plotly-plot .plotly .xtick text,
.js-plotly-plot .plotly .ytick text,
.js-plotly-plot .plotly .legend text,
.js-plotly-plot .plotly .gtitle,
.js-plotly-plot .plotly .xtitle,
.js-plotly-plot .plotly .ytitle {
    fill: #24332C !important;
    font-size: 13px !important;
}

/* ---------- HEADINGS ---------- */
h1, h2, h3, h4 { color: var(--forest); }

/* ---------- HERO ---------- */
.hero {
    background: linear-gradient(120deg, var(--forest), var(--green));
    padding: 40px 42px;
    border-radius: 18px;
    color: white;
    margin-bottom: 26px;
}
.hero h1 {
    color: #FFFFFF !important;
    font-family: 'Playfair Display', serif;
    font-size: 52px;
    line-height: 1.12;
    letter-spacing: -0.8px;
    margin: 0 0 10px 0;
}
.hero p {
    color: #E6EFE8 !important;
    font-size: 18px;
    line-height: 1.55;
    max-width: 920px;
    margin: 0 0 14px 0;
}
.hero span { color: #F7F6EF !important; }
.hero .tag { color: #000000 !important; }

/* ---------- SECTION TITLES ---------- */
.section-title {
    font-family: 'Playfair Display', serif;
    color: var(--forest) !important;
    font-size: 34px;
    line-height: 1.2;
    letter-spacing: -0.5px;
    margin: 22px 0 0 0;
}
.sub-heading {
    font-family: 'Playfair Display', serif;
    color: var(--forest) !important;
    font-size: 24px;
    line-height: 1.28;
    letter-spacing: -0.3px;
    margin: 26px 0 8px 0;
}
.subtitle {
    color: var(--muted) !important;
    font-size: 15px;
    line-height: 1.55;
    margin-bottom: 18px;
}

/* ---------- CARDS ---------- */
.info-card {
    background: white;
    border: 1px solid #E5E8DF;
    border-left: 4px solid #6B947A;
    border-radius: 12px;
    padding: 18px 22px;
    margin-bottom: 14px;
}
.info-card h4 {
    color: var(--forest) !important;
    font-size: 17.5px;
    line-height: 1.35;
    letter-spacing: -0.15px;
    margin: 0 0 6px 0;
}
.info-card p {
    color: #4C5A52 !important;
    font-size: 15.5px;
    line-height: 1.55;
    margin-bottom: 0;
}

/* ---------- TAG / CHIP / FOOTER ---------- */
.tag {
    display: inline-block;
    background: #E8EFE8;
    color: #000000 !important;
    border-radius: 20px;
    padding: 5px 12px;
    font-size: 13px;
    font-weight: 500;
    letter-spacing: 0.2px;
    margin: 4px 4px 0 0;
}
.chip {
    display: block;
    text-align: center;
    background: #FFFFFF;
    border: 1px solid #E5E8DF;
    border-radius: 10px;
    padding: 11px 10px;
    font-size: 13.5px;
    font-weight: 600;
    color: #173D32;
    margin-bottom: 9px;
}
.footer {
    border-top: 1px solid #DCE3D9;
    margin-top: 44px;
    padding-top: 18px;
    color: var(--muted) !important;
    font-size: 13px;
    line-height: 1.6;
    text-align: center;
}

/* ---------- RESPONSIVE ---------- */
@media (max-width: 768px) {
    .hero { padding: 26px 24px; border-radius: 14px; }
    .hero h1 { font-size: 32px; letter-spacing: -0.5px; }
    .hero p { font-size: 16px; }
    .section-title { font-size: 24px; }
    .sub-heading { font-size: 19px; }
    .info-card h4 { font-size: 16px; }
    .info-card p { font-size: 14.5px; }
    [data-testid="stMetricValue"] > div,
    [data-testid="stMetricValue"] { font-size: 26px !important; }
    [data-testid="stMetricLabel"] > div,
    [data-testid="stMetricLabel"] { font-size: 13px !important; }
    [data-testid="stMetric"] { padding: 14px 16px !important; }
}
</style>
""", unsafe_allow_html=True)


_INFO_VARIANT_BORDERS = {
    "info": "#6B947A", "warning": "#C7A66A",
    "success": "#28634E", "note": "#A8B0AB",
}


def info_card(title, description, variant: str = "info"):
    border = _INFO_VARIANT_BORDERS.get(variant, _INFO_VARIANT_BORDERS["info"])
    st.markdown(
        f'<div class="info-card" style="border-left-color:{border};">'
        f'<h4>{esc(title)}</h4><p>{safe_md(description)}</p></div>',
        unsafe_allow_html=True,
    )


def section_title(title, subtitle=None):
    st.markdown(f'<h2 class="section-title">{esc(title)}</h2>', unsafe_allow_html=True)
    if subtitle:
        st.markdown(f'<p class="subtitle">{esc(subtitle)}</p>', unsafe_allow_html=True)


def sub_title(title, subtitle=None):
    st.markdown(f'<h3 class="sub-heading">{esc(title)}</h3>', unsafe_allow_html=True)
    if subtitle:
        st.markdown(f'<p class="subtitle">{esc(subtitle)}</p>', unsafe_allow_html=True)


def page_hero(title: str, subtitle: str, tags: list[str] | None = None):
    tags_html = "".join(f'<span class="tag">{esc(t)}</span>' for t in (tags or []))
    st.markdown(
        f'<div class="hero">'
        f'<h1>{esc(title)}</h1>'
        f'<p>{esc(subtitle)}</p>'
        f'{tags_html}'
        f'</div>',
        unsafe_allow_html=True,
    )


# ==========================================================
# 4. PAGE FUNCTIONS
# ==========================================================

def _page_dashboard():
    page_hero(
        "ROMARINOVA",
        "Computational exploration of Tunisian rosemary for type 2 diabetes research.",
        tags=["ETHNOPHARMACOLOGY", "BIOINFORMATICS", "NATURAL PRODUCTS"],
    )

    section_title("Research Dashboard", "Project overview.")

    pr_auc_display = fmt(_RF_PR_AUC) if _RF_PR_AUC is not None else "—"
    _model_label = _RF_MODEL or "n/a"
    _pr_auc_help = (
        f"Precision-Recall AUC on the held-out test split (model: {_model_label}). "
        f"Higher = better."
    )

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        metric("Screened Compounds", meta_val("screened_compound_n"),
               help_text="Unique compounds in the screened database.")
    with c2:
        metric("Shortlist (Pairs)", meta_val("N_star"),
               help_text="Compound–variety pairs in the final shortlist (N*).")
    with c3:
        metric("Unique Shortlisted", meta_val("U_star"),
               help_text="Unique compounds in the final shortlist (U*).")
    with c4:
        metric("ML Test PR-AUC", pr_auc_display, help_text=_pr_auc_help)

    st.write("")
    left, right = st.columns([1.2, 1])
    with left:
        sub_title("The research question")
        info_card(
            "Can Tunisian rosemary provide promising molecules?",
            "Does rosemary's chemical diversity include compounds that interact "
            "with carbohydrate-digesting enzymes?",
        )
        info_card(
            "Why α-glucosidase?",
            "A key enzyme in carbohydrate digestion; its inhibition is a validated "
            "strategy for controlling postprandial glucose.",
        )
    with right:
        sub_title("Scientific objectives")
        for item in [
            "Explore the chemical diversity of Tunisian rosemary.",
            "Organize compounds into a searchable library.",
            "Identify relevant biological targets.",
            "Explore molecule–target interactions computationally.",
            "Organize evidence for future validation.",
        ]:
            st.markdown(f"- {item}")

    st.markdown("---")
    section_title("Model Evaluation", "PR-AUC across candidate models and splits.")

    if not eval_table.empty and "pr_auc" in eval_table.columns:
        _eval_df = eval_table.dropna(subset=["pr_auc"]).copy()
        if not _eval_df.empty:
            _test = _eval_df[_eval_df["split"] == "test"]
            _order = (
                list(_test.sort_values("pr_auc", ascending=False)["model"])
                if not _test.empty
                else list(_eval_df["model"].drop_duplicates())
            )
            _eval_df["model"] = pd.Categorical(
                _eval_df["model"], categories=_order, ordered=True)
            fig_eval = px.bar(
                _eval_df, x="model", y="pr_auc", color="split", barmode="group",
                color_discrete_map={
                    "train": "#E5E8DF", "val": "#6B947A", "validation": "#6B947A",
                    "test": "#173D32",
                },
                labels={"pr_auc": "PR-AUC", "model": "", "split": "Split"},
            )
            fig_eval.update_xaxes(tickangle=-25)
            style_plotly(
                fig_eval, title="Model comparison — PR-AUC by split",
                yaxis=dict(title="PR-AUC", color="#24332C", gridcolor="#EDEFE9"),
                height=CHART_MEDIUM,
                legend=LEGEND_TOP_RIGHT,
            )
            st.plotly_chart(fig_eval, key="dash_model_cmp", **_CHART_WIDTH_KW)
            st.caption("Test PR-AUC is the honest metric; train is inflated by memorisation.")
    else:
        st.info("Evaluation table unavailable.")

    st.markdown("---")
    section_title("The ML Boundary")

    _n_i = _to_int_or_none(meta_val("screened_compound_n", None))
    _count_i = _to_int_or_none(meta_val("ml_rank_top_band_count", None))
    _val = meta_val("ml_rank_top_band_value", None)

    if (_n_i is not None and _n_i > 0 and _count_i is not None and _count_i > 1
            and _val is not None and _count_i <= _n_i):
        _pct = 100.0 * _count_i / _n_i
        _title = f"{_pct:.1f}% ML Tie-Band"
        _body = (f"{_count_i} of {_n_i} compounds share the same top ML_rank_score "
                 f"({fmt(_val)}). The model cannot distinguish between them.")
        _variant = "warning"
    elif _n_i is not None and _count_i == 1 and _val is not None:
        _title = "No ML Tie-Band Detected"
        _body = (f"The top ML_rank_score ({fmt(_val)}) belongs to a single molecule.")
        _variant = "success"
    elif _count_i is not None and _n_i is not None and _count_i > _n_i:
        _title = "ML Tie-Band (inconsistent counts)"
        _body = (f"`ml_rank_top_band_count` ({_count_i}) exceeds "
                 f"`screened_compound_n` ({_n_i}). Rerun the pipeline.")
        _variant = "warning"
    else:
        _title = "ML Tie-Band (metadata incomplete)"
        _body = "Tie-band count is unavailable in the pipeline metadata."
        _variant = "note"

    info_card(_title, _body, variant=_variant)

    sub_title("ML Score Distribution")
    _uniq_screened = (
        screened.dropna(subset=["compound_id"]).drop_duplicates("compound_id")
        if "compound_id" in screened.columns else screened
    )

    _ml_series = pd.Series(dtype=float)
    if "ML_rank_score" in _uniq_screened.columns and not _uniq_screened.empty:
        _ml_series = pd.to_numeric(
            _uniq_screened["ML_rank_score"], errors="coerce").dropna()

    if _ml_series.empty:
        st.info("ML score data unavailable.")
    else:
        _n_display = (
            int(_uniq_screened["compound_id"].nunique())
            if "compound_id" in _uniq_screened.columns else len(_ml_series)
        )
        _total = len(_ml_series)
        _n_unique = int(_ml_series.nunique())
        _tb_val: float | None
        try:
            _tb_val = float(_val) if _val is not None else None
        except (TypeError, ValueError):
            _tb_val = None

        if _total < _n_display:
            st.caption(f"{_n_display - _total} compound(s) excluded (no ML score).")

        if _n_unique <= 10:
            _counts = (
                _ml_series.round(6).value_counts().sort_index()
                .rename_axis("ML_rank_score").reset_index(name="count")
            )
            _counts["pct"] = 100.0 * _counts["count"] / _total
            _bar_colors = [
                "#C7A66A" if (_tb_val is not None and abs(v - _tb_val) < 1e-6)
                else "#28634E"
                for v in _counts["ML_rank_score"]
            ]
            fig = go.Figure(go.Bar(
                x=[f"{v:.4f}" for v in _counts["ML_rank_score"]],
                y=_counts["pct"],
                marker_color=_bar_colors,
                text=[f"{p:.1f}%" for p in _counts["pct"]],
                textposition="outside",
                customdata=_counts["count"],
                name="Compounds",
                hovertemplate=(
                    "ML_rank_score: %{x}<br>"
                    "Share: %{y:.1f}%<br>"
                    "Compounds: %{customdata}<extra></extra>"
                ),
            ))
            _y_max = float(_counts["pct"].max())
            fig.update_layout(showlegend=False, bargap=0.5)
            style_plotly(
                fig,
                title=f"ML_rank_score across {_n_display} compounds",
                xaxis=dict(title="ML Rank Score", color="#24332C", type="category"),
                yaxis=dict(title="% of Compounds", color="#24332C",
                           gridcolor="#EDEFE9", range=[0, _y_max * 1.25]),
                height=CHART_LARGE,
            )
            st.plotly_chart(fig, key="dash_ml_hist_cat", **_CHART_WIDTH_KW)
            if _tb_val is not None:
                _shared = _count_i if _count_i is not None else "?"
                st.caption(f"Gold = tie-band value **{fmt(_tb_val)}** "
                           f"shared by **{_shared}** compound(s).")
        else:
            fig = px.histogram(
                _uniq_screened, x="ML_rank_score", nbins=20,
                histnorm="percent",
                color_discrete_sequence=["#28634E"],
            )
            if _tb_val is not None:
                fig.add_vline(
                    x=_tb_val, line_dash="dash", line_color="#C7A66A", line_width=2,
                    annotation_text=(
                        f"Tie-band · {_count_i if _count_i is not None else '?'} compounds"
                    ),
                    annotation_position="top",
                )
            style_plotly(
                fig,
                title=f"ML_rank_score across {_n_display} compounds",
                xaxis=dict(title="ML Rank Score", color="#24332C", gridcolor="#EDEFE9"),
                yaxis=dict(title="% of Compounds", color="#24332C",
                           gridcolor="#EDEFE9"),
                height=CHART_LARGE,
            )
            st.plotly_chart(fig, key="dash_ml_hist", **_CHART_WIDTH_KW)

        if _count_i is not None and _count_i > 1 and _val is not None:
            st.caption(
                f"**{_count_i}** compounds share the top ML score of {fmt(_val)}."
            )

    sub_title("Score Component Distributions")
    _score_uniq = (
        comp_sl.dropna(subset=["compound_id"]).drop_duplicates("compound_id")
        if "compound_id" in comp_sl.columns else comp_sl
    )
    _component_specs = [
        ("ML_rank_score",     "ML rank",       "#173D32"),
        ("similarity_support", "Similarity",   "#28634E"),
        ("docking_score",     "Docking score", "#C7A66A"),
        ("final_score",       "Final score",   "#6B947A"),
    ]
    _present = [
        (col, lbl, colr) for col, lbl, colr in _component_specs
        if col in _score_uniq.columns and
        pd.to_numeric(_score_uniq[col], errors="coerce").dropna().shape[0] > 0
    ]
    if _present:
        _cols = st.columns(len(_present))
        for i, (col, lbl, colr) in enumerate(_present):
            with _cols[i]:
                _vals = pd.to_numeric(_score_uniq[col], errors="coerce").dropna()
                _fig_c = px.histogram(x=_vals, nbins=15, color_discrete_sequence=[colr])
                style_plotly(
                    _fig_c, title=lbl,
                    xaxis=dict(title="", color="#24332C", gridcolor="#EDEFE9"),
                    yaxis=dict(title="Count", color="#24332C", gridcolor="#EDEFE9"),
                    height=CHART_SMALL, showlegend=False,
                    margin=dict(l=20, r=20, t=40, b=30),
                )
                st.plotly_chart(_fig_c, key=f"dash_dist_{col}", **_CHART_WIDTH_KW)
        st.caption(
            "A flat ML rank confirms the tie-band; similarity and docking carry the ranking."
        )


def _page_atlas():
    page_hero(
        "Chemical Atlas",
        "Rosemary's phytochemical diversity across Tunisian varieties, with "
        "pipeline scores for each compound.",
        tags=["STRUCTURE", "DIVERSITY", "APPLICABILITY DOMAIN"],
    )

    if "variety" in screened.columns:
        _df_var = screened.dropna(subset=["variety", "compound_id"]).copy()
        if not _df_var.empty:
            sub_title("Botanical Variety Coverage",
                      "Unique compounds detected per variety.")

            _per_var = (
                _df_var.drop_duplicates(["compound_id", "variety"])
                .groupby("variety")["compound_id"].nunique()
                .reset_index(name="n_compounds")
                .sort_values("n_compounds", ascending=False)
            )
            if "known_active_flag" in _df_var.columns:
                _actives = (
                    _df_var[_df_var["known_active_flag"].fillna(False)]
                    .drop_duplicates(["compound_id", "variety"])
                    .groupby("variety")["compound_id"].nunique()
                    .reindex(_per_var["variety"]).fillna(0).astype(int).values
                )
            else:
                _actives = [0] * len(_per_var)

            _va, _vb = st.columns(2)
            with _va:
                _fig_v = px.bar(
                    _per_var, x="variety", y="n_compounds",
                    color_discrete_sequence=["#28634E"], text="n_compounds",
                )
                _fig_v.update_traces(textposition="outside")
                style_plotly(
                    _fig_v, title="Unique compounds per variety",
                    xaxis=dict(title="", color="#24332C"),
                    yaxis=dict(title="Compounds", color="#24332C", gridcolor="#EDEFE9"),
                    height=CHART_SMALL, showlegend=False,
                )
                st.plotly_chart(_fig_v, key="atlas_var_count", **_CHART_WIDTH_KW)
            with _vb:
                _share = pd.DataFrame({
                    "variety": _per_var["variety"],
                    "Known actives": _actives,
                    "Other": _per_var["n_compounds"].values - _actives,
                }).melt(id_vars="variety", var_name="class", value_name="count")
                _fig_v2 = px.bar(
                    _share, x="variety", y="count", color="class", barmode="stack",
                    color_discrete_map={"Known actives": "#C7A66A", "Other": "#DCE3D9"},
                )
                style_plotly(
                    _fig_v2, title="Known actives vs other",
                    xaxis=dict(title="", color="#24332C"),
                    yaxis=dict(title="Compounds", color="#24332C", gridcolor="#EDEFE9"),
                    height=CHART_SMALL,
                    legend=LEGEND_TOP_RIGHT,
                )
                st.plotly_chart(_fig_v2, key="atlas_var_active", **_CHART_WIDTH_KW)

    sub_title("Full Compound Database")

    if screened.empty:
        st.info("No compounds in the dataset.")
        return

    _clean = screened.dropna(subset=["compound_id"])
    unique_comps = _clean.drop_duplicates("compound_id").copy()

    if "variety" in _clean.columns:
        _vc = _clean.groupby("compound_id")["variety"].nunique().rename("variety_count")
        unique_comps = unique_comps.merge(
            _vc, left_on="compound_id", right_index=True, how="left")
        unique_comps["variety_count"] = (
            pd.to_numeric(unique_comps["variety_count"], errors="coerce")
            .fillna(0).astype(int)
        )

    if {"docking_missing", "docking_affinity_kcal_per_mol"}.issubset(_clean.columns):
        _ok = _clean[~_clean["docking_missing"].fillna(True).astype(bool)]
        _bd = (_ok.groupby("compound_id")["docking_affinity_kcal_per_mol"].min()
               .rename("best_docking"))
        unique_comps = unique_comps.merge(
            _bd, left_on="compound_id", right_index=True, how="left")

    if unique_comps.empty:
        st.info("No compounds to display.")
        return

    id_to_name: dict = {}
    for cid, name in zip(unique_comps["compound_id"], unique_comps["compound_name"]):
        if _is_missing(name) or not str(name).strip():
            id_to_name[cid] = f"Unnamed ({cid})"
        else:
            id_to_name[cid] = str(name).strip()

    compound_ids = list(id_to_name.keys())

    if st.session_state.get("atlas_compound") not in compound_ids:
        st.session_state["atlas_compound"] = compound_ids[0]

    selected_id = st.selectbox(
        "🔍 Select a molecule",
        options=compound_ids,
        format_func=lambda cid: id_to_name.get(cid, f"Unnamed ({cid})"),
        key="atlas_compound",
    )

    all_rows = screened[screened["compound_id"] == selected_id]
    if all_rows.empty:
        st.warning("Selected compound has no rows in the screened table.")
        return
    comp_data = all_rows.iloc[0]

    dq_issues = []
    for _col, _label in (("canonical_smiles", "SMILES"),
                         ("inchikey", "InChIKey")):
        if _col in all_rows.columns:
            _stripped = all_rows[_col].dropna().astype(str).str.strip()
            _variants = _stripped[_stripped != ""].unique()
            if len(_variants) > 1:
                dq_issues.append(f"{_label} ({len(_variants)} values)")
    if dq_issues:
        st.warning(f"Varieties disagree on {'; '.join(dq_issues)} — showing first row.")

    col1, col2 = st.columns([1, 2])
    with col1:
        st.markdown("#### 2D Structure")
        smiles = row_get(comp_data, "canonical_smiles", None)
        img = mol_image(smiles)
        if img is not None:
            st.image(img, caption=f"{id_to_name.get(selected_id)} — 2D depiction.")
        else:
            st.warning("Could not parse SMILES.")

    with col2:
        st.markdown("#### Molecule Details")
        st.markdown("**Canonical SMILES**")
        _smiles_val = row_get(comp_data, "canonical_smiles", None)
        st.code(_smiles_val if isinstance(_smiles_val, str) and _smiles_val else "—",
                language=None)

        st.markdown("**InChIKey**")
        inchikey = row_get(comp_data, "inchikey", None)
        st.code(inchikey if isinstance(inchikey, str) and inchikey else "—",
                language=None)

        c1, c2, c3 = st.columns(3)
        with c1:
            metric("ML Rank", fmt(row_get(comp_data, "ML_rank_score", None)),
                   help_text="Percentile rank from the Random Forest model.")
        with c2:
            metric("Similarity", fmt(row_get(comp_data, "similarity_support", None)),
                   help_text="Max Tanimoto similarity to ChEMBL training actives.")
        with c3:
            best_dock, best_var = best_docking_row(all_rows)
            if best_dock is None:
                metric("Best Docking", "No successful dock")
            else:
                metric("Best Docking", f"{best_dock:.2f} kcal/mol",
                       help_text="Best raw Vina affinity — not a binding free energy.")
                if (best_var and "variety" in all_rows.columns
                        and all_rows["variety"].nunique() > 1):
                    st.caption(f"From: {best_var}")

        if isinstance(inchikey, str) and inchikey and inchikey != "—":
            st.link_button(
                "Explore in PubChem",
                f"https://pubchem.ncbi.nlm.nih.gov/#query={quote(inchikey)}",
            )

    st.markdown("---")

    sub_title("Structure-Activity Landscape")
    if "similarity_support" in unique_comps.columns and "best_docking" in unique_comps.columns:
        scatter_df = unique_comps.dropna(subset=["similarity_support", "best_docking"])
        if not scatter_df.empty:
            _plot = scatter_df.copy()
            if "ad_flag" in _plot.columns:
                _plot["ad_flag"] = _plot["ad_flag"].fillna("unknown").astype(str)
            else:
                _plot["ad_flag"] = "unknown"
            fig = px.scatter(
                _plot, x="similarity_support", y="best_docking",
                color="ad_flag", symbol="ad_flag",
                hover_name="compound_name",
                labels={
                    "similarity_support": "Similarity Support",
                    "best_docking": "Best Docking (kcal/mol)",
                    "ad_flag": "Applicability domain",
                },
                color_discrete_map={
                    AD_WITHIN: "#28634E",
                    AD_EXTRAPOLATION: "#C7A66A",
                    "unknown": "#A8B0AB",
                },
                symbol_map={
                    AD_WITHIN: "circle",
                    AD_EXTRAPOLATION: "diamond",
                    "unknown": "x",
                },
            )
            fig.update_yaxes(autorange="reversed")
            fig.update_layout(legend={**LEGEND_TOP_RIGHT, "title": "Applicability domain"})
            style_plotly(
                fig,
                xaxis=dict(title="Similarity Support", color="#24332C", gridcolor="#EDEFE9"),
                yaxis=dict(title="Best Docking (kcal/mol) — lower is better",
                           color="#24332C", gridcolor="#EDEFE9"),
                height=CHART_MEDIUM,
            )
            st.plotly_chart(fig, key="atlas_scatter", **_CHART_WIDTH_KW)
            st.caption(
                "Upper-right = high similarity + favorable docking (y-axis inverted)."
            )
            if len(scatter_df) < len(unique_comps):
                st.caption(f"{len(scatter_df)} of {len(unique_comps)} compounds plotted.")
        else:
            st.info("No compounds with both similarity and docking data.")
    else:
        st.info("Similarity or docking data unavailable.")

    sub_title("Similarity Distribution by Applicability Domain")
    if "similarity_support" in unique_comps.columns:
        _sim_ad = unique_comps.dropna(subset=["similarity_support"]).copy()
        if "ad_flag" in _sim_ad.columns:
            _sim_ad["ad_flag"] = _sim_ad["ad_flag"].fillna("unknown").astype(str)
        else:
            _sim_ad["ad_flag"] = "unknown"
        if not _sim_ad.empty:
            _fig_sim = px.histogram(
                _sim_ad, x="similarity_support", color="ad_flag",
                nbins=20, barmode="overlay", histnorm="percent",
                color_discrete_map={
                    AD_WITHIN: "#28634E",
                    AD_EXTRAPOLATION: "#C7A66A",
                    "unknown": "#A8B0AB",
                },
                labels={"similarity_support": "Similarity Support",
                        "ad_flag": "Applicability domain"},
            )
            _fig_sim.update_traces(opacity=0.75)
            style_plotly(
                _fig_sim, title="Where rosemary compounds sit in similarity space",
                xaxis=dict(title="Similarity Support", color="#24332C", gridcolor="#EDEFE9"),
                yaxis=dict(title="% of Compounds", color="#24332C",
                           gridcolor="#EDEFE9"),
                height=CHART_MEDIUM,
                legend=LEGEND_TOP_RIGHT,
            )
            st.plotly_chart(_fig_sim, key="atlas_sim_ad", **_CHART_WIDTH_KW)

            _n_extrap = int((_sim_ad["ad_flag"] == AD_EXTRAPOLATION).sum())
            _n_total = len(_sim_ad)
            if _n_total and _n_extrap >= _n_total:
                st.caption(f"All {_n_total} compounds sit outside the training domain.")
            else:
                st.caption(f"{_n_extrap} of {_n_total} compounds outside the training domain.")

    sub_title("Lipinski / Veber Chemical Space")
    if "canonical_smiles" in unique_comps.columns:
        _smiles_series = unique_comps["canonical_smiles"].fillna("").astype(str)
        _smiles_tuple = tuple(_smiles_series.tolist())
        _desc_rows = _compute_descriptors(_smiles_tuple)
        _desc_df = pd.DataFrame(_desc_rows, index=unique_comps.index)
        _merged = unique_comps.join(_desc_df).dropna(subset=["MW", "LogP"])
        if not _merged.empty:
            _color_arg = (
                "developability_score" if "developability_score" in _merged.columns
                else None
            )
            _fig_lip = px.scatter(
                _merged, x="MW", y="LogP",
                color=_color_arg,
                color_continuous_scale=["#C7A66A", "#6B947A", "#173D32"],
                hover_name="compound_name",
                labels={
                    "MW": "Molecular weight (Da)",
                    "LogP": "cLogP",
                    "developability_score": "Developability",
                },
            )
            _fig_lip.add_vline(
                x=500, line_dash="dash", line_color="#A8B0AB", line_width=1.5,
                annotation_text="RO5 MW ≤ 500", annotation_position="top",
            )
            _fig_lip.add_hline(
                y=5, line_dash="dash", line_color="#A8B0AB", line_width=1.5,
                annotation_text="RO5 LogP ≤ 5", annotation_position="right",
            )
            style_plotly(
                _fig_lip, title="Molecular weight vs cLogP (RO5 reference)",
                xaxis=dict(title="Molecular weight (Da)", color="#24332C",
                           gridcolor="#EDEFE9"),
                yaxis=dict(title="cLogP", color="#24332C", gridcolor="#EDEFE9"),
                height=CHART_LARGE,
            )
            st.plotly_chart(_fig_lip, key="atlas_lipinski", **_CHART_WIDTH_KW)

            _n_above = int(((_merged["MW"] > 500) | (_merged["LogP"] > 5)).sum())
            st.caption(
                f"{_n_above} of {len(_merged)} compounds exceed at least one RO5 cutoff "
                f"(typically large diterpenes)."
            )
        else:
            st.info("No valid SMILES for descriptor calculation.")

    sub_title("Compound Table")
    filter_cols = st.columns(3)
    with filter_cols[0]:
        _ad_opts = (
            sorted(unique_comps["ad_flag"].dropna().astype(str).unique())
            if "ad_flag" in unique_comps.columns else []
        )
        ad_filter = st.multiselect(
            "Applicability domain", options=_ad_opts, default=_ad_opts,
            key="atlas_ad_filter",
        )
    with filter_cols[1]:
        known_only = st.checkbox("Known actives only", value=False, key="atlas_known_only")
    with filter_cols[2]:
        dev_range: tuple[float, float] | None = None
        if "developability_score" in unique_comps.columns:
            _dv = pd.to_numeric(
                unique_comps["developability_score"], errors="coerce").dropna()
            if not _dv.empty:
                _dv_min, _dv_max = float(_dv.min()), float(_dv.max())
                if _dv_min < _dv_max:
                    dev_range = st.slider(
                        "Developability score",
                        min_value=_dv_min, max_value=_dv_max,
                        value=(_dv_min, _dv_max),
                        key="atlas_dev_range",
                    )
                else:
                    st.caption(f"Developability score: {_dv_min:.2f} (uniform)")
                    dev_range = (_dv_min, _dv_max)

    _filtered = unique_comps.copy()
    if ad_filter and "ad_flag" in _filtered.columns:
        _filtered = _filtered[_filtered["ad_flag"].astype(str).str.strip().isin(ad_filter)]
    if known_only and "known_active_flag" in _filtered.columns:
        _filtered = _filtered[_filtered["known_active_flag"].fillna(False).astype(bool)]
    if dev_range is not None and "developability_score" in _filtered.columns:
        _filtered = _filtered[
            (pd.to_numeric(_filtered["developability_score"], errors="coerce") >= dev_range[0]) &
            (pd.to_numeric(_filtered["developability_score"], errors="coerce") <= dev_range[1])
        ]

    st.caption(f"Showing {len(_filtered)} of {len(unique_comps)} compounds.")

    display_cols = select_existing(
        _filtered,
        ["compound_id", "compound_name", "variety_count", "ad_flag", "structure_status",
         "known_active_flag", "best_docking", "developability_score"],
    )
    col_cfg = {
        "variety_count": st.column_config.NumberColumn(
            label="Varieties",
            help="Number of botanical varieties where this compound was detected."),
        "ad_flag": st.column_config.TextColumn(
            label="AD",
            help="`within_AD` = inside training domain; `extrapolation` = outside."),
        "known_active_flag": st.column_config.CheckboxColumn(
            label="Known active", help="Matches a ChEMBL active InChIKey."),
        "best_docking": st.column_config.NumberColumn(
            label="Best docking", format="%.2f",
            help="Best raw Vina estimate — not a binding free energy."),
        "developability_score": st.column_config.ProgressColumn(
            label="Developability", format="%.2f", min_value=0.0, max_value=1.0,
            help="Composite score from Lipinski, Veber, PAINS, and aggregator filters."),
    }

    _sort_cols: list[str] = []
    _sort_ascending: list[bool] = []
    if "best_docking" in _filtered.columns:
        _sort_cols.append("best_docking")
        _sort_ascending.append(True)
    if "developability_score" in _filtered.columns:
        _sort_cols.append("developability_score")
        _sort_ascending.append(False)
    if _sort_cols:
        _sorted = _filtered.sort_values(
            _sort_cols, ascending=_sort_ascending, na_position="last",
        ).reset_index(drop=True)
    else:
        _sorted = _filtered.reset_index(drop=True)

    if _sorted.empty:
        st.info("No compounds match the current filters.")
    else:
        st.dataframe(
            _sorted[display_cols], hide_index=True,
            column_config=col_cfg, key="atlas_table",
            **_DATAFRAME_WIDTH_KW,
        )
        st.download_button(
            "⬇️ Download filtered table (CSV)",
            data=_sorted[display_cols].to_csv(index=False).encode("utf-8"),
            file_name="romarinova_atlas_filtered.csv",
            mime="text/csv",
            key="atlas_download",
            help="CSV of the current filtered and sorted table.",
        )


def _page_shortlist():
    page_hero(
        "α-Glucosidase Screening",
        "The final consensus shortlist generated by Romarinova AI.",
        tags=["ALPHA-GLUCOSIDASE", "CONSENSUS SCORING", "SENSITIVITY"],
    )

    coverage_gap = meta_val("coverage_gap")
    coverage_gap_display = "—" if coverage_gap == "None" else coverage_gap
    sensitivity = str(meta_val("sensitivity_label", "unknown")).upper()
    _n = meta_val("N_star")
    _u = meta_val("U_star")

    st.info(
        f"**N∗ = {_n}** pairs  |  **U∗ = {_u}** unique compounds  |  "
        f"Coverage gap: {coverage_gap_display}  |  Sensitivity: {sensitivity}"
    )

    df = comp_sl.copy()

    if df.empty:
        st.warning("No compounds in the shortlist.")
        return

    if {"final_rank", "ligand_only_rank"}.issubset(df.columns):
        _lo = pd.to_numeric(df["ligand_only_rank"], errors="coerce")
        _fi = pd.to_numeric(df["final_rank"], errors="coerce")
        df["rank_delta"] = (_lo - _fi)

    if "final_rank" in df.columns:
        df = df.sort_values("final_rank", na_position="last").reset_index(drop=True)

    if "compound_id" in df.columns:
        _before = len(df)
        df = df.drop_duplicates("compound_id", keep="first").reset_index(drop=True)
        _n_dup = _before - len(df)
        if _n_dup > 0:
            st.caption(
                f"Collapsed {_n_dup} duplicate row(s); kept best-ranked per compound."
            )

    sub_title("Score Contribution")

    plot_df = (
        comp_sl.dropna(subset=["compound_id"]).drop_duplicates("compound_id").copy()
        if "compound_id" in comp_sl.columns else comp_sl.copy()
    )
    if "final_rank" in plot_df.columns:
        plot_df = plot_df.sort_values("final_rank", na_position="last").reset_index(drop=True)

    _weights = meta.get("scoring_weights")
    if not isinstance(_weights, dict):
        _weights = _DEFAULT_WEIGHTS
    try:
        w_ml = float(_weights.get("ML", _DEFAULT_WEIGHTS["ML"]))
        w_sim = float(_weights.get("similarity", _DEFAULT_WEIGHTS["similarity"]))
        w_dock = float(_weights.get("docking", _DEFAULT_WEIGHTS["docking"]))
    except (TypeError, ValueError):
        w_ml, w_sim, w_dock = (_DEFAULT_WEIGHTS["ML"], _DEFAULT_WEIGHTS["similarity"],
                               _DEFAULT_WEIGHTS["docking"])

    _required = ["final_score", "ML_rank_score", "similarity_support", "docking_score"]
    if len(plot_df) > 0 and all(c in plot_df.columns for c in _required):
        _n_before = len(plot_df)
        plot_df = plot_df.dropna(subset=_required).reset_index(drop=True)
        _n_dropped = _n_before - len(plot_df)
        if _n_dropped > 0:
            st.caption(f"{_n_dropped} compound(s) excluded (incomplete scores).")

        if plot_df.empty:
            st.info("No complete rows for the score-contribution chart.")
        else:
            plot_df["ML_Contribution"] = w_ml * plot_df["ML_rank_score"]
            plot_df["Sim_Contribution"] = w_sim * plot_df["similarity_support"]
            plot_df["Dock_Contribution"] = w_dock * plot_df["docking_score"]
            plot_df["Base_Contribution"] = plot_df["final_score"] - (
                plot_df["ML_Contribution"]
                + plot_df["Sim_Contribution"]
                + plot_df["Dock_Contribution"]
            )

            _bad_base = plot_df["Base_Contribution"] < -0.001
            if _bad_base.any():
                st.caption(
                    f"{int(_bad_base.sum())} compound(s) show negative base contribution "
                    "— check that scoring weights match the pipeline."
                )
                plot_df = plot_df[~_bad_base].reset_index(drop=True)

            if plot_df.empty:
                st.warning(
                    "All compounds excluded by negative-base-contribution filter — "
                    "check scoring weights in `shortlist_metadata.json`."
                )
            else:
                ctl_a, ctl_b = st.columns([2, 1])
                with ctl_a:
                    if len(plot_df) > 3:
                        _top_n = st.slider(
                            "Compounds displayed",
                            min_value=3, max_value=len(plot_df),
                            value=min(15, len(plot_df)), key="score_top_n",
                        )
                    else:
                        _top_n = len(plot_df)
                        st.caption(f"Showing all {_top_n} compound(s).")
                with ctl_b:
                    _orient = st.radio(
                        "Bar orientation", ["Vertical", "Horizontal"],
                        horizontal=True, key="score_orient",
                    )

                _disp = plot_df.sort_values("final_score", ascending=False).head(_top_n)

                def _build_traces(orient: str):
                    if orient == "Horizontal":
                        fig_local = go.Figure()
                        fig_local.add_trace(go.Bar(
                            y=_disp["compound_name"], x=_disp["ML_Contribution"],
                            name=f"ML ({w_ml*100:.0f}%)", marker_color="#173D32",
                            orientation="h"))
                        fig_local.add_trace(go.Bar(
                            y=_disp["compound_name"], x=_disp["Sim_Contribution"],
                            name=f"Similarity ({w_sim*100:.0f}%)", marker_color="#6B947A",
                            orientation="h"))
                        fig_local.add_trace(go.Bar(
                            y=_disp["compound_name"], x=_disp["Dock_Contribution"],
                            name=f"Docking ({w_dock*100:.0f}%)", marker_color="#C7A66A",
                            orientation="h"))
                        fig_local.add_trace(go.Bar(
                            y=_disp["compound_name"], x=_disp["Base_Contribution"],
                            name="Developability & Novelty", marker_color="#E5E8DF",
                            orientation="h"))
                        fig_local.update_layout(
                            barmode="stack",
                            height=max(CHART_SMALL, 40 * len(_disp) + 120),
                        )
                        style_plotly(
                            fig_local, title="Score Contribution by Component",
                            xaxis=dict(title="Contribution to final_score",
                                       color="#24332C", gridcolor="#EDEFE9"),
                            yaxis=dict(title="", color="#24332C", automargin=True),
                            legend=LEGEND_TOP_RIGHT,
                        )
                        return fig_local
                    fig_local = go.Figure()
                    fig_local.add_trace(go.Bar(
                        x=_disp["compound_name"], y=_disp["ML_Contribution"],
                        name=f"ML ({w_ml*100:.0f}%)", marker_color="#173D32"))
                    fig_local.add_trace(go.Bar(
                        x=_disp["compound_name"], y=_disp["Sim_Contribution"],
                        name=f"Similarity ({w_sim*100:.0f}%)", marker_color="#6B947A"))
                    fig_local.add_trace(go.Bar(
                        x=_disp["compound_name"], y=_disp["Dock_Contribution"],
                        name=f"Docking ({w_dock*100:.0f}%)", marker_color="#C7A66A"))
                    fig_local.add_trace(go.Bar(
                        x=_disp["compound_name"], y=_disp["Base_Contribution"],
                        name="Developability & Novelty", marker_color="#E5E8DF"))
                    fig_local.update_layout(barmode="stack")
                    fig_local.update_xaxes(tickangle=-35)
                    style_plotly(
                        fig_local, title="Score Contribution by Component",
                        yaxis=dict(title="Contribution to final_score",
                                   color="#24332C", gridcolor="#EDEFE9"),
                        height=CHART_MEDIUM,
                        margin=dict(l=40, r=40, t=60, b=90),
                        legend=LEGEND_TOP_RIGHT,
                    )
                    return fig_local

                fig = _build_traces(_orient)
                st.plotly_chart(fig, key="shortlist_contrib", **_CHART_WIDTH_KW)
                st.caption(
                    "Developability & Novelty = final_score − ML − Similarity − Docking."
                )

                _ml_unique_chart = int(
                    pd.to_numeric(plot_df["ML_rank_score"], errors="coerce")
                    .dropna().nunique()
                )
                if _ml_unique_chart == 1:
                    st.caption(
                        "Flat ML block = tie-band; similarity and docking drive the ranking."
                    )
                else:
                    st.caption(f"{_ml_unique_chart} distinct ML scores contribute to ordering.")

    sub_title("Multi-Criteria Ranking View")
    _pc_df = (
        comp_sl.dropna(subset=["compound_id"]).drop_duplicates("compound_id").copy()
        if "compound_id" in comp_sl.columns else comp_sl.copy()
    )
    _pc_needed = ["ML_rank_score", "similarity_support", "docking_score", "final_score"]
    if "developability_score" in _pc_df.columns:
        _pc_needed.append("developability_score")
    if all(c in _pc_df.columns for c in _pc_needed) and not _pc_df.empty:
        _pc_df = _pc_df.dropna(subset=_pc_needed).copy()
        if "final_rank" in _pc_df.columns:
            _pc_df = _pc_df.sort_values("final_rank").head(20)

        if not _pc_df.empty:
            _pc_dims = [
                dict(label="ML rank",    values=_pc_df["ML_rank_score"]),
                dict(label="Similarity", values=_pc_df["similarity_support"]),
                dict(label="Docking",    values=_pc_df["docking_score"]),
            ]
            if "developability_score" in _pc_df.columns:
                _pc_dims.append(dict(label="Developability",
                                     values=_pc_df["developability_score"]))
            _pc_dims.append(dict(label="Final score", values=_pc_df["final_score"]))

            _pc_color_vals = (
                _pc_df["final_rank"] if "final_rank" in _pc_df.columns
                else _pc_df["final_score"]
            )
            fig_pc = go.Figure(go.Parcoords(
                line=dict(
                    color=_pc_color_vals,
                    colorscale=[[0, "#173D32"], [1, "#C7A66A"]],
                    showscale=True,
                    colorbar=dict(title="Final rank", len=0.6),
                ),
                dimensions=_pc_dims,
            ))
            style_plotly(
                fig_pc, title=f"Parallel coordinates — top {len(_pc_df)} compounds",
                height=CHART_LARGE,
                margin=dict(l=60, r=80, t=60, b=40),
            )
            st.plotly_chart(fig_pc, key="shortlist_parcoords", **_CHART_WIDTH_KW)
            st.caption(
                "Lines = compounds; dark green = best rank. Widest axes differentiate most."
            )

    col1, col2 = st.columns(2)

    locked_cols = select_existing(
        df, ["final_rank", "compound_name", "final_score",
             "ML_rank_score", "similarity_support", "docking_score", "rank_delta"])
    indep_cols = select_existing(
        df, ["ligand_only_rank", "compound_name", "ligand_only_score",
             "ML_rank_score", "similarity_support", "docking_score", "rank_delta"])

    col_cfg = {
        "final_score": st.column_config.NumberColumn(
            format="%.3f", help="Weighted consensus score."),
        "ligand_only_score": st.column_config.NumberColumn(
            format="%.3f", help="Weighted score with ML weight = 0.00."),
        "ML_rank_score": st.column_config.NumberColumn(
            format="%.3f", help="Validation-rank percentile. Not a probability."),
        "docking_score": st.column_config.NumberColumn(
            format="%.3f",
            help="Transformed Vina affinity — not comparable to raw Atlas values."),
        "rank_delta": st.column_config.NumberColumn(
            label="Δ rank", format="%+d",
            help="ligand_only_rank − final_rank. Positive = ML helped this compound."),
    }

    with col1:
        st.markdown("#### Locked Ranking (final_score)")
        st.dataframe(
            df.sort_values("final_rank", na_position="last")[locked_cols],
            hide_index=True, column_config=col_cfg, key="df_locked",
            **_DATAFRAME_WIDTH_KW,
        )

    with col2:
        st.markdown("#### ML-Independent Ranking (ligand_only_score)")
        st.dataframe(
            df.sort_values("ligand_only_rank", na_position="last")[indep_cols],
            hide_index=True, column_config=col_cfg, key="df_indep",
            **_DATAFRAME_WIDTH_KW,
        )

    st.caption("Δ rank: positive = ML helped; negative = ML hurt.")

    with st.expander("Consensus Analysis", expanded=False):
        st.caption(
            "Compares locked `final_score` ranking against ML-independent "
            "`ligand_only_score` ranking."
        )

        _ml_vals = pd.to_numeric(comp_sl["ML_rank_score"], errors="coerce").dropna()
        _ml_unique = int(_ml_vals.nunique())
        overlap = meta_val("pair_overlap", "unknown")

        if _ml_vals.empty:
            _card_title = "ML Score Data Unavailable"
            _card_body = (f"Sensitivity: **{sensitivity}**, pair overlap **{overlap}**. "
                          f"No ML_rank_score values in compound shortlist.")
            _card_variant = "note"
        elif _ml_unique == 1 and len(_ml_vals) > 1:
            _ml_val_str = fmt(_ml_vals.iloc[0])
            _card_title = "The ML Tie-Band"
            _card_body = (f"Sensitivity **{sensitivity}**, pair overlap **{overlap}**. "
                          f"All shortlisted compounds share ML_rank_score "
                          f"({_ml_val_str}); ranking is driven by **similarity** "
                          f"and **docking**.")
            _card_variant = "warning"
        elif _ml_unique == 1 and len(_ml_vals) == 1:
            _card_title = "Single Compound — No Tie Analysis"
            _card_body = (f"Shortlist contains one compound (ML_rank_score "
                          f"{fmt(_ml_vals.iloc[0])}). Tie analysis requires ≥2.")
            _card_variant = "note"
        else:
            _card_title = "ML Rank Contribution"
            _card_body = (f"Sensitivity **{sensitivity}**, pair overlap **{overlap}**. "
                          f"{_ml_unique} distinct ML scores — removing ML may change "
                          f"the ordering. See Δ rank column.")
            _card_variant = "info"

        info_card(_card_title, _card_body, variant=_card_variant)

    st.download_button(
        "⬇️ Download compound-level shortlist (CSV)",
        data=df.to_csv(index=False).encode("utf-8"),
        file_name="romarinova_shortlist_compounds.csv",
        mime="text/csv",
        key="download_shortlist_compounds",
        help="CSV of the compound-level shortlist, sorted by final_rank.",
    )


def _page_targets():
    page_hero(
        "Biological Target Explorer",
        "α-glucosidase: biological relevance and 3D structure.",
        tags=["STRUCTURAL BIOLOGY", "DOCKING", "INTERPRETATION"],
    )

    pdb_id = dock_cfg.get("receptor_pdb_id", "3L4W")
    col1, col2 = st.columns([1, 1.5])

    with col1:
        info_card(
            "α-Glucosidase (MGAM)",
            "Key enzyme in carbohydrate digestion; a validated target for "
            "postprandial glucose control.",
        )
        info_card(
            "Target Details",
            f"Docking target: human maltase-glucoamylase C-terminal domain "
            f"(PDB: {pdb_id}).",
        )
        st.link_button(
            "PubMed — α-glucosidase inhibitors",
            "https://pubmed.ncbi.nlm.nih.gov/?term=alpha-glucosidase+inhibitors",
        )
        st.caption("In silico only — biochemical inhibition assays required.")

    with col2:
        sub_title("3D Receptor Structure")
        pdb_rel = dock_cfg.get("receptor_pdb_path")
        pdb_path = ((PROJECT_ROOT / pdb_rel).resolve() if pdb_rel
                    else PROJECT_ROOT / "data" / f"{pdb_id}.pdb")

        if not pdb_path.exists():
            st.info(f"PDB file `{pdb_path.name}` not found. Run "
                    f"`scripts/prepare_receptor.py` to download the receptor.")
            st.link_button(
                f"Fetch PDB {pdb_id} from RCSB",
                f"https://www.rcsb.org/structure/{pdb_id}",
            )
        else:
            cx, cy, cz = parse_grid_center(dock_cfg.get("grid_center", [0, 0, 0]))
            ctl_a, ctl_b = st.columns(2)
            with ctl_a:
                show_surface = st.toggle("Molecular surface", value=False,
                                         key="pdb_surface",
                                         help="Translucent van der Waals surface.")
            with ctl_b:
                show_site = st.toggle("Docking site", value=True, key="pdb_site",
                                      help="Highlight the docking grid centre.")
            try:
                pdb_mtime = pdb_path.stat().st_mtime
            except OSError:
                pdb_mtime = 0.0
            try:
                html_str = render_pdb_html(
                    str(pdb_path), cx, cy, cz, 600, 450,
                    show_surface, show_site, pdb_mtime,
                )
                if html_str is None:
                    st.info("PDB file could not be read.")
                else:
                    if hasattr(st, "iframe"):
                        st.iframe(html_str, height=450)
                    else:
                        components.html(html_str, height=450)
                    st.caption("Gold sphere = docking grid centre (catalytic site).")
            except Exception as e:
                log.exception("render_pdb_html failed")
                st.error(f"Could not render 3D structure: `{e}`")


def _page_pipeline():
    page_hero(
        "Bioinformatics & Docking",
        "Workflow from ML screening to consensus shortlist.",
        tags=["WORKFLOW", "REPRODUCIBILITY", "DOCKING"],
    )

    pr_auc = fmt(_RF_PR_AUC) if _RF_PR_AUC is not None else "—"
    pdb_id = dock_cfg.get("receptor_pdb_id", "3L4W")
    _n_act = _to_int_or_none(meta.get("n_actives_train"))
    _n_screened_n = _to_int_or_none(meta_val("screened_compound_n", None))
    _n_star = meta_val("N_star")
    _u_star = meta_val("U_star")

    section_title("Pipeline at a Glance", "Parameters defining this run.")

    m1, m2, m3, m4 = st.columns(4)
    with m1:
        metric("Classifier", _RF_MODEL or "Random Forest",
               help_text="Selected on validation PR-AUC.")
    with m2:
        metric("Test PR-AUC", pr_auc,
               help_text="Precision-recall AUC on the held-out test split.")
    with m3:
        metric("Training Actives",
               _n_act if _n_act is not None else "—",
               help_text="Active compounds in the ChEMBL MGAM training set.")
    with m4:
        metric("Receptor", f"PDB {pdb_id}",
               help_text="Human maltase-glucoamylase (MGAM) C-terminal domain.")

    st.write("")

    section_title("Consensus Score Composition",
                  "How components combine into `final_score`.")

    _weights = meta.get("scoring_weights")
    _weights_from_meta = isinstance(_weights, dict)
    if not _weights_from_meta:
        _weights = _DEFAULT_WEIGHTS

    _weight_rows = [
        ("ML ranking", _weights.get("ML", _DEFAULT_WEIGHTS["ML"])),
        ("Similarity", _weights.get("similarity", _DEFAULT_WEIGHTS["similarity"])),
        ("Developability",
         _weights.get("developability", _DEFAULT_WEIGHTS["developability"])),
        ("Novelty", _weights.get("novelty", _DEFAULT_WEIGHTS["novelty"])),
        ("Docking", _weights.get("docking", _DEFAULT_WEIGHTS["docking"])),
    ]
    _weight_df = pd.DataFrame(_weight_rows, columns=["component", "weight"])

    _w_max = float(_weight_df["weight"].max()) if not _weight_df.empty else 0.4
    fig_w = px.bar(
        _weight_df, x="weight", y="component", orientation="h",
        text=[f"{float(w)*100:.0f}%" for w in _weight_df["weight"]],
        color_discrete_sequence=["#28634E"],
    )
    fig_w.update_traces(textposition="outside")
    style_plotly(
        fig_w, title="Component weights in final_score",
        xaxis=dict(title="Weight", color="#24332C", gridcolor="#EDEFE9",
                   range=[0, max(0.4, _w_max * 1.25)]),
        yaxis=dict(title="", color="#24332C", automargin=True),
        height=CHART_SMALL, showlegend=False,
    )
    st.plotly_chart(fig_w, key="pipeline_weights", **_CHART_WIDTH_KW)

    _formula = " + ".join(
        f"{float(w)*100:.0f}% · {n}" for n, w in _weight_rows
    )
    st.markdown(f"**final_score = {_formula}**")

    if "ad_flag" in screened.columns:
        _uniq_screened = (
            screened.dropna(subset=["compound_id"])
            .drop_duplicates("compound_id")
        )
        _n_extrap = _count_extrapolated_unique_compounds()
        _n_total = len(_uniq_screened)
        if _n_total > 0:
            _ad_pct = 100 * _n_extrap / _n_total
            st.caption(
                f"{_n_extrap} of {_n_total} compounds ({_ad_pct:.1f}%) outside the "
                f"training domain (extrapolation)."
            )

    section_title("Docking Configuration",
                  "AutoDock Vina parameters.")

    _dock_labels = {
        "receptor_pdb_id": "Receptor PDB",
        "receptor_pdb_path": "Receptor file",
        "grid_center": "Grid centre (x, y, z)",
        "grid_size": "Grid box (Å)",
        "grid_box": "Grid box (Å)",
        "exhaustiveness": "Exhaustiveness",
        "num_modes": "Number of binding modes",
        "n_poses": "Number of poses",
        "seed": "Random seed",
        "scoring": "Scoring function",
        "scoring_function": "Scoring function",
        "cpu": "CPU count",
    }
    _dock_rows = []
    for _key, _label in _dock_labels.items():
        if _key in dock_cfg:
            _val = dock_cfg[_key]
            if isinstance(_val, (list, tuple)):
                _val = ", ".join(str(v) for v in _val)
            _dock_rows.append({"Parameter": _label, "Value": str(_val)})
    if _dock_rows:
        st.dataframe(
            pd.DataFrame(_dock_rows), hide_index=True,
            key="pipeline_dock_cfg", **_DATAFRAME_WIDTH_KW,
        )
    else:
        st.info("Docking configuration keys not exposed in "
                "`docking_config.json`.")

    st.markdown("---")
    sub_title("Pipeline Workflow")

    steps = [
        ("01", "ML Screening",
         f"Random Forest on ChEMBL MGAM (test PR-AUC {pr_auc})."),
        ("02", "Similarity",
         "Max Tanimoto to training actives."),
        ("03", "Developability",
         "Lipinski / Veber / PAINS / aggregator filters."),
        ("04", "Docking",
         f"Vina against PDB {pdb_id}."),
        ("05", "Consensus",
         "Weighted blend of ML, similarity, developability, novelty, docking."),
        ("06", "Applicability Domain",
         "ECFP4 Jaccard distance to the training set."),
    ]

    st.markdown("""
    <style>
    .pipeline-flow {
        display: flex;
        align-items: flex-start;
        justify-content: center;
        margin: 30px 0 14px 0;
        padding: 0 12px;
    }
    .pipeline-step {
        display: flex;
        flex-direction: column;
        align-items: center;
        flex: 0 0 auto;
        min-width: 92px;
    }
    .pipeline-node {
        width: 54px;
        height: 54px;
        border-radius: 50%;
        background: #173D32;
        color: #F7F6EF;
        display: flex;
        align-items: center;
        justify-content: center;
        font-family: 'DM Sans', sans-serif;
        font-weight: 700;
        font-size: 17px;
        letter-spacing: 0.6px;
        box-shadow: 0 0 0 2px #28634E, 0 2px 8px rgba(23, 61, 50, 0.18);
        margin-bottom: 12px;
        position: relative;
        z-index: 2;
    }
    .pipeline-label {
        font-size: 12px;
        font-weight: 600;
        color: #173D32;
        text-align: center;
        line-height: 1.35;
        max-width: 108px;
    }
    .pipeline-connector {
        flex: 1 1 auto;
        height: 2px;
        background: linear-gradient(90deg, #6B947A, #C7A66A);
        margin: 26px 4px 0 4px;
        min-width: 12px;
    }
    @media (max-width: 768px) {
        .pipeline-flow { padding: 0; }
        .pipeline-node { width: 42px; height: 42px; font-size: 14px; }
        .pipeline-label { font-size: 11px; max-width: 78px; }
        .pipeline-connector { margin-top: 20px; min-width: 6px; }
    }
    @media (max-width: 420px) {
        .pipeline-node { width: 36px; height: 36px; font-size: 12px; }
        .pipeline-label { font-size: 10px; max-width: 62px; }
        .pipeline-connector { margin-top: 17px; min-width: 4px; }
    }
    </style>
    """, unsafe_allow_html=True)

    stepper_parts = []
    for i, (num, title, _) in enumerate(steps):
        stepper_parts.append(
            f'<div class="pipeline-step">'
            f'<div class="pipeline-node">{esc(num)}</div>'
            f'<div class="pipeline-label">{esc(title)}</div>'
            f'</div>'
        )
        if i < len(steps) - 1:
            stepper_parts.append('<div class="pipeline-connector"></div>')

    st.markdown(
        f'<div class="pipeline-flow">{"".join(stepper_parts)}</div>',
        unsafe_allow_html=True,
    )

    st.markdown("---")

    sub_title("Stage Details")

    for i in range(0, len(steps), 2):
        col_a, col_b = st.columns(2)
        with col_a:
            num, title, desc = steps[i]
            info_card(f"{num} · {title}", desc, variant="info")
        if i + 1 < len(steps):
            with col_b:
                num, title, desc = steps[i + 1]
                info_card(f"{num} · {title}", desc, variant="info")


def _render_diag_section(df_section: pd.DataFrame, cols_candidates: list[str],
                         empty_msg: str, key_prefix: str) -> None:
    if df_section.empty:
        st.info(empty_msg)
        return
    cols = select_existing(df_section, cols_candidates)
    if not cols:
        st.warning("Diagnostics section has no recognised columns. "
                   f"Available: {', '.join(map(str, df_section.columns))}")
        return

    _search = st.text_input("Filter", key=f"{key_prefix}_search",
                            placeholder="Search any column…")
    _filtered = df_section
    if _search:
        mask = pd.Series(False, index=df_section.index)
        for c in cols:
            mask |= df_section[c].astype(str).str.contains(
                _search, case=False, na=False, regex=False)
        _filtered = df_section[mask]

    if _filtered.empty:
        st.caption("No rows match the filter.")
        return

    st.dataframe(_filtered[cols], hide_index=True, key=f"{key_prefix}_table",
                 **_DATAFRAME_WIDTH_KW)
    st.download_button(
        "⬇️ Download (CSV)",
        data=_filtered[cols].to_csv(index=False).encode("utf-8"),
        file_name=f"romarinova_{key_prefix}.csv",
        mime="text/csv", key=f"{key_prefix}_dl",
    )


def _page_evidence():
    page_hero(
        "Evidence Explorer",
        "Diagnostics, traceability, and out-of-scope candidates.",
        tags=["DIAGNOSTICS", "TRACEABILITY", "LIMITATIONS"],
    )

    def diag_section(name: str) -> pd.DataFrame:
        if "section" not in diag.columns:
            return pd.DataFrame()
        return diag[diag["section"] == name]

    _filtered_df = diag_section(DIAG_FILTERED)
    _docking_df = diag_section(DIAG_DOCKING)
    _missing_df = diag_section(DIAG_MISSING)
    _dis_df = diag_section(DIAG_DISAGREEMENT)

    filt_count = len(_filtered_df)
    dock_count = len(_docking_df)
    miss_count = len(_missing_df)
    dis_count = len(_dis_df)

    tabs = st.tabs([
        f"Filtered ({filt_count})",
        f"Docking Failures ({dock_count})",
        f"Missing Structures ({miss_count})",
        f"Disagreements ({dis_count})",
        "Limitations",
    ])

    with tabs[0]:
        sub_title("Filtered but Interesting")
        st.write("Failed ≥1 developability rule; kept for experimental follow-up.")

        if not _filtered_df.empty:
            _rule_col = next(
                (c for c in ("failed_rules", "reason", "failure_reason")
                 if c in _filtered_df.columns),
                None,
            )
            if _rule_col:
                _exploded = (
                    _filtered_df[_rule_col].fillna("").astype(str)
                    .str.split(r"[,;|]")
                    .explode().str.strip()
                )
                _exploded = _exploded[_exploded != ""]
                if not _exploded.empty:
                    _counts = (
                        _exploded.value_counts().rename_axis("rule")
                        .reset_index(name="count")
                        .sort_values("count", ascending=True)
                    )
                    fig_rule = px.bar(
                        _counts, x="count", y="rule", orientation="h",
                        color_discrete_sequence=["#C7A66A"], text="count",
                    )
                    fig_rule.update_traces(textposition="outside")
                    style_plotly(
                        fig_rule, title="Failures by developability rule",
                        xaxis=dict(title="Compounds", color="#24332C",
                                   gridcolor="#EDEFE9"),
                        yaxis=dict(title="", color="#24332C", automargin=True),
                        height=max(260, 35 * len(_counts) + 80),
                        showlegend=False,
                    )
                    st.plotly_chart(fig_rule, key="evidence_rule_bars",
                                    **_CHART_WIDTH_KW)
                    st.caption(
                        "Compounds can fail multiple rules; totals can exceed row count."
                    )

        _render_diag_section(
            _filtered_df,
            ["compound_id", "compound_name", "failed_rules", "final_score"],
            "No filtered candidates in diagnostics.",
            "diag_filtered",
        )

    with tabs[1]:
        sub_title("Docking Failures")

        if not _docking_df.empty:
            _reason_col = next(
                (c for c in ("failure_reason", "reason")
                 if c in _docking_df.columns),
                None,
            )
            if _reason_col:
                _reasons = (
                    _docking_df[_reason_col].fillna("unknown").astype(str)
                    .str.strip().replace("", "unknown")
                )
                _counts = (
                    _reasons.value_counts().rename_axis("reason")
                    .reset_index(name="count")
                    .sort_values("count", ascending=True)
                )
                if not _counts.empty:
                    fig_reason = px.bar(
                        _counts, x="count", y="reason", orientation="h",
                        color_discrete_sequence=["#28634E"], text="count",
                    )
                    fig_reason.update_traces(textposition="outside")
                    style_plotly(
                        fig_reason, title="Docking failures by reason",
                        xaxis=dict(title="Compounds", color="#24332C",
                                   gridcolor="#EDEFE9"),
                        yaxis=dict(title="", color="#24332C", automargin=True),
                        height=max(220, 35 * len(_counts) + 80),
                        showlegend=False,
                    )
                    st.plotly_chart(fig_reason, key="evidence_reason_bars",
                                    **_CHART_WIDTH_KW)

        _render_diag_section(
            _docking_df,
            ["compound_id", "compound_name", "failure_reason", "reason"],
            "All compounds docked successfully.",
            "diag_docking",
        )

    with tabs[2]:
        sub_title("Missing Structures")
        missing = _missing_df
        if missing.empty:
            _n = _to_int_or_none(meta_val("screened_compound_n", None))
            if _n:
                st.info(f"All {_n} compounds had defined 3D structures.")
            else:
                st.info("No missing structures recorded.")
        else:
            _render_diag_section(
                missing,
                ["compound_id", "compound_name", "reason", "failure_reason"],
                "None.", "diag_missing",
            )

    with tabs[3]:
        sub_title("Disagreement Patterns")
        _render_diag_section(
            _dis_df,
            ["compound_id", "compound_name", "pattern", "notes"],
            "No disagreement patterns recorded.",
            "diag_disagree",
        )

    with tabs[4]:
        sub_title("Limitations & Caveats")
        info_card(
            "No clinical claims",
            "Exploratory hypotheses only. Experimental validation and clinical trials required.",
            variant="warning",
        )
        info_card(
            "Docking is a hypothesis, not evidence",
            f"Vina against PDB {dock_cfg.get('receptor_pdb_id', '3L4W')} "
            f"produces scoring-function estimates, not binding free energies; "
            f"rigid receptors cannot capture induced fit.",
        )

        _n_screened = _to_int_or_none(meta_val("screened_compound_n", None))
        _n_extrap = _count_extrapolated_unique_compounds()
        if _n_screened and _n_extrap >= _n_screened:
            _ad_body = f"All {_n_screened} compounds sit outside the training domain."
        elif _n_screened:
            _ad_body = f"{_n_extrap} of {_n_screened} compounds sit outside the training domain."
        else:
            _ad_body = f"{_n_extrap} compound(s) sit outside the training domain."
        info_card("Applicability domain exceeded", _ad_body)

        info_card(
            "Data sources",
            "ChEMBL 37 (CC BY-SA 4.0), RDKit (BSD 3-Clause), PDB (Public Domain).",
        )
        st.caption("Full justification: **08 · Methodology**.")


def _page_methodology():
    page_hero(
        "Methodology & Justification",
        "Design decisions, evidence, and reproducibility notes.",
        tags=["REPRODUCIBILITY", "MODEL SELECTION", "APPLICABILITY DOMAIN"],
    )

    _n_act = _to_int_or_none(meta.get("n_actives_train"))
    _n_screened_meth = _to_int_or_none(meta_val("screened_compound_n", None)) or 0
    _n_extrap_meth = _count_extrapolated_unique_compounds()

    if _n_screened_meth and _n_extrap_meth >= _n_screened_meth:
        ad_conclusion = (
            f"All {_n_screened_meth} rosemary compounds are extrapolations."
        )
    elif _n_screened_meth:
        ad_conclusion = (
            f"{_n_extrap_meth} of {_n_screened_meth} rosemary compounds "
            f"are extrapolations."
        )
    else:
        ad_conclusion = "Rosemary compounds are flagged as extrapolations."

    if _n_act is not None:
        cv_intro = f"Training set is small ({_n_act} actives)."
    else:
        cv_intro = "Training set is small."

    section_title("Design Decisions", "Rationale behind each choice.")

    _decisions = [
        ("01 · Random Forest & ECFP4",
         "Random Forest with ECFP4 features (2048 bits, radius 2), selected over "
         "XGBoost by CV PR-AUC.",
         "info"),
        ("02 · Document-Aware Splitting",
         "Document-aware 70/15/15 split prevents paper-level leakage.",
         "success"),
        ("03 · Class-Aware Cross-Validation",
         f"{cv_intro} 4-fold class-aware CV with ≥2 actives per fold.",
         "info"),
        ("04 · Developability Filters",
         "Lipinski, Veber, PAINS, and aggregator filters remove non-drug-like "
         "candidates before ranking.",
         "warning"),
        ("05 · Rigid-Receptor Docking",
         "Vina against PDB 3L4W as an orthogonal structural filter; large "
         "diterpenes are correctly penalised.",
         "warning"),
        ("06 · Applicability Domain",
         "95th-percentile Jaccard distance to training set flags extrapolations. "
         + ad_conclusion,
         "success"),
    ]

    for i in range(0, len(_decisions), 2):
        col_a, col_b = st.columns(2)
        with col_a:
            t, d, v = _decisions[i]
            info_card(t, d, variant=v)
        if i + 1 < len(_decisions):
            with col_b:
                t, d, v = _decisions[i + 1]
                info_card(t, d, variant=v)

    st.markdown("---")

    section_title(
        "Model Selection Evidence",
        "Validation PR-AUC across candidate models; selected model highlighted.",
    )

    if not eval_table.empty and {"model", "split", "pr_auc"}.issubset(eval_table.columns):
        _val_df = eval_table[
            eval_table["split"].astype(str).str.lower().isin(["val", "validation"])
        ].dropna(subset=["pr_auc"])
        if _val_df.empty:
            _val_df = eval_table[eval_table["split"] == "test"].dropna(subset=["pr_auc"])

        if not _val_df.empty:
            _val_df = _val_df.sort_values("pr_auc", ascending=True)
            _colors = [
                "#28634E" if (_RF_MODEL and str(m) == _RF_MODEL) else "#A8B0AB"
                for m in _val_df["model"]
            ]
            fig_ms = go.Figure(go.Bar(
                x=_val_df["pr_auc"],
                y=_val_df["model"].astype(str),
                orientation="h",
                marker_color=_colors,
                text=[f"{v:.3f}" for v in _val_df["pr_auc"]],
                textposition="outside",
            ))
            style_plotly(
                fig_ms,
                title="PR-AUC comparison across candidate models",
                xaxis=dict(title="PR-AUC", color="#24332C", gridcolor="#EDEFE9"),
                yaxis=dict(title="", color="#24332C", automargin=True),
                height=max(CHART_SMALL, 45 * len(_val_df) + 100),
                showlegend=False,
                margin=dict(l=40, r=80, t=60, b=40),
            )
            st.plotly_chart(fig_ms, key="meth_model_cmp", **_CHART_WIDTH_KW)
            st.caption("PR-AUC preferred for rare actives; green bar = selected model.")
        else:
            st.info("No validation/test PR-AUC rows available.")
    else:
        st.info("Evaluation table unavailable — model-selection chart omitted.")

    st.markdown("---")

    section_title("Key References")
    col1, col2 = st.columns(2)
    with col1:
        info_card(
            "Primary sources",
            "ChEMBL 37 (Mendez et al., 2019); RDKit (Landrum, 2023); AutoDock Vina "
            "(Trott & Olson, 2010); ECFP4 (Rogers & Hahn, 2010); PAINS "
            "(Baell & Holloway, 2010).",
            variant="note",
        )
    with col2:
        info_card(
            "Reproducibility",
            "All consumed artefacts are regenerated by a single pipeline run.",
            variant="success",
        )


def _page_about():
    page_hero(
        "About ROMARINOVA",
        "Bridging Tunisian plant biodiversity and computational research.",
        tags=["OPEN SCIENCE", "ATTRIBUTION"],
    )

    _version = str(meta.get("version") or "1.1.0")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Powered by:** Romarinova AI")
        st.markdown("**Event:** EthnoHack 2026")
        st.markdown("**Track:** Natural Products for Drug Discovery")
    with c2:
        st.markdown(f"**Platform Version:** {_version}")
        st.markdown(f"**Build Date:** {BUILD_DATE}")
        st.markdown(f"**Contact:** [{CONTACT_URL}]({CONTACT_URL})")

    st.markdown("---")
    sub_title("Licenses & Attribution")
    st.markdown("- **ChEMBL:** CC BY-SA 4.0")
    st.markdown("- **RDKit:** BSD 3-Clause")
    st.markdown("- **AutoDock Vina:** Apache 2.0")
    st.markdown("- **PDB:** Public Domain")
    st.markdown("- **GBIF occurrence records:** CC BY-NC 4.0 / CC BY 4.0 per record")


# ==========================================================
# 4b. BIODIVERSITY MAP — helpers and page
# ==========================================================

_VARIETY_COLOURS = {
    "typicus":        [40, 99, 78, 200],
    "troglodytorum":  [199, 166, 106, 220],
    "laxiflorus":     [107, 148, 122, 200],
}
_PREDICTED_RING = [199, 166, 106, 255]
_OBSERVED_RING = [23, 61, 50, 255]

_BASIS_COLOURS = {
    "HUMAN_OBSERVATION":   [40, 99, 78, 200],
    "PRESERVED_SPECIMEN":  [199, 166, 106, 200],
    "OTHER":               [168, 176, 171, 160],
}

_CONCERN_ISSUES = (
    "PRESUMED_NEGATED_LONGITUDE",
    "COORDINATE_ROUNDED",
    "COUNTRY_DERIVED_FROM_COORDINATES",
)


def _flag_concerns(issue_str: str) -> str:
    if not issue_str:
        return ""
    flags = [f for f in _CONCERN_ISSUES if f in issue_str]
    return " · ".join(flags)


def _derive_variety(row) -> str | None:
    v = str(row.get("infraspecificEpithet") or "").strip().lower()
    if v in ("troglodytorum", "trogloditarum"):
        return "troglodytorum"
    if v in ("typicus", "typica"):
        return "typicus"
    if v == "laxiflorus":
        return "laxiflorus"
    for col in ("verbatimScientificName", "scientificName"):
        text = str(row.get(col) or "").lower()
        if "var. typicus" in text or "var.typicus" in text:
            return "typicus"
        if "troglodyt" in text or "troglodit" in text:
            return "troglodytorum"
        if "laxiflor" in text:
            return "laxiflorus"
    return None


def _predict_variety(lat, lon):
    """
    Assign a candidate variety from latitude/longitude, based on published
    distribution descriptions. THIS IS AN INFERENCE, NOT AN OBSERVATION.

    Zones:
      • Cap Bon (lat > 36.5, lon > 10.7): var. laxiflorus.
      • Southern Dahar / Matmata (lat < 34.5): var. troglodytorum.
      • Central-west highlands (34.5 ≤ lat < 35.5, lon < 9.5): var. troglodytorum.
      • Northern Tell (default): var. typicus.
    """
    if pd.isna(lat) or pd.isna(lon):
        return None
    if lat > 36.5 and lon > 10.7:
        return "laxiflorus"
    if lat < 34.5:
        return "troglodytorum"
    if 34.5 <= lat < 35.5 and lon < 9.5:
        return "troglodytorum"
    return "typicus"


def _build_tooltip(row) -> str:
    parts: list[str] = []
    sci = row.get("scientificName")
    parts.append(f"<b>{esc(sci)}</b>")

    variety = row.get("variety")
    source = row.get("variety_source")
    if isinstance(variety, str) and variety:
        if source == "observed":
            parts.append(f"Observed variety: <i>{esc(variety)}</i>")
        else:
            parts.append(
                f"<span style='color:#E0B670'>Predicted variety: "
                f"<i>{esc(variety)}</i> (inferred)</span>"
            )
    else:
        parts.append("Species-level (no variety assignment)")

    if not _is_missing(row.get("recordedBy")):
        parts.append(f"Recorded by: {esc(row.get('recordedBy'))}")
    if not _is_missing(row.get("stateProvince")):
        parts.append(f"Governorate: {esc(row.get('stateProvince'))}")
    if not _is_missing(row.get("locality")):
        loc = str(row.get("locality"))
        if len(loc) > 80:
            loc = loc[:77] + "…"
        parts.append(f"Locality: {esc(loc)}")
    if not _is_missing(row.get("year")):
        parts.append(f"Year: {int(row['year'])}")
    if not _is_missing(row.get("basisOfRecord")):
        parts.append(f"Basis: {esc(row.get('basisOfRecord'))}")
    if not _is_missing(row.get("institutionCode")):
        parts.append(f"Institution: {esc(row.get('institutionCode'))}")
    if row.get("is_country_centroid"):
        parts.append(
            "<span style='color:#E07B39'>⚠ country centroid — imprecise</span>"
        )
    if source == "predicted":
        parts.append(
            "<span style='color:#C7A66A;font-size:12px'>"
            "~ variety inferred from latitude, not observed</span>"
        )
    return "<br/>".join(parts)

def _read_occurrences_csv(path: Path) -> pd.DataFrame:
    """Robust reader for GBIF exports (tab, comma, semicolon, auto)."""
    attempts = [
        dict(sep="\t", engine="python"),
        dict(sep=",", engine="python"),
        dict(sep=";", engine="python"),
        dict(sep=None, engine="python"),
    ]
    last_exc: Exception | None = None

    for opts in attempts:
        try:
            df = pd.read_csv(
                path,
                dtype=str,
                encoding="utf-8-sig",
                on_bad_lines="skip",
                **opts,
            )
            if df.shape[1] >= 2:
                log.info(
                    "Loaded %s with sep=%r (%d rows, %d cols)",
                    path.name, opts.get("sep"), len(df), df.shape[1],
                )
                return df
        except Exception as exc:
            last_exc = exc
            continue

    if last_exc is not None:
        raise last_exc
    return pd.read_csv(
        path, dtype=str,
        encoding="utf-8-sig", on_bad_lines="skip",
    )

@st.cache_data(show_spinner="Loading GBIF occurrence records…")
def load_occurrences() -> pd.DataFrame:
    path = RESULTS_DIR / "rosemary_occurrences_tn.csv"
    if not path.exists():
        return pd.DataFrame()

    df = _read_occurrences_csv(path)

    for col in ("decimalLatitude", "decimalLongitude",
                "coordinateUncertaintyInMeters", "year", "month", "day"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["decimalLatitude", "decimalLongitude"]).copy()

    df["is_country_centroid"] = (
        (df["decimalLatitude"].round(2) == 34.0)
        & (df["decimalLongitude"].round(2) == 9.0)
    )

    df["variety_observed"] = df.apply(_derive_variety, axis=1)
    df["variety_predicted"] = df.apply(
        lambda r: _predict_variety(r["decimalLatitude"], r["decimalLongitude"]),
        axis=1,
    )
    df["variety"] = df["variety_observed"].fillna(df["variety_predicted"])
    df["variety_source"] = df["variety_observed"].apply(
        lambda v: "observed" if isinstance(v, str) and v else "predicted"
    )

    df["basisOfRecord"] = (
        df["basisOfRecord"].fillna("UNKNOWN") if "basisOfRecord" in df.columns
        else "UNKNOWN"
    )

    df["tooltip"] = df.apply(_build_tooltip, axis=1)

    return df


@st.cache_data(show_spinner=False)
def _occurrences_without_coords() -> pd.DataFrame:
    """Variety-level records lacking coordinates (herbarium sheets)."""
    path = RESULTS_DIR / "rosemary_occurrences_tn.csv"
    if not path.exists():
        return pd.DataFrame()
    df = _read_occurrences_csv(path)
    df["variety"] = df.apply(_derive_variety, axis=1)
    has_no_coords = df["decimalLatitude"].isna() | df["decimalLongitude"].isna()
    return df[has_no_coords & df["variety"].notna()].copy()


def _page_map():
    page_hero(
        "Biodiversity Map",
        "GBIF occurrences of Tunisian Salvia rosmarinus, coloured by inferred "
        "botanical variety. Full caveats at page bottom.",
        tags=["GBIF", "INFERRED VARIETIES", "SAMPLING EFFORT"],
    )

    df = load_occurrences()
    if df.empty:
        st.warning(
            "No occurrence file found at `results/rosemary_occurrences_tn.csv`. "
            "Download a GBIF CSV for *Salvia rosmarinus* (`country = TN`) "
            "and save it there."
        )
        st.link_button(
            "Open GBIF search — Tunisia",
            "https://www.gbif.org/occurrence/search?country=TN&taxon_key=10902460",
        )
        return


    f1, f2, f3 = st.columns([1, 1.2, 1.3])

    with f1:
        _basis_opts = sorted(df["basisOfRecord"].dropna().unique().tolist())
        basis_filter = st.multiselect(
            "Basis of record", options=_basis_opts, default=_basis_opts,
            key="map_basis_filter",
        )
    with f2:
        _var_opts = sorted(df["variety"].dropna().unique().tolist())
        var_filter = st.multiselect(
            "Variety (predicted)", options=_var_opts, default=_var_opts,
            key="map_var_filter",
        )
    with f3:
        _yrs = df["year"].dropna()
        year_range = None
        if not _yrs.empty:
            year_range = st.slider(
                "Collection year",
                min_value=int(_yrs.min()), max_value=int(_yrs.max()),
                value=(int(_yrs.min()), int(_yrs.max())),
                key="map_year_range",
            )

    _view = df.copy()
    if basis_filter:
        _view = _view[_view["basisOfRecord"].isin(basis_filter)]
    if var_filter:
        _view = _view[_view["variety"].isin(var_filter)]
    if year_range is not None:
        _view = _view[
            _view["year"].isna()
            | ((_view["year"] >= year_range[0]) & (_view["year"] <= year_range[1]))
        ]
    _view = _view[~_view["is_country_centroid"]]

    if _view.empty:
        st.warning("No occurrences match the current filters.")
        return

    _view["fill_colour"] = _view["variety"].map(_VARIETY_COLOURS).apply(
        lambda c: c if isinstance(c, list) else [168, 176, 171, 180]
    )
    _view["ring_colour"] = _view["variety_source"].apply(
        lambda s: _OBSERVED_RING if s == "observed" else _PREDICTED_RING
    )

    m1, m2, m3, m4 = st.columns(4)
    with m1:
        metric("Records on map", f"{len(_view)}",
               help_text="Records after filters, with usable coordinates.")
    with m2:
        metric("Variety determinations",
               f"{int((df['variety_source'] == 'observed').sum())}",
               help_text="Herbarium sheets with a variety label (none have GPS).")
    with m3:
        metric("Inferred assignments",
               f"{int((df['variety_source'] == 'predicted').sum())}",
               help_text="Points whose variety is inferred from coordinates.")
    with m4:
        metric("Institutions",
               f"{df['institutionCode'].dropna().nunique() if 'institutionCode' in df else 0}",
               help_text="Publishing herbaria and observation platforms.")

    _show_heatmap = st.toggle(
        "Density heatmap", value=True, key="map_heatmap_toggle",
        help="Sampling effort, not true abundance.",
    )

    layers = []
    if _show_heatmap:
        layers.append(pdk.Layer(
            "HeatmapLayer", data=_view,
            get_position=["decimalLongitude", "decimalLatitude"],
            radiusPixels=60, intensity=1.0, threshold=0.05, aggregation="SUM",
        ))
    layers.append(pdk.Layer(
        "ScatterplotLayer", data=_view,
        get_position=["decimalLongitude", "decimalLatitude"],
        get_fill_color="fill_colour",
        stroked=True,
        get_line_color="ring_colour",
        get_line_width=2,
        radius_min_pixels=6,
        radius_max_pixels=28,
        radius_scale=1,
        radius=3500,
        pickable=True, auto_highlight=True,
    ))

    view_state = pdk.ViewState(
        latitude=float(_view["decimalLatitude"].mean()),
        longitude=float(_view["decimalLongitude"].mean()),
        zoom=6.2, pitch=0,
    )

    tooltip = {
        "html": "{tooltip}",
        "style": {
            "backgroundColor": "#173D32", "color": "#F7F6EF",
            "fontFamily": "DM Sans, sans-serif", "fontSize": "13.5px",
            "borderRadius": "10px", "padding": "12px",
        },
    }

    st.pydeck_chart(
        pdk.Deck(layers=layers, initial_view_state=view_state, tooltip=tooltip),
        **_CHART_WIDTH_KW,
    )

    legend_items = [
        ("var. typicus (inferred)",       "#28634E"),
        ("var. troglodytorum (inferred)", "#C7A66A"),
        ("var. laxiflorus (inferred)",    "#6B947A"),
    ]

    legend_html = "".join(
        f'<span style="display:inline-block;margin-right:18px;">'
        f'<span style="display:inline-block;width:14px;height:14px;'
        f'border-radius:50%;background:{c};border:2px solid #C7A66A;'
        f'vertical-align:middle;margin-right:7px;"></span>{esc(label)}</span>'
        for label, c in legend_items
    )
    ring_note = (
        '<span style="display:inline-block;margin-left:24px;color:#8A6F3B;">'
        '<span style="display:inline-block;width:14px;height:14px;'
        'border-radius:50%;border:2px solid #C7A66A;background:transparent;'
        'vertical-align:middle;margin-right:6px;"></span>'
        'gold ring = inferred from coordinates'
        '</span>'
    )

    st.markdown(
        f'<div style="font-size:13px;color:#4A5550;margin-top:10px;">'
        f'{legend_html}{ring_note}</div>',
        unsafe_allow_html=True,
    )

    st.caption(
        "**Sampling effort, not abundance.** Bright cells mark collection hotspots, "
        "not rosemary density. Two-thirds of modern records fall within 50 km of "
        "Tunis; empty cells ≠ absence."
    )

    _n_observed = int((df["variety_source"] == "observed").sum())
    _n_predicted = int((df["variety_source"] == "predicted").sum())

    st.markdown(
        f"""
        <div style="background:#FFF8E6;border-left:5px solid #C7A66A;
                    border-radius:12px;padding:14px 20px;margin-bottom:18px;
                    color:#24332C;font-size:14.5px;line-height:1.5;">
            <b style="color:#173D32;">⚠ Variety colours are inferred, not observed.</b>
            All {_n_predicted} mappable points use a coordinate heuristic; the
            {_n_observed} actual determinations lack GPS. Details in
            <i>“Known data issues”</i>.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("---")
    sub_title("Temporal and Institutional Coverage")

    _yr_df = (
        _view.dropna(subset=["year"])
        .assign(year=lambda d: d["year"].astype(int))
        .groupby(["year", "basisOfRecord"]).size().reset_index(name="count")
    )
    if not _yr_df.empty:
        fig_yr = px.bar(
            _yr_df, x="year", y="count", color="basisOfRecord", barmode="stack",
            color_discrete_map={
                "HUMAN_OBSERVATION": "#28634E",
                "PRESERVED_SPECIMEN": "#C7A66A",
            },
            labels={"year": "Year", "count": "Records", "basisOfRecord": "Basis"},
        )
        style_plotly(
            fig_yr, title="Records by collection year and basis of record",
            xaxis=dict(title="Year", color="#24332C", gridcolor="#EDEFE9"),
            yaxis=dict(title="Records", color="#24332C", gridcolor="#EDEFE9"),
            height=CHART_MEDIUM, legend=LEGEND_TOP_RIGHT,
        )
        st.plotly_chart(fig_yr, key="map_year_chart", **_CHART_WIDTH_KW)
        st.caption(
            "Bimodal: historic herbaria (1907–1980) and a post-2007 surge from "
            "iNaturalist / observation.org, clustered around Tunis, Nabeul, Sousse."
        )

    st.markdown("---")
    _variety_df = _occurrences_without_coords()
    if _variety_df.empty:
        st.info("No variety-level determinations found.")
    else:
        with st.expander(
            f"Variety-level determinations without coordinates "
            f"({len(_variety_df)})",
        ):
            st.caption(
                "The only true variety records in this export — none carry GPS. "
                "The holotype of var. *troglodytorum* (Maire & Weiller, 1938) is "
                "among them, so neither its locality nor its taxonomic status can "
                "be mapped — a recurring limitation of historical botany in "
                "North Africa."
            )
            _vcols = select_existing(_variety_df, [
                "gbifID", "variety", "scientificName", "recordedBy",
                "year", "locality", "institutionCode", "catalogNumber",
                "typeStatus",
            ])
            st.dataframe(
                _variety_df[_vcols], hide_index=True,
                key="map_variety_table", **_DATAFRAME_WIDTH_KW,
            )

    st.markdown("---")
    with st.expander("⚠ Known data issues and modelling choices"):
        st.markdown(
            "- **Variety colours are heuristics.** No GBIF record here pairs a "
            "variety label with GPS. The variety layer is a bioclimatic "
            "inference, not evidence.\n"
            "- **Conservation context.** The two endemic varieties occupy "
            "different bioclimates; var. *troglodytorum* is considered "
            "threatened. Read records alongside protected-area coverage.\n"
            "- **Country centroids.** Two MNHNL records use (34.0, 9.0) as a "
            "placeholder — always excluded.\n"
            "- **BOM.** The raw GBIF CSV begins with a UTF-8 byte-order mark; "
            "handled via `encoding='utf-8-sig'`.\n"
            "- **Misspelled epithet.** `infraspecificEpithet` contains "
            "`trogloditarum`; normalised to `troglodytorum`.\n"
            "- **A fourth variety omitted.** `var. lavandulaceum` overlaps "
            "almost entirely with *troglodytorum* and cannot be resolved by "
            "coordinate — excluded.\n"
            "- **Herbaria only for variety.** All iNaturalist records are "
            "species-level; amateur platforms do not capture infra-specific "
            "taxonomy.\n"
            "- **Spatial bias.** Two-thirds of modern records fall within "
            "50 km of Tunis; the south is under-recorded.\n"
            "- **How to validate.** A specimen collected within a point's "
            "buffer, with a taxonomic determination, would confirm or refute "
            "the heuristic — the natural next experiment."
        )


# ==========================================================
# 5. NAVIGATION (st.navigation / st.Page)
# ==========================================================

_StopException: type | None = None
for _mod_path in ("streamlit.runtime.scriptrunner", "streamlit.scriptrunner"):
    try:
        _mod = __import__(_mod_path, fromlist=["StopException"])
        _StopException = getattr(_mod, "StopException", None)
        if _StopException is not None:
            break
    except ImportError:
        continue


def _is_stop_exception(e: BaseException) -> bool:
    if _StopException is not None and isinstance(e, _StopException):
        return True
    return type(e).__name__ == "StopException"


def _safe_render(page_fn: Callable[[], None], name: str):
    try:
        page_fn()
    except Exception as e:
        if _is_stop_exception(e):
            raise
        log.exception("Page %s failed", name)
        st.error(f"An error occurred rendering **{name}**: `{e}`")
        with st.expander("Traceback"):
            st.code(traceback.format_exc())


_PAGE_DASHBOARD = st.Page(
    lambda: _safe_render(_page_dashboard, "Dashboard"),
    title="01 · Dashboard",
    icon=":material/dashboard:",
    url_path="dashboard",
    default=True,
)
_PAGE_MAP = st.Page(
    lambda: _safe_render(_page_map, "Biodiversity Map"),
    title="02 · Biodiversity Map",
    icon=":material/public:",
    url_path="map",
)
_PAGE_ATLAS = st.Page(
    lambda: _safe_render(_page_atlas, "Chemical Atlas"),
    title="03 · Chemical Atlas",
    icon=":material/biotech:",
    url_path="atlas",
)
_PAGE_SHORTLIST = st.Page(
    lambda: _safe_render(_page_shortlist, "α-Glucosidase Screening"),
    title="04 · α-Glucosidase Screening",
    icon=":material/medication:",
    url_path="shortlist",
)
_PAGE_TARGETS = st.Page(
    lambda: _safe_render(_page_targets, "Biological Targets"),
    title="05 · Biological Targets",
    icon=":material/account_tree:",
    url_path="targets",
)
_PAGE_PIPELINE = st.Page(
    lambda: _safe_render(_page_pipeline, "Bioinformatics & Docking"),
    title="06 · Bioinformatics & Docking",
    icon=":material/route:",
    url_path="pipeline",
)
_PAGE_EVIDENCE = st.Page(
    lambda: _safe_render(_page_evidence, "Scientific Evidence"),
    title="07 · Scientific Evidence",
    icon=":material/library_books:",
    url_path="evidence",
)
_PAGE_METHODOLOGY = st.Page(
    lambda: _safe_render(_page_methodology, "Methodology"),
    title="08 · Methodology",
    icon=":material/science:",
    url_path="methodology",
)
_PAGE_ABOUT = st.Page(
    lambda: _safe_render(_page_about, "About ROMARINOVA"),
    title="09 · About ROMARINOVA",
    icon=":material/info:",
    url_path="about",
)


@st.fragment
def _sidebar_search():
    query = st.text_input(
        "🔍 Search compound",
        placeholder="Name, ID, or InChIKey",
        key="_sidebar_search_query",
    )
    if not query or len(query.strip()) < 2:
        return

    base = screened.dropna(subset=["compound_id"]).drop_duplicates("compound_id")
    q = query.strip().lower()

    mask = base["compound_name"].astype(str).str.lower().str.contains(
        q, na=False, regex=False,
    )
    mask |= base["compound_id"].astype(str).str.lower().str.contains(
        q, na=False, regex=False,
    )
    if "inchikey" in base.columns:
        mask |= base["inchikey"].astype(str).str.lower().str.contains(
            q, na=False, regex=False,
        )

    total = int(mask.sum())
    if total == 0:
        st.caption(f"No match for “{query.strip()}”.")
        return

    top = base[mask].head(6)
    if total > len(top):
        st.caption(f"{total} matches — top {len(top)} shown.")
    else:
        st.caption(f"{total} match{'es' if total != 1 else ''}.")

    current = st.session_state.get("atlas_compound")
    for _, row in top.iterrows():
        cid = row["compound_id"]
        raw_name = row["compound_name"]
        name = (str(raw_name).strip()
                if not _is_missing(raw_name) and str(raw_name).strip()
                else f"Unnamed ({cid})")
        trunc = name if len(name) <= 26 else name[:23] + "…"
        label = f"{'● ' if cid == current else ''}→ {trunc}"
        if st.button(label, key=f"_srch_btn_{cid}", **_BUTTON_WIDTH_KW):
            st.session_state["atlas_compound"] = cid
            st.switch_page(_PAGE_ATLAS)


_PAGES_LIST = [
    _PAGE_DASHBOARD, _PAGE_MAP, _PAGE_ATLAS, _PAGE_SHORTLIST, _PAGE_TARGETS,
    _PAGE_PIPELINE, _PAGE_EVIDENCE, _PAGE_METHODOLOGY, _PAGE_ABOUT,
]

pg = st.navigation(_PAGES_LIST, position="hidden")


with st.sidebar:
    logo_path = find_logo()
    if logo_path:
        _logo_uri = _logo_data_uri(str(logo_path))
        if _logo_uri:
            st.markdown(
                f'<img class="romarinova-logo" src="{_logo_uri}" '
                f'alt="ROMARINOVA logo"/>',
                unsafe_allow_html=True,
            )
        else:
            st.markdown("## 🌿 ROMARINOVA")
    else:
        st.markdown("## 🌿 ROMARINOVA")

    st.markdown("---")

    for _page in _PAGES_LIST:
        st.page_link(_page)

    st.markdown("---")

    _sidebar_search()

    st.markdown("---")
    st.caption("Powered by Romarinova AI · EthnoHack 2026 · Track 2")

    st.markdown("---")

    _issues: list[str] = []
    if _filled_from_derivation:
        _issues.append(
            "Pipeline metadata missing or inconsistent with artefacts "
            f"({', '.join(sorted(_filled_from_derivation))}). "
            "Values corrected at runtime."
        )
    _issues.extend(_CONSISTENCY_WARNINGS)
    _issues.extend(_ENUM_ISSUES)

    if _issues:
        st.markdown(
            f'<div class="health-pill warn">⚠️ {len(_issues)} data-health '
            f'note{"s" if len(_issues) != 1 else ""}</div>',
            unsafe_allow_html=True,
        )
        with st.expander("Data-health details"):
            for msg in _issues:
                st.caption(f"• {msg}")
            st.caption(
                "Rerun `scripts/fix_shortlist_metadata.py` to regenerate a clean "
                "`shortlist_metadata.json`."
            )

    if st.button("🔄 Reload data", key="reload_btn", **_BUTTON_WIDTH_KW):
        for _fn in (load_data, mol_image, render_pdb_html, _logo_data_uri,
                    _artefact_mtime, _compute_descriptors,
                    load_occurrences, _occurrences_without_coords):
            try:
                _fn.clear()
            except Exception:
                pass
        st.session_state["_reload_toast_pending"] = True
        st.rerun()

    if st.session_state.pop("_reload_toast_pending", False):
        st.toast("Data reloaded from disk", icon="✅")


pg.run()


# ==========================================================
# 6. FOOTER
# ==========================================================

_artefact_ts = _artefact_mtime() or "unknown"

st.markdown(f"""
<div class="footer">
    <b>ROMARINOVA</b> · Scientific Platform
    <br>
    Powered by Romarinova AI · EthnoHack 2026
    <br>
    Artefacts updated: {esc(_artefact_ts)} · Build: {esc(BUILD_DATE)}
    <br>
    Hypotheses only — experimental validation required.
</div>
""", unsafe_allow_html=True)