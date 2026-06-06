"""Central config loader for Python (atlas) scripts.

Usage:
    from config.config import CFG, path, obj, SEED
    adata = sc.read_h5ad(obj("atlas_final"))      # resolves data_root + object relpath
    out = path("output_root", "03_epithelial_nmf")
    np.random.seed(SEED)

Resolves ${VAR:-default} expansions against the environment so DATA_ROOT /
OUTPUT_ROOT can be overridden without editing config.yml.
"""
from __future__ import annotations
import os
import re
from pathlib import Path

import yaml

_CFG_PATH = Path(__file__).resolve().parent / "config.yml"
_ENV_RE = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(?::-([^}]*))?\}")


def _expand(value):
    if isinstance(value, str):
        def repl(m):
            return os.environ.get(m.group(1), m.group(2) or "")
        return _ENV_RE.sub(repl, value)
    if isinstance(value, dict):
        return {k: _expand(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand(v) for v in value]
    return value


with open(_CFG_PATH) as fh:
    CFG = _expand(yaml.safe_load(fh))

SEED: int = int(CFG.get("seed", 42))
DATA_ROOT = Path(CFG["paths"]["data_root"]).expanduser()
OUTPUT_ROOT = Path(CFG["paths"]["output_root"]).expanduser()


def obj(key: str) -> str:
    """Absolute path to an entry-point object by its key in config.objects."""
    return str(DATA_ROOT / CFG["objects"][key])


def path(root_key: str, *parts) -> str:
    """Build a path under a configured root ('data_root' or 'output_root'); makes parent dirs for outputs."""
    base = Path(CFG["paths"][root_key]).expanduser()
    p = base.joinpath(*map(str, parts))
    # auto-create dirs for any OUTPUT-side root (output_root, figures_dir, …); never for
    # data_root (read-only input). If the path looks like a directory (no file extension)
    # create it; if it looks like a file create its parent.
    if root_key != "data_root":
        target = p.parent if p.suffix else p
        target.mkdir(parents=True, exist_ok=True)
    return str(p)
