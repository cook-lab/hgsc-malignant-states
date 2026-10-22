#!/usr/bin/env python3
"""
environment.yml from the Python environment that ran the analyses (environment E2, env/README.md).

PURPOSE
    Pins every third-party package the repository's Python imports, plus everything those packages
    require (recursive Requires-Dist, markers evaluated for this machine, extras skipped), at the
    versions installed in the environment that runs this script. Run it with the python of the
    conda base env that ran all final-run Python (Python 3.12.9; the scientific stack there was
    installed with pip, so it is pinned under `pip:`).
    Imports are read from the tracked .py files (ast) and notebooks (code cells); repo-local modules
    and the standard library are skipped. The full export of the same env (every package, including
    tools the repo never imports) is written by env/export_conda_envs.sh.

INPUTS
    the running interpreter's environment
    the repo's tracked *.py / *.ipynb files (git ls-files), minus env/

OUTPUTS
    environment.yml (env name epitype-py)
    stdout: imports not installed in this environment; pinned packages whose requirements the pinned
            set does not meet

MANUSCRIPT PANEL(S): none (environment record)
RUNTIME TIER: fast (seconds)

USAGE (from the repo root)
    <python of the analysis env> env/make_environment_yml.py
"""
import ast, datetime, json, platform, re, subprocess, sys
import importlib.metadata as md
from pathlib import Path
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "environment.yml"
# Packages other packages load on demand, so the code never imports them
ON_DEMAND = ["leidenalg"]   # scanpy.tl.leiden (igraph comes with it)


def tracked(*patterns):
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--", *patterns],
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.split("\0") if p and not p.startswith("env/")]


def imports_in(src):
    """Top-level module names imported by a block of Python source."""
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return set(re.findall(r"^\s*(?:import|from)\s+([A-Za-z_]\w*)", src, flags=re.M))
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            mods.add(node.module.split(".")[0])
        elif (isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant)
              and isinstance(node.args[0].value, str)
              and ((isinstance(node.func, ast.Attribute) and node.func.attr == "import_module")
                   or (isinstance(node.func, ast.Name) and node.func.id == "__import__"))):
            mods.add(node.args[0].value.split(".")[0])
    return mods


# ---- modules the repo imports ----------------------------------------------------------
py_files, nb_files = tracked("*.py"), tracked("*.ipynb")
local = {Path(f).stem for f in py_files} | {p for f in py_files for p in Path(f).parent.parts}
used = {}
for f in py_files:
    for m in imports_in((REPO / f).read_text(encoding="utf-8", errors="replace")):
        used.setdefault(m, set()).add(f)
for f in nb_files:
    nb = json.loads((REPO / f).read_text(encoding="utf-8"))
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        src = "".join(cell.get("source", []))
        src = "\n".join(l for l in src.splitlines() if not l.lstrip().startswith(("%", "!")))
        for m in imports_in(src):
            used.setdefault(m, set()).add(f)
third_party = sorted(m for m in used if m not in sys.stdlib_module_names
                     and m not in local and m != "__future__")

# ---- map to installed distributions ----------------------------------------------------
mod2dist = md.packages_distributions()
dists = {canonicalize_name(d.metadata["Name"]): d for d in md.distributions()}
roots, not_installed = {canonicalize_name(n) for n in ON_DEMAND if canonicalize_name(n) in dists}, []
for m in third_party:
    names = mod2dist.get(m)
    if names:
        roots.update(canonicalize_name(n) for n in names)
    else:
        not_installed.append(m)

# where the missing ones are installed: other envs of this conda install
elsewhere = {}
for site in sorted(Path(sys.prefix, "envs").glob("*/lib/python3*/site-packages")):
    env = site.parts[-4]
    for d in md.distributions(path=[str(site)]):
        top = (d.read_text("top_level.txt") or "").split() or [d.metadata["Name"].replace("-", "_")]
        for m in not_installed:
            if m in top and m not in elsewhere:
                elsewhere[m] = f"{d.metadata['Name']}=={d.version} in conda env {env}"

# ---- recursive requirements ----------------------------------------------------------------
pinned, queue = set(), sorted(roots)
while queue:
    name = queue.pop()
    if name in pinned or name not in dists:
        continue
    pinned.add(name)
    for r in dists[name].requires or []:
        req = Requirement(r)
        if req.marker and not req.marker.evaluate({"extra": ""}):
            continue
        queue.append(canonicalize_name(req.name))

pip_roots = [f"{dists[n].metadata['Name']}=={dists[n].version}" for n in sorted(pinned) if n in roots]
pip_deps = [f"{dists[n].metadata['Name']}=={dists[n].version}" for n in sorted(pinned) if n not in roots]

# ---- consistency of the pinned set ------------------------------------------------------
conflicts = []
for name in sorted(pinned):
    for r in dists[name].requires or []:
        req = Requirement(r)
        if req.marker and not req.marker.evaluate({"extra": ""}):
            continue
        dep = canonicalize_name(req.name)
        if dep in pinned and req.specifier and not req.specifier.contains(dists[dep].version, prereleases=True):
            conflicts.append(f"{name} {dists[name].version} requires {req}, pinned {dists[dep].version}")

# ---- write ------------------------------------------------------------------------------
py = platform.python_version()
lines = [
    f"# Written by env/make_environment_yml.py on {datetime.date.today()} -- do not edit by hand.",
    f"# The Python environment that ran the analyses (E2 in env/README.md): Python {py},",
    "# every third-party package the repo's Python imports plus their requirements, at the installed",
    "# versions. That env's scientific stack was installed with pip, hence the pip: block (a few",
    "# packages there came from conda; pip installs the same versions). Full export of the env:",
    "# env/E2_base_full.yml. Atlas QC and doublet scoring ran in E3 (env/E3_ovarian_atlas_full.yml);",
    "# scVI/scANVI integration on a GPU cluster with scvi-tools 1.4.1, the version pinned below.",
    f"# Pins were taken on {platform.system()} {platform.machine()}.",
    "# Imported by the repo but not in this env:",
]
lines += [f"#   {m:14s} {elsewhere.get(m, 'not installed in any local env')}" for m in not_installed]
lines += [f"#   NOTE: {c}" for c in conflicts]
lines += ["name: epitype-py", "channels:", "  - conda-forge", "dependencies:", f"  - python={py}"]
lines += ["  - pip", "  - pip:", "      # imported by the repo, or loaded on demand by it (" + ", ".join(ON_DEMAND) + ")"]
lines += [f"      - {p}" for p in sorted(pip_roots, key=str.lower)]
lines += ["      # required by the packages above"]
lines += [f"      - {p}" for p in sorted(pip_deps, key=str.lower)]
OUT.write_text("\n".join(lines) + "\n")

print(f"Python {py} at {sys.prefix}")
print(f"{len(third_party)} third-party modules imported; {len(roots)} distributions with on-demand ones; "
      f"{len(pinned)} pinned with requirements")
print("imported but not installed here:")
for m in not_installed:
    print(f"  {m:14s} {', '.join(sorted(used[m])[:3])}  [{elsewhere.get(m, 'not installed in any local env')}]")
print("pinned set conflicts:", "none" if not conflicts else "")
for c in conflicts:
    print("  " + c)
