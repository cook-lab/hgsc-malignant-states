"""Study codes for samples and patients (de-identification).

Everything published uses study codes: WT01-WT08 (whole tissues), TMA-P### (TMA patients) and FTE##
(fallopian-tube cores). Data held in the lab may still carry the source identifiers; the mapping lives in a
private key file that is never stored in this repository (path: $STUDY_CODE_KEY, by default
$OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv).

code_map(arm) returns {identifier: study code} for one arm and also maps every code to itself, so scripts
run unchanged on the lab's data and on the deposited, already-coded data. Without the key file only the
identity part is available: any name shaped like a code of that arm maps to itself (`in`, [], .get and
pandas .map all work).

    from config.study_codes import code_map, code_of
    code_of("Whole tissue", name)        # name may be a source ID or a code
    code_of("TMA A", patient_id)         # TMA patient numbers are looked up within their TMA
    code_of("FTE core", core_id)         # FTE cores are looked up by core id
"""
import os
import re

import pandas as pd

KEY_PATH = os.environ.get("STUDY_CODE_KEY", os.path.join(
    os.environ.get("OUTPUT_ROOT", "output"), "private_keys", "SuppTable8_code_key_PRIVATE.csv"))
CODE_RX = {"Whole tissue": r"WT\d\d", "TMA A": r"TMA-P\d{3}", "TMA B": r"TMA-P\d{3}", "FTE core": r"FTE\d\d"}
_KEY = None


def key():
    """The private key as a DataFrame (arm, study_code, source_id); empty when the file is absent."""
    global _KEY
    if _KEY is None:
        _KEY = (pd.read_csv(KEY_PATH, dtype=str) if os.path.exists(KEY_PATH)
                else pd.DataFrame(columns=["arm", "study_code", "source_id"]))
    return _KEY


def _source_key(arm, source_id):
    # FTE cores are keyed by core id ("patient N / core M" -> "M"); other arms by the identifier itself
    return source_id.split("core ")[-1] if arm == "FTE core" else source_id


class CodeMap(dict):
    """{identifier: code} from the key; any code of the arm also maps to itself, with or without the key."""
    def __init__(self, arm, m):
        super().__init__(m)
        self.rx = re.compile(CODE_RX[arm])

    def _is_code(self, x):
        return isinstance(x, str) and self.rx.fullmatch(x) is not None

    def __contains__(self, x):
        return dict.__contains__(self, x) or self._is_code(x)

    def __missing__(self, x):      # used by [] and by pandas .map
        if self._is_code(x):
            return x
        raise KeyError(x)

    def get(self, x, default=None):
        return self[x] if x in self else default


def code_map(arm):
    k = key()
    k = k[k.arm == arm]
    m = {_source_key(arm, s): c for s, c in zip(k.source_id, k.study_code)}
    m.update({c: c for c in list(m.values())})
    return CodeMap(arm, m)


def code_of(arm, x):
    x = str(x)
    if re.fullmatch(r"\d+\.0", x):          # numbers read back as floats
        x = x[:-2]
    m = code_map(arm)
    if x in m:
        return m[x]
    if re.fullmatch(CODE_RX[arm], x):
        return x
    raise KeyError(f"no study code for {arm} {x!r} (is the private key at $STUDY_CODE_KEY?)")
