#!/usr/bin/env python3
"""
Correct metastatic_site in the master atlas so it reflects anatomic_site
(authoritative) for vazquez_garcia_2022.

Within vazquez the only discrepancy is 7 samples with anatomic_site == 'ascites'
that were coded metastatic_site == 'metastasis'. anatomic_site is correct, so
these are reassigned metastatic_site -> 'ascites'. Every other vazquez sample
already matches its anatomic site.

The edit is surgical: obs/metastatic_site is a categorical stored as int8
codes; we flip the affected codes from 'metastasis' -> 'ascites' in place. The
expression matrix (X) is never touched.

A backup of the ORIGINAL metastatic_site codes + categories is written first, so
the change is reversible.

Run:  python fix_metastatic_site_vazquez.py           # dry-run, reports only
      python fix_metastatic_site_vazquez.py --apply   # perform the edit
"""
import argparse
import os
import numpy as np
import h5py

ROOT = "/Volumes/CookLab/Sarah/202601_final_atlas"
ATLAS = os.path.join(ROOT, "hgsc_atlas_final.h5ad")
BACKUP = os.path.join(ROOT, "sarahs_side_quests/data/metastatic_site_ORIGINAL_backup.npz")

STUDY = "vazquez_garcia_2022"
ANAT = "ascites"
FROM = "metastasis"
TO = "ascites"


def cats_codes(h5grp):
    cats = np.array([c.decode() if isinstance(c, bytes) else c
                     for c in h5grp["categories"][:]])
    codes = h5grp["codes"][:]
    return cats, codes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="actually write the change (default: dry-run)")
    args = ap.parse_args()

    mode = "r+" if args.apply else "r"
    f = h5py.File(ATLAS, mode)
    try:
        s_cats, s_codes = cats_codes(f["obs/study"])
        a_cats, a_codes = cats_codes(f["obs/anatomic_site"])
        m_cats, m_codes = cats_codes(f["obs/metastatic_site"])

        def code_of(cats, val):
            hit = np.where(cats == val)[0]
            if len(hit) != 1:
                raise SystemExit(f"category {val!r} not found uniquely in {cats}")
            return hit[0]

        s_i, a_i, from_i, to_i = (code_of(s_cats, STUDY), code_of(a_cats, ANAT),
                                  code_of(m_cats, FROM), code_of(m_cats, TO))

        mask = (s_codes == s_i) & (a_codes == a_i) & (m_codes == from_i)
        n = int(mask.sum())
        print(f"cells matching [{STUDY} & anat={ANAT} & met={FROM}]: {n:,}")
        if n == 0:
            print("Nothing to change (already corrected?). Exiting.")
            return

        # sanity: confirm these are exactly the expected samples via sample_id
        sid_cats, sid_codes = cats_codes(f["obs/sample_id"])
        aff_samples = sorted(set(sid_cats[sid_codes[mask]]))
        print(f"affected samples ({len(aff_samples)}): {aff_samples}")

        if not args.apply:
            print("\n[dry-run] re-run with --apply to write. No changes made.")
            return

        # backup original codes + categories
        np.savez_compressed(BACKUP, codes=m_codes, categories=m_cats,
                            note=f"original obs/metastatic_site before "
                                 f"{STUDY} {ANAT} {FROM}->{TO} fix")
        print(f"backup written -> {BACKUP}")

        # surgical edit
        new_codes = m_codes.copy()
        new_codes[mask] = to_i
        f["obs/metastatic_site/codes"][...] = new_codes
        f.flush()
        print(f"applied: {n:,} cells recoded {FROM} -> {TO}")

        # verify
        _, chk = cats_codes(f["obs/metastatic_site"])
        remaining = int(((s_codes == s_i) & (a_codes == a_i) &
                         (chk == from_i)).sum())
        now_asc = int(((s_codes == s_i) & (a_codes == a_i) &
                       (chk == to_i)).sum())
        print(f"verify: vazquez-ascites cells now '{TO}'={now_asc:,}, "
              f"still '{FROM}'={remaining}")
    finally:
        f.close()


if __name__ == "__main__":
    main()
