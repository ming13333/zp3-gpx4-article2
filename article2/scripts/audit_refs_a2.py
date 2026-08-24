# -*- coding: utf-8 -*-
"""
A2 reference read-only audit (2026-08-19)
==================================
Goal: do not modify the file; only output the current manuscript (v0.7, body text has been changed but the reference section is the same as original v0.7)'s
  citation integrity audit + the in-text citation location of each reference, for manual/follow-up decision to "replace mismatched entries".
Do not reorder or delete. List each entry:
  1) whether each number in the reference section exists
  2) whether each entry is cited in the body (list of citation locations)
  3) whether body citation numbers are dangling (citation without an entry)
  4) whether the order is consecutive
"""
import re

SRC = r"C:\D\workbuddy\research\extracellularGPX4immunosuppression\article2\manuscripts\Article2_A2_Cancersformat_v0.7.md"
with open(SRC, encoding="utf-8") as f:
    text = f.read()

marker = "## References"
ridx = text.find(marker)
body = text[:ridx]
refs = text[ridx:]

# ---- 1. Parse reference section entry numbers ----
ref_nums, ref_lines = [], {}
for ln in refs.split("\n"):
    m = re.match(r"^(\d{2})\.\s", ln.strip())
    if m:
        n = int(m.group(1))
        ref_nums.append(n)
        ref_lines[n] = ln.strip()

# ---- 2. Scan all citation tokens in the body, record the occurrence position (line) of each number ----
occur = {}   # n -> list of lines
for i, line in enumerate(body.split("\n"), 1):
    for m in re.finditer(r"\[([0-9,\-\s]+)\]", line):
        inner = m.group(1)
        for tok in re.split(r"\s*,\s*", inner.strip()):
            tok = tok.strip()
            if re.fullmatch(r"\d+", tok):
                n = int(tok); occur.setdefault(n, []).append(i)
            elif re.fullmatch(r"\d+-\d+", tok):
                lo, hi = tok.split("-")
                for n in range(int(lo), int(hi) + 1):
                    occur.setdefault(n, []).append(i)

cited = set(occur)
# ---- 3. Audit ----
print("=" * 70)
print("Reference list numbering (should be consecutive 1..33):")
print(sorted(ref_nums))
seq_ok = ref_nums == list(range(34))
print(f"Consecutive 1..33: {seq_ok}")
print()
print("=" * 70)
print("Dangling references (cited in text but missing from reference list):")
dangling = sorted(cited - set(ref_nums))
print(sorted(cited - set(ref_nums)) if dangling else "  none")
print()
print("Orphan entries (in reference list but never cited in text):")
orphan = sorted(set(ref_nums) - cited)
print(sorted(orphan) if orphan else "  none")
print()
print("=" * 70)
print("Occurrences of each reference in text (line numbers):")
for n in range(1, 34):
    occ = occur.get(n, [])
    status = "CITED" if n in cited else "ORPHAN(not cited)"
    keep = "KEEP" if n not in (22, 23) else "?MISMATCH"
    flag = ""
    if n == 22: flag = "  <== renal cancer ECM, suspected mismatch with 'prognostic significance'"
    if n == 23: flag = "  <== cerebral venous thrombosis, suspected mismatch with 'prognostic significance'"
    print(f"  [{n:2d}] {status:16s} cited at {len(occ)} places  {flag}")
