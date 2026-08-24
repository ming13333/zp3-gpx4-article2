# -*- coding: utf-8 -*-
"""
A2 reference renumbering script (2026-08-19)
==================================
Background: In v0.7, Ref 22 (Ahluwalia renal cancer ECM) and Ref 23 (Song cerebral venous thrombosis)
Unrelated to the claim of "prognostic significance of ZP3 in liver cancer/pancreatic cancer" (the real support is Ref 20=Bi HCC, 21=Lyu pancreatic).
Approach:
  1. Change the introduction sentence "...prognostic significance in hepatocellular and pancreatic cancer [22,23]"
     to cite [20,21] (both support the claim; duplicate citations are legitimate).
  2. Delete Ref 22, 23 (they become orphan citations).
  3. Renumber 24-33 throughout to 22-31 (mapping below).
  4. Bidirectional validation: every in-text [n] has a corresponding reference entry, and 1..31 are consecutive, with no omissions or orphans.
Output: article2/manuscripts/Article2_A2_CancersFormat_v0.8.md
"""
import re
import sys

SRC = r"C:\D\workbuddy\research\extracellular_GPX4_immunosuppression\article2\manuscripts\Article2_A2_CancersFormat_v0.7.md"
DST = r"C:\D\workbuddy\research\extracellular_GPX4_immunosuppression\article2\manuscripts\Article2_A2_CancersFormat_v0.8.md"

# Renumbering mapping (old -> new); 22,23 are deleted
MAP = {i: i for i in range(1, 22)}          # 1-21 -> 1-21
MAP.update({24: 22, 25: 23, 26: 24, 27: 25,
            28: 26, 29: 27, 30: 28, 31: 29,
            32: 30, 33: 31})                 # 24-33 -> 22-31

with open(SRC, encoding="utf-8") as f:
    text = f.read()

# ---------- 1. First change the [22,23] sentence to [20,21] ----------
old_sentence = "carries prognostic significance in hepatocellular and pancreatic cancer [22,23]"
new_sentence = "carries prognostic significance in hepatocellular and pancreatic cancer [20,21]"
assert old_sentence in text, "intro sentence not found"
text = text.replace(old_sentence, new_sentence)

# ---------- 2. Process citation numbers in the body ----------
# body = part before the References heading
marker = "## References"
ref_idx = text.find(marker)
assert ref_idx != -1, "References section marker not found"
body, refs = text[:ref_idx], text[ref_idx:]

def replace_citation(m):
    token = m.group(0)
    inner = token[1:-1]
    parts = re.split(r"\s*,\s*", inner)
    new_parts = []
    for p in parts:
        p = p.strip()
        if re.fullmatch(r"\d+", p):
            n = int(p)
            assert n in MAP, f"citation number {n} not in map"
            new_parts.append(str(MAP[n]))
        elif "-" in p and re.fullmatch(r"\d+-\d+", p):
            lo, hi = p.split("-")
            lo_i, hi_i = int(lo), int(hi)
            vals = [MAP[i] for i in range(lo_i, hi_i + 1)]
            assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)), f"range {p} not monotonic"
            new_parts.append(f"{vals[0]}-{vals[-1]}")
        else:
            raise ValueError(f"unparsed citation token: {p}")
    return "[" + ",".join(new_parts) + "]"

new_body = re.sub(r"\[\d+(?:\s*,\s*\d+)*\]", replace_citation, body)
# Handling ranges [a-b] together would be missed by the above regex (the above pattern doesn't contain '-').
# Therefore handle ranges still containing '-' separately:
def replace_range(m):
    return "[" + ",".join(str(MAP[int(n)]) for n in m.group(1).split("-")) + "]"
new_body = re.sub(r"\[(\d+-\d+)\]", replace_range, new_body)

# ---------- 3. Delete Ref 22,23 entries and renumber References labels ----------
ref_lines = refs.split("\n")
kept = []
for ln in ref_lines:
    m = re.match(r"^(\d{2})\.\s", ln)
    if m:
        num = int(m.group(1))
        if num in (22, 23):
            continue  # delete
        kept.append((num, ln))
    else:
        # Non-entry lines (e.g. headings, blank lines) are kept as-is, but heading lines need deduplication
        kept.append((None, ln))

# Rebuild References: renumber labels in order
out_lines = []
pending_num = 1
for orig_num, ln in kept:
    if orig_num is None:
        out_lines.append(ln)
    else:
        new_label = MAP[orig_num]
        new_ln = re.sub(r"^\d{2}\.\s", f"{new_label:02d}. ", ln, count=1)
        out_lines.append(new_ln)

new_refs = "\n".join(out_lines)

# ---------- 4. Verify bidirectional consistency ----------
# Extract all citation numbers from the body
body_nums = set()
for m in re.finditer(r"\[(\d+(?:-\d+)?(?:\s*,\s*\d+(?:-\d+)?)*)\]", new_body):
    inner = m.group(1)
    for tok in re.split(r"\s*,\s*", inner):
        tok = tok.strip()
        if "-" in tok and re.fullmatch(r"\d+-\d+", tok):
            lo, hi = tok.split("-")
            body_nums.update(range(int(lo), int(hi) + 1))
        elif re.fullmatch(r"\d+", tok):
            body_nums.add(int(tok))
# Extract reference entry numbers
ref_nums = set()
for ln in new_refs.split("\n"):
    m = re.match(r"^(\d{2})\.\s", ln)
    if m:
        ref_nums.add(int(m.group(1)))
# The title/abstract section at the beginning of the body also contains citations; include them in body_nums
all_body = new_body
# The abstract section is in Simple Summary / Abstract, already inside new_body.

missing_refs = sorted(ref_nums - body_nums)   # has entry but not cited in body = orphan
missing_body = sorted(body_nums - ref_nums)   # cited in body but no entry = dangling
expected = set(range(1, 32))

print(f"Body citation number set: {sorted(body_nums)}")
print(f"Reference number set: {sorted(ref_nums)}")
print(f"Orphan citations (entries exist but not cited): {missing_refs}")
print(f"Dangling citations (cited but no entry): {missing_body}")
print(f"Numbers 1..31 complete: {body_nums == expected == ref_nums}")

if missing_refs or missing_body or body_nums != expected or ref_nums != expected:
    print("!! Validation failed, not writing v0.8")
    sys.exit(1)

with open(DST, "w", encoding="utf-8") as f:
    f.write(new_body + "\n" + new_refs)
print(f"OK -> wrote {DST}")
print(f"Body length: {len(new_body)}, reference section length: {len(new_refs)}")
