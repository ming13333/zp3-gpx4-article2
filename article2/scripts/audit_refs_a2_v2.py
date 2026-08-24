# -*- coding: utf-8 -*-
"""
A2 references read-only audit v2 (2026-08-19)
====================================
Fix v1's en-dash (–) issue: manuscript citation ranges use en-dash (U+2013) rather than ASCII '-',
This caused v1 to falsely report many 'orphan citations'. This version treats both ASCII '-' and en-dash '–' as range separators.
Read-only, does not modify any files.
"""
import re

SRC = r"C:\D\workbuddy\Research\ExtracellularGPX4Immunosuppression\article2\manuscripts\Article2_A2_CancersFormat_v0.7.md"
with open(SRC, encoding="utf-8") as f:
    text = f.read()

marker = "## References"
ridx = text.find(marker)
body = text[:ridx]
refs = text[ridx:]

# Reference section numbering
ref_nums = []
for ln in refs.split("\n"):
    m = re.match(r"^(\d{2})\.\s", ln.strip())
    if m:
        ref_nums.append(int(m.group(1)))

# In-text citations: split only by comma/whitespace; each token is recognized as a single number or a range (- or –)
def tokenize(inner):
    out = []
    for tok in re.split(r"[,;，\s]+", inner.strip()):
        tok = tok.strip()
        if tok:
            out.append(tok)
    return out

def expand(tok):
    """Return the set of numbers expanded from a token. tok looks like '5', '5-9', '5–9'."""
    if re.fullmatch(r"\d+", tok):
        return {int(tok)}
    if re.fullmatch(r"\d+[–-]\d+", tok):
        lo, hi = re.split(r"[–-]", tok)
        return set(range(int(lo), int(hi) + 1))
    return set()

occur = {}
for i, line in enumerate(body.split("\n"), 1):
    for m in re.finditer(r"\[([0-9,\-–\s]+)\]", line):
        inner = m.group(1)
        for tok in tokenize(inner):
            for n in expand(tok):
                occur.setdefault(n, []).append(i)

cited = set(occur)
print("Reference section entry count:", len(ref_nums))
print(f"Continuous {ref_nums[0]}..{ref_nums[-1]}: {ref_nums == list(range(1,34))}")
print()
print("Dangling references (cited in text but no entry):", sorted(cited - set(ref_nums)) or "none")
print("Orphan entries (have entry but not cited in text):", sorted(set(ref_nums) - cited) or "none")
print()
print("Per-reference citation locations:")
for n in range(1, 34):
    occ = occur.get(n, [])
    line_of = ",".join(str(x) for x in occ[:5])
    print(f"  [{n:2d}] referenced {len(occ):2d} times (line {line_of})")
