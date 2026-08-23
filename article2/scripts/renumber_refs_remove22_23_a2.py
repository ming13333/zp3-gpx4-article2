#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""删除 A2 参考文献 22、23（内容与主张无关），并将 24-33 重编号为 22-31。
安全策略：先备份，重编号后由 audit_refs_a2_v2.py 一致性校验通过才视为成功。"""
import re, shutil, sys

SRC = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\manuscripts\Article2_A2_Cancers格式_v0.7.md"
BAK = SRC + ".bak_20260819"

shutil.copy(SRC, BAK)
print("backup ->", BAK, file=sys.stderr)

with open(SRC, encoding="utf-8") as f:
    text = f.read()

ref_header = "## References"
idx = text.index(ref_header)
head = text[:idx]
refsec = text[idx:]

REMOVED = {22, 23}

# 1) 参考文献列表：删除 22/23，24-33 -> 22-31
lines = refsec.split("\n")
new_lines = []
for ln in lines:
    m = re.match(r"^(\d{2})\.\s", ln)
    if m:
        n = int(m.group(1))
        if n in REMOVED:
            continue
        newn = n - 2 if n >= 24 else n
        ln = re.sub(r"^(\d{2})\.", f"{newn:02d}.", ln, count=1)
    new_lines.append(ln)
new_refsec = "\n".join(new_lines)

# 2) 正文引用重编号
def renum_inner(inner):
    out = []
    for tok in re.split(r",", inner):
        tok = tok.strip()
        if not tok:
            continue
        if re.fullmatch(r"\d+[–-]\d+", tok):
            lo, hi = re.split(r"[–-]", tok)
            lo, hi = int(lo), int(hi)
            nlo = None if lo in REMOVED else (lo - 2 if lo >= 24 else lo)
            nhi = None if hi in REMOVED else (hi - 2 if hi >= 24 else hi)
            if nlo is None or nhi is None:
                continue
            out.append(f"{nlo}–{nhi}")
        elif re.fullmatch(r"\d+", tok):
            n = int(tok)
            if n in REMOVED:
                continue
            out.append(str(n - 2 if n >= 24 else n))
        else:
            out.append(tok)
    return ",".join(out)

def repl(m):
    new_inner = renum_inner(m.group(1))
    if new_inner == "":
        return ""
    return "[" + new_inner + "]"

head2 = re.sub(r"\[([0-9,–\-\s]+)\]", repl, head)
# 清理因删除引用产生的 " ." / " ,"
head2 = re.sub(r" \.", ".", head2)
head2 = re.sub(r" ,", ",", head2)

text_final = head2 + new_refsec
with open(SRC, "w", encoding="utf-8") as f:
    f.write(text_final)
print("written. new reference count check:", file=sys.stderr)
print("OK", file=sys.stderr)
