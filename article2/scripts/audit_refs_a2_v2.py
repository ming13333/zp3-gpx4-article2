# -*- coding: utf-8 -*-
"""
A2 参考文献只读审计 v2 (2026-08-19)
====================================
修复 v1 的 en-dash (–) 问题: 稿件引用区间用 en-dash (U+2013) 而非 ASCII '-',
导致 v1 误报大量 '孤儿引用'。本版把 ASCII '-' 与 en-dash '–' 均视为区间分隔符。
只读, 不修改任何文件。
"""
import re

SRC = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\manuscripts\Article2_A2_Cancers格式_v0.7.md"
with open(SRC, encoding="utf-8") as f:
    text = f.read()

marker = "## References"
ridx = text.find(marker)
body = text[:ridx]
refs = text[ridx:]

# 参考文献区编号
ref_nums = []
for ln in refs.split("\n"):
    m = re.match(r"^(\d{2})\.\s", ln.strip())
    if m:
        ref_nums.append(int(m.group(1)))

# 正文引用: 只按逗号/空白分隔, 每个 token 识别为 单编号 或 区间(- 或 –)
def tokenize(inner):
    out = []
    for tok in re.split(r"[,;，\s]+", inner.strip()):
        tok = tok.strip()
        if tok:
            out.append(tok)
    return out

def expand(tok):
    """返回一个 token 展开的编号集合。tok 形如 '5', '5-9', '5–9'."""
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
print("参考文献区编号数:", len(ref_nums))
print(f"连续 {ref_nums[0]}..{ref_nums[-1]}: {ref_nums == list(range(1,34))}")
print()
print("悬空引用 (正文引用但无条目):", sorted(cited - set(ref_nums)) or "无")
print("孤儿条目 (有条目但正文未引用):", sorted(set(ref_nums) - cited) or "无")
print()
print("逐条引用位置:")
for n in range(1, 34):
    occ = occur.get(n, [])
    line_of = ",".join(str(x) for x in occ[:5])
    print(f"  [{n:2d}] 引用{len(occ):2d}处 (行 {line_of})")
