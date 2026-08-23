# -*- coding: utf-8 -*-
"""
A2 参考文献纯只读审计 (2026-08-19)
==================================
目标: 不修改文件, 只输出当前稿件 (v0.7, 正文字段已改但引用区同原 v0.7) 的
  引用完整性审计 + 每条引用的正文引用位置, 供人工/后续按"替换不匹配条目"决策。
不重排、不删除。逐条列出:
  1) 参考文献区每条编号是否存在
  2) 每条是否被正文引用 (引用位置列表)
  3) 正文引用编号是否有悬空 (引用无条目)
  4) 顺序是否连续
"""
import re

SRC = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\manuscripts\Article2_A2_Cancers格式_v0.7.md"
with open(SRC, encoding="utf-8") as f:
    text = f.read()

marker = "## References"
ridx = text.find(marker)
body = text[:ridx]
refs = text[ridx:]

# ---- 1. 解析参考文献区条目编号 ----
ref_nums, ref_lines = [], {}
for ln in refs.split("\n"):
    m = re.match(r"^(\d{2})\.\s", ln.strip())
    if m:
        n = int(m.group(1))
        ref_nums.append(n)
        ref_lines[n] = ln.strip()

# ---- 2. 扫描正文中所有引用 token, 记录每个编号出现位置(所在行) ----
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
# ---- 3. 审计 ----
print("=" * 70)
print("参考文献区编号 (应连续 1..33):")
print(sorted(ref_nums))
seq_ok = ref_nums == list(range(34))
print(f"连续 1..33: {seq_ok}")
print()
print("=" * 70)
print("悬空引用 (正文引用但参考文献区无此条目):")
dangling = sorted(cited - set(ref_nums))
print(sorted(cited - set(ref_nums)) if dangling else "  无")
print()
print("孤儿条目 (参考文献区有条目但正文从未引用):")
orphan = sorted(set(ref_nums) - cited)
print(sorted(orphan) if orphan else "  无")
print()
print("=" * 70)
print("每条引用在正文的引用次数 (行号):")
for n in range(1, 34):
    occ = occur.get(n, [])
    status = "CITED" if n in cited else "ORPHAN(未引用)"
    keep = "KEEP" if n not in (22, 23) else "?MISMATCH"
    flag = ""
    if n == 22: flag = "  <== 肾癌ECM, 疑不匹配'预后意义'"
    if n == 23: flag = "  <== 脑静脉血栓, 疑不匹配'预后意义'"
    print(f"  [{n:2d}] {status:16s} 引用{len(occ)}处  {flag}")
