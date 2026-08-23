# -*- coding: utf-8 -*-
"""
A2 参考文献重编号脚本 (2026-08-19)
==================================
背景: v0.7 中 Ref 22 (Ahluwalia 肾癌ECM) 与 Ref 23 (Song 脑静脉血栓)
与"ZP3 在肝癌/胰腺癌的预后意义"主张无关 (真正的支撑是 Ref 20=Bi HCC, 21=Lyu 胰腺)。
做法:
  1. 将引言句 "...prognostic significance in hepatocellular and pancreatic cancer [22,23]"
     改为引用 [20,21] (两者支持该主张, 属重复引用, 合法)。
  2. 删除 Ref 22, 23 (变为孤引用)。
  3. 全篇把编号 24-33 重排为 22-31 (映射见下)。
  4. 双向校验: 每个正文 [n] 均有对应参考文献条目, 且 1..31 连续、无遗漏、无孤儿。
输出: article2/manuscripts/Article2_A2_Cancers格式_v0.8.md
"""
import re
import sys

SRC = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\manuscripts\Article2_A2_Cancers格式_v0.7.md"
DST = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\manuscripts\Article2_A2_Cancers格式_v0.8.md"

# 重编号映射 (旧 -> 新); 22,23 被删除
MAP = {i: i for i in range(1, 22)}          # 1-21 -> 1-21
MAP.update({24: 22, 25: 23, 26: 24, 27: 25,
            28: 26, 29: 27, 30: 28, 31: 29,
            32: 30, 33: 31})                 # 24-33 -> 22-31

with open(SRC, encoding="utf-8") as f:
    text = f.read()

# ---------- 1. 先把 [22,23] 那句改为 [20,21] ----------
old_sentence = "carries prognostic significance in hepatocellular and pancreatic cancer [22,23]"
new_sentence = "carries prognostic significance in hepatocellular and pancreatic cancer [20,21]"
assert old_sentence in text, "intro sentence not found"
text = text.replace(old_sentence, new_sentence)

# ---------- 2. 处理正文中的引用编号 ----------
# 正文 = References 标题之前的部分
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
# 处理含区间的 [a-b] 一并纳入上面的正则会漏 (上面 pattern 不含 '-').
# 故对仍含 '-' 的区间单独处理:
def replace_range(m):
    return "[" + ",".join(str(MAP[int(n)]) for n in m.group(1).split("-")) + "]"
new_body = re.sub(r"\[(\d+-\d+)\]", replace_range, new_body)

# ---------- 3. 删除 Ref 22,23 条目并重排 References 标签 ----------
ref_lines = refs.split("\n")
kept = []
for ln in ref_lines:
    m = re.match(r"^(\d{2})\.\s", ln)
    if m:
        num = int(m.group(1))
        if num in (22, 23):
            continue  # 删除
        kept.append((num, ln))
    else:
        # 非条目行 (如标题、空行) 原样保留, 但需去重标题行
        kept.append((None, ln))

# 重建 References: 按顺序重排编号标签
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

# ---------- 4. 校验双向一致性 ----------
# 提取正文所有引用编号
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
# 提取参考文献条目编号
ref_nums = set()
for ln in new_refs.split("\n"):
    m = re.match(r"^(\d{2})\.\s", ln)
    if m:
        ref_nums.add(int(m.group(1)))
# 正文开头的标题/摘要区也含引用, 一并纳入 body_nums
all_body = new_body
# 摘要区在 Simple Summary / Abstract, 已在 new_body 内。

missing_refs = sorted(ref_nums - body_nums)   # 有条目但正文未引用 = 孤儿
missing_body = sorted(body_nums - ref_nums)   # 正文引用但无条目 = 悬空
expected = set(range(1, 32))

print(f"正文引用编号集合: {sorted(body_nums)}")
print(f"参考文献编号集合: {sorted(ref_nums)}")
print(f"孤儿引用(条目存在未引用): {missing_refs}")
print(f"悬空引用(引用无条目): {missing_body}")
print(f"编号 1..31 完整: {body_nums == expected == ref_nums}")

if missing_refs or missing_body or body_nums != expected or ref_nums != expected:
    print("!! 校验失败, 不写出 v0.8")
    sys.exit(1)

with open(DST, "w", encoding="utf-8") as f:
    f.write(new_body + "\n" + new_refs)
print(f"OK -> 已写出 {DST}")
print(f"正文长度: {len(new_body)}, 引用区长度: {len(new_refs)}")
