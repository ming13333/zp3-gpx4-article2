# A2 降 AI 痕迹反思报告（Article2_A2_Cancers格式_v0.4）

> 基于 academic-paper-humanizer（Academic 模式）四维度检测 + 全文机器扫描（2026-08-17）。

---

## 0. 合规前提（先于一切）

降 AI 痕迹的目标是让文本**更自然、更准确、更可验证**，像成熟的作者写作，而**不是**隐藏 AI 使用或欺骗检测器。MDPI 要求披露生成式 AI 的用途，A2 的 Acknowledgments 已含 AI 辅助语言编辑声明，与期刊政策一致，**务必保留**。全文用 iThenticate 检查的是抄袭，不是 AI；审稿人识破"AI 味"靠的是模板感与不可验证的引用，这两点正是本文要解决的。

---

## 1. 四维度诊断结果

### 1.1 语言层 —— 整体干净 ✅（1 处小问题）

机器扫描结果（正文 1–135 行）：

- Moreover / Furthermore / Notably / Importantly / Overall / Additionally 等模板连接词：**0 处**
- demonstrate / reveal / highlight / underscore / pivotal / crucial 等"学术化动词"滥用：**0 处**（参考文献标题中的 reveals 为原文，不算）
- **唯一问题**：Simple Summary（L13）有 2 处 em-dash："one candidate gene—ZP3—can serve as such a predictor"。AI 偏好用破折号做插入语，改为逗号或括号即可。

### 1.2 逻辑层 —— 强项 ✅

- Introduction：背景 → 明确 Gap（单基因 ICB 标志物跨队列不可重复、ZP3 无多队列评估）→ 预设问题，链条完整，无"大背景→薄 Gap"跳步。
- Discussion：六点 limitations + 明确边界声明（"cannot directly validate or refute the proposed molecular axis"、"association-level tests cannot establish causality"），**无"结论强于数据"**。这是本文最不像 AI 的部分，保持。

### 1.3 学术细节层 —— 强项，但有一个致命风险点 ⚠️

**强项**：数字密度极高且具体（n=27/298/105、P=0.90/0.86/0.0023、r 与 95% CI、中位数 3.78 vs 3.35、MDE |r|=0.22、TOST P=0.045、Entrez ID 8277、FMOne mutation burden）——这是最强的"人类感"来源，保持。

**致命风险：参考文献可验证性**。75 条中 **2026 年条目约 30+ 条**，含大量"看起来合理但难以核实"的条目（如 Ref 45 "Supty S, et al. ...Pathol. Res. Pract. 2026, 286, 156634"、Ref 52 "Zhang X, et al. ...Front. Immunol. 2026, 17, 1404802"、Ref 69 "Xu X, et al. ...J. Craniofac. Surg. 2026"）。这是全文最大的 AI 痕迹雷区：**审稿人抽查引用真伪，一条虚构 DOI 即可毁掉可信度**。

### 1.4 风格层 —— 最明显的 AI 痕迹 ❌

**Results 3.1 / 3.2 / 3.3 三小节句式几乎同构**，同一模板套三次：

> In [队列] (n=..., [癌种] patients treated with [药物]; cohort reported in [ref]), ZP3 expression did/did not differ between responders (n=...) and non-responders (n=...) (Mann–Whitney U test, P=...; rank-biserial r=..., 95% CI ...; median ...)。

**图注 Figure 1–3 同样同构**（"Box plots of ZP3 expression in responders (CR/PR, n=X) versus non-responders (PD, n=Y) in [cohort]..."）。这是"模板化"的最典型信号。

---

## 2. 高优先级修改清单

| # | 位置 | 问题 | 动作 |
|---|---|---|---|
| 1 | Results 3.1–3.3 | 三节同构句式 | **重写打散**（见 §4 示例）：三段分别从"效应量视角 / 等效检验视角 / 方向与稳健性预告视角"切入 |
| 2 | 参考文献（75 条） | ~30+ 条 2026 年条目未逐条核验 | **逐条 Crossref/PubMed 核验**（批量 API 或 nature-ref-verifier）；核验通过的补 DOI，无法核实的删除或替换为可核验文献 |
| 3 | Simple Summary L13 | 2 处 em-dash | 改为逗号/括号 |
| 4 | Figure captions 1–3 | 图注同构 | 措辞差异化（Fig1/Fig2 微调，Fig3 已略异） |
| 5 | Discussion Limitations | "Several limitations apply. First...Sixth..." 机械分条 | 合并为自然衔接的段落（见 §4 示例），保留全部信息 |

## 3. 低优先级（可不动）

- Abstract 的结构化分节（Background/Objectives → ...）是 Cancers 要求，保留。
- "This cohort is larger than GSE78220 and provides a more precise estimate"（3.2 末句）——自然，保留。
- Author Contributions 中 "writing—original draft preparation" 是 CRediT 标准术语（MDPI 模板原文），不是 AI 痕迹。
- 6 处参考文献的 "reveals" 是文献原标题，不动。

## 4. 改写示例（可直接采用）

### 4.1 Results 3.1（原 82 词，从"效应量"切入）

**原文：**
> In GSE78220 (n=27 pre-treatment melanoma patients treated with anti-PD-1 pembrolizumab; cohort reported in [27]), ZP3 expression at baseline did not differ between responders (CR/PR, n=15) and non-responders (PD, n=12) (Mann–Whitney U test, P=0.90; rank-biserial r=0.03, 95% CI −0.41 to 0.48; median 3.78 in responders versus 3.35 in non-responders). Mean ZP3 was 6.14 (CR), 4.70 (PR), and 5.56 (PD), and one-way ANOVA across the three groups was not statistically significant (P=0.86). The point estimate of the effect size is close to zero with a wide confidence interval that does not exclude moderate effects in either direction, indicating a null result rather than a precise demonstration of equivalence in this small cohort.

**改写：**
> In GSE78220 (anti-PD-1 pembrolizumab; n=27 pre-treatment samples; Hugo et al. [27]), the ZP3 distributions of the 15 responders (CR/PR) and 12 non-responders (PD) overlapped almost completely: median 3.78 versus 3.35, Mann–Whitney P=0.90, rank-biserial r=0.03 (95% CI −0.41 to 0.48). Means moved little across the three response categories (6.14 CR, 4.70 PR, 5.56 PD) and one-way ANOVA gave P=0.86. With r near zero but a wide interval, this cohort can neither confirm nor exclude a moderate association; the result is a null, not a demonstrated equivalence.

### 4.2 Results 3.2（原 140 词，从"等效检验"切入，压缩）

**改写：**
> IMvigor210 provided the largest test set: 298 evaluable urothelial carcinoma patients treated with atezolizumab (68 responders, 230 non-responders; 50 of the 348 profiled patients lacked a response annotation and were excluded; cohort as described in [28–30]). ZP3 did not separate the two groups (P=0.86; median 441.0 versus 460.5), and the estimate, r=0.015 (95% CI −0.14 to 0.17), is the tightest of the three cohorts. Two further analyses frame this null: a TOST equivalence test with a small-effect bound of |r|=0.15 was positive (P=0.045), and the minimum detectable effect at 80% power was |r|=0.22. A small-to-moderate association, had it existed, would therefore likely have been detected, making the negative result informative rather than underpowered.

### 4.3 Results 3.3（原 54 词，从"唯一阳性 + 主动预告"切入）

**改写：**
> GSE91061 was the exception. Among 105 evaluable melanoma patients treated with nivolumab [31], non-responders carried markedly higher ZP3 than responders (mean 2.17 versus 0.33; Mann–Whitney P=0.0023; r=0.42, 95% CI 0.15 to 0.69). Because this was the only nominally significant association across the three cohorts, we examined its robustness before interpreting it (Section 3.5).

### 4.4 Discussion Limitations（原六点式，改为自然衔接）

**改写：**
> Several limitations bound these conclusions. We evaluated bulk-tumor ZP3 transcript abundance, not cell-surface ZP3 protein, receptor occupancy, or the proposed GPX4–ZP3 interaction; an association-level design cannot validate or refute that molecular axis. GSE78220 is small (n=27), and its wide confidence interval leaves moderate effects possible in either direction. GSE91061 lacked complete BRAF/NRAS and TMB annotations in the molecular-subtype file, and prior treatment lines were unavailable, so full multivariable adjustment was not possible [8,57]. Molecular subtypes were expression-inferred rather than TCGA-official, and subtype-stratified estimates rest on sparse subgroups. Positive controls reached significance only in the largest cohort (IMvigor210), so pipeline sensitivity in the smaller cohorts is not guaranteed. Finally, the public-data, observational design supports association-level statements only; no causal claim is made.

（信息零丢失：原文六点全部保留。）

---

## 5. 合规提醒（最终检查）

1. **引用核验是唯一的"造假级"风险**：75 条须逐条可查。优先用 Crossref API 批量核验 DOI/卷期/页码，2026 年条目重点核查；核验记录建议留档（投稿包内附）。
2. AI 披露语句保持现状（Acknowledgments），与 MDPI 政策一致；不要因降痕而删除披露。
3. 所有数字、统计量、队列名、基因名在改写中保持字节级一致（示例已保证）。
4. 改写后建议重跑本报告的 §1.1 扫描（连接词/动词/破折号计数归零）作为验收。

---

## 附：验收标准（改写完成后）

- [ ] em-dash 计数 = 0（Author Contributions 的 CRediT 术语除外）
- [ ] Moreover/Furthermore/Notably/Overall 等模板连接词 = 0
- [ ] 3.1/3.2/3.3 首句句式互不相同
- [ ] 75 条参考文献 100% 可通过 Crossref/PubMed 检索到 DOI
- [ ] 摘要 193 词、全部 P 值/CI/样本量与 v0.4 一致
