# A2 审稿质疑整改 Checklist（资深审稿人模拟审查 → 逐条整改）

**日期：** 2026-08-17
**稿件：** `article2/manuscripts/Article2_A2_Cancers格式_v0.7.md`（基于 v0.6 修订）
**状态：** 全部 17 条质疑已逐条整改并验证 ✅

---

## 一、科学性质疑（5 条）

| # | 质疑 | 整改动作 | 位置 | 状态 |
|---|---|---|---|---|
| 1 | "预测性生物标志物"结论证据不足（无 untreated control） | 全文将 predictive biomarker 收敛为 response association；结论明确"单臂队列无法做 treatment-interaction 预测推断"；Simple Summary/Abstract/Conclusions/Discussion 全部同步 | 标题、Simple Summary、Abstract、§1、§4、§5、cover letter | ✅ |
| 2 | 三队列不同质（癌种/药物/平台/反应定义不同） | Discussion 增加显式声明：队列异质性支持"无泛化性 ZP3–response 关联"的结论，但无法排除未采样情境的 context-specific 关联 | §4 Discussion 首段 | ✅ |
| 3 | IMvigor210 排除 50 例可能造成选择偏倚 | Methods §2.1 补充：50 例因缺二元反应标注而排除，缺失由数据可得性驱动而非反应相关标准，但无法完全排除非随机缺失 | §2.1 | ✅ |
| 4 | GSE91061 on-treatment 敏感性分析伪重复（105 样本≠105 患者） | 实证统计：105 样本=62 唯一患者，42 患者双 visit（脚本验证）；Methods/Results/图注全部更新为患者级口径；新增 patient-level 聚合敏感性（n=62, P=0.18） | §2.1、§2.4、§3.3、§3.5、Fig5、Supp Table S5 | ✅ |
| — | 代码级 bug：gse91061_zp3_analysis_v2.py 缺 `import os` | 修复并实跑验证（Pre n=49, r=−0.195 与冻结表一致） | `article2/scripts/gse91061_zp3_analysis_v2.py` | ✅ |

## 二、统计方法质疑（5 条）

| # | 质疑 | 整改动作 | 位置 | 状态 |
|---|---|---|---|---|
| 5 | rank-biserial CI 近似法需稳健性论证 | 新增 bootstrap 95% CI（B=20000，分层抽样）：GSE78220 −0.48~0.43 / IMvigor210 −0.18~0.15 / GSE91061 −0.58~0.22，与近似 CI 高度一致；Methods/Results/图注/Supp Table S4 更新 | §2.3、§3 Table 1 下、Fig1-3 图注、`a2_robustness_bootstrap.csv` | ✅ |
| 6 | TOST P=0.045 边界性，不宜过度解读 | 措辞改为 "borderline evidence compatible with equivalence within the pre-specified small-effect bound"，不再声称 "supported equivalence" | §2.6、§3.2、cover letter | ✅ |
| 7 | MDE/post-hoc power 不宜作为充分 power 证据 | 明确 MDE 为 interpretive context；承认 GSE78220/GSE91061 无法排除中等效应 | §2.6、§3.2 | ✅ |
| 8 | 三个主要检验未做多重校正 | 补充说明：以 cohort-level replication 为目标而非选择单一显著关联，故不校正；探索性 P 值描述性解读 | §2.3 | ✅ |
| 9 | GSE91061 仅 10 事件，logistic 不稳定 | 新增 Firth 惩罚 logistic：univariable OR=0.82 (0.38–1.77) P=0.61；+IFN-γ P=0.83；+Supp P=0.67——全部不显著，与标准 logistic 一致；Methods/Results/Supp Table S3 更新 | §2.4、§3.5、`a2_gse91061_firth.csv` | ✅ |

## 三、ZP3 测量与生物学解释（3 条）

| # | 质疑 | 整改动作 | 位置 | 状态 |
|---|---|---|---|---|
| 10 | bulk RNA 不代表功能性蛋白/受体 | 标题加 "Bulk Transcript Abundance"；Abstract/§1/§4 明确 transcript 层面结论不外推至 GPX4–ZP3 轴 | 标题、Abstract、§4 limitation 1 | ✅ |
| 11 | GSE91061 Entrez ID 修复需可复核证据 | Methods §2.2 补充：Entrez 7784=ZP3 经 NCBI Gene 与平台注释核验；8277=SP5 未使用；脚本含行匹配断言 | §2.2 | ✅ |
| 12 | 三队列表达值不可直接比较 | Table 1 下、Fig1-3 图注全部加 "cohort-specific scale, not cross-cohort comparable" | Table 1、Fig1-3 图注 | ✅ |

## 四、阳性对照（2 条）

| # | 质疑 | 整改动作 | 位置 | 状态 |
|---|---|---|---|---|
| 13 | "pre-specified positive controls" 需时间证据 | 措辞统一为 "specified in the analysis plan（仓库内版本化）"；Results §3.4 保留预指定表述 | §2.5、§3.4 | ✅ |
| 14 | 阳性对照不能证明 ZP3 分析"敏感" | Discussion limitation 6 扩展：TMB/IHC 为不同测量模态、效应量更大，仅证明"数据集可复现已知关联"，不能保证每个阴性结果都来自充足 power | §4 limitation 6 | ✅ |

## 五、定位与引用（3 条）

| # | 质疑 | 整改动作 | 位置 | 状态 |
|---|---|---|---|---|
| 15 | 标题略强 | 改为 "ZP3 Expression Shows No Consistent Association with Immune Checkpoint Blockade Response in Three Public Cohorts: An Association-Level Evaluation of Bulk Transcript Abundance"；cover letter 同步 | 标题、cover letter | ✅ |
| 16 | Discussion 超推论措辞 | "signal was driven by..." → "signal was observed only in analyses that included on-treatment measurements"；"reflecting immune-context features" → "consistent with, but not sufficient to establish..." | §3.5、§4 | ✅ |
| 17 | 参考文献背景过宽（34–75 大量 glioma/ferroptosis 泛引） | 删除 References 34–75 共 42 条与 A2 主题无关文献，保留 1–33；正文引用全部指向 1–33 且一一对应（脚本校验通过）；Discussion 中 [34-36]/[37-40]/[41-51]/[52-56]/[58-69]/[70-74]/[75] 相应改写或删除 | References、§4 | ✅ |

---

## 新增产物（已登记 verify_manifest）

| 文件 | 说明 |
|---|---|
| `article2/scripts/robustness_ci_a2.py` | bootstrap CI + Firth logistic + patient-level 聚合（A2 登记） |
| `article2/results/a2_robustness_bootstrap.csv` | 三队列 bootstrap 95% CI |
| `article2/results/a2_gse91061_firth.csv` | Firth 惩罚 logistic 结果（3 模型） |
| `article2/results/a2_gse91061_patient_level.csv` | GSE91061 患者级敏感性（n=62, P=0.1753） |
| `article2/figures/A2_Fig5.svg` | 左面板新增 Firth×3 + Patient-level 共 4 根柱 |

## 验证记录

- [x] `gse91061_zp3_analysis_v2.py` 修复 import os 后实跑：Pre n=49, r=−0.195（与冻结表一致）
- [x] `robustness_ci_a2.py` 实跑：bootstrap CI 与近似 CI 一致；Firth 全不显著；patient-level P=0.175
- [x] 17 张 SVG 全量再生成，XML 0 invalid；A2_Fig5 含新内容（Firth/Patient-level/62 patients/0.18）
- [x] Abstract 198 词（Cancers ≤200）
- [x] References 33 条，正文引用编号 1–33 全部有对应条目（脚本校验 0 缺失）
- [x] `verify_manifest.py --repo .`：23 项通过 0 问题
- [x] cover letter 标题与数字同步（P=0.45–0.58 standard; Firth P=0.61–0.83; patient-level P=0.18）
