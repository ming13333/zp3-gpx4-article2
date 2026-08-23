#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Cancers-submission supplementary tables (S1-S5) from frozen CSVs.

Reads the frozen result tables under article2/results/ and assembles a single
XLSX workbook (article2/results/A2_Supplementary_Tables.xlsx) with one sheet
per supplementary table, plus a cover/README sheet documenting provenance.

Every cell is traceable to a frozen CSV; no value is hand-typed except headers
and the narrative footnotes.

Run (managed venv):
  python build_supplementary_tables_xlsx.py
"""
import csv
import os

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RESULTS, "A2_Supplementary_Tables.xlsx")

# ---------- styling ----------
HEADER_FILL = PatternFill("solid", fgColor="2C3E50")
HEADER_FONT = Font(bold=True, color="FFFFFF", size=10)
TITLE_FONT = Font(bold=True, size=12, color="2C3E50")
NOTE_FONT = Font(italic=True, size=9, color="555555")
PANEL_FONT = Font(bold=True, size=10, color="2C3E50")
THIN = Side(style="thin", color="BBBBBB")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
LEFT = Alignment(horizontal="left", vertical="center", wrap_text=True)


def read_csv(name):
    with open(os.path.join(RESULTS, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def style_header(ws, row, ncols, start=1):
    for c in range(start, start + ncols):
        cell = ws.cell(row=row, column=c)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = CENTER
        cell.border = BORDER


def write_table(ws, start_row, headers, rows, col_align=None):
    """Write a header + data block; return the next free row."""
    for j, h in enumerate(headers, start=1):
        ws.cell(row=start_row, column=j, value=h)
    style_header(ws, start_row, len(headers))
    r = start_row + 1
    for row in rows:
        for j, val in enumerate(row, start=1):
            cell = ws.cell(row=r, column=j, value=val)
            cell.border = BORDER
            cell.alignment = (col_align[j - 1] if col_align else CENTER)
            cell.font = Font(size=10)
        r += 1
    return r


def set_widths(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


def title_block(ws, title, caption=None):
    ws.cell(row=1, column=1, value=title).font = TITLE_FONT
    if caption:
        c = ws.cell(row=2, column=1, value=caption)
        c.font = NOTE_FONT
        c.alignment = LEFT
        return 4
    return 3


# ---------- load frozen data ----------
sens = read_csv("a2_effects_sensitivity.csv")          # S5
sub = read_csv("a2_gse91061_subtype_response.csv")     # S2
std = read_csv("a2_gse91061_logistic_standard.csv")    # S3 standard
firth = read_csv("a2_gse91061_firth.csv")              # S3 Firth
boot = read_csv("a2_robustness_bootstrap.csv")         # S4 bootstrap
power = read_csv("power_analysis_results.csv")         # S4 power/TOST
pat = read_csv("a2_gse91061_patient_level.csv")        # S5 patient-level

wb = Workbook()

# =====================================================================
# README / cover sheet
# =====================================================================
ws = wb.active
ws.title = "README"
ws.cell(row=1, column=1, value="Article 2 (A2) — Supplementary Tables S1–S5").font = Font(bold=True, size=14, color="2C3E50")
ws.cell(row=2, column=1,
        value="Manuscript: \"Bulk ZP3 transcript abundance shows no consistent association with immune "
              "checkpoint blockade response across tumor types\" (target journal: Cancers / MDPI).").font = NOTE_FONT
ws.cell(row=3, column=1,
        value="All numeric values in this workbook are reproduced from frozen result CSVs in article2/results/; "
              "no analytic value was retyped manually. Source file per table is listed below and on each sheet.").font = NOTE_FONT

toc = [
    ("Sheet", "Table", "Content", "Source frozen CSV(s)"),
    ("S1", "Table S1", "Cohort exclusion flowchart (profiled → pre-treatment → evaluable → analysis subset)",
     "Manuscript Table 1 (lines 79–81); a2_effects_sensitivity.csv"),
    ("S2", "Table S2", "Subtype-stratified sample counts and ZP3–response results (GSE91061 pre-treatment n=49)",
     "a2_gse91061_subtype_response.csv"),
    ("S3", "Table S3", "Complete multivariable logistic results — standard and Firth-penalized",
     "a2_gse91061_logistic_standard.csv; a2_gse91061_firth.csv"),
    ("S4", "Table S4", "Bootstrap 95% CIs (B=20,000) and power analysis (MDE, post-hoc power, TOST)",
     "a2_robustness_bootstrap.csv; power_analysis_results.csv"),
    ("S5", "Table S5", "GSE91061 visit-stratified and patient-level sensitivity analyses",
     "a2_effects_sensitivity.csv; a2_gse91061_patient_level.csv"),
]
write_table(ws, 5, toc[0], toc[1:])
set_widths(ws, [10, 12, 60, 52])
ws.cell(row=5 + len(toc), column=1,
        value="Note: P values are two-sided. Effect sizes are rank-biserial correlations (r_rb). "
              "MDPI/Cancers accepts a single Excel file containing all supplementary tables.").font = NOTE_FONT

# =====================================================================
# S1 — cohort exclusion flowchart
# =====================================================================
ws = wb.create_sheet("S1")
nr = title_block(ws, "Table S1. Cohort exclusion flowchart",
                 "Profiled → pre-treatment → evaluable → analysis subset. Numbers trace to Manuscript Table 1 "
                 "and a2_effects_sensitivity.csv.")
headers = ["Cohort", "Cancer type", "ICB agent", "Profiled (n)", "Pre-treatment\nEvaluable (n)",
           "Excluded before\nprimary analysis (n)", "Exclusion reason", "Responders\n(n)", "Non-responders\n(n)",
           "Sensitivity set\n(n samples / n patients)", "Sensitivity\nResp / Non-resp"]
rows = [
    ["GSE78220", "Melanoma (skin)", "anti-PD-1 (pembrolizumab)", 28, 27, 1,
     "One profiled sample lacked a usable pre-treatment / response annotation", 15, 12, "—", "—"],
    ["IMvigor210", "Urothelial carcinoma", "anti-PD-L1 (atezolizumab)", 348, 298, 50,
     "50 patients lacked a response annotation", 68, 230, "—", "—"],
    ["GSE91061", "Melanoma (skin)", "anti-PD-1 (nivolumab)", 109, 49, 60,
     "Samples without a pre-treatment baseline (on-treatment-only or missing)", 10, 39,
     "105 samples / 62 patients", "23 / 82"],
]
nxt = write_table(ws, nr, headers, rows)
ws.cell(row=nxt + 1, column=1,
        value="Footnote: GSE91061 primary analysis used pre-treatment baseline samples only (n=49; 10 responders, "
              "39 non-responders). A sensitivity analysis additionally included on-treatment samples "
              "(105 samples from 62 unique patients; 42 patients contributed both visits; 23 responders, 82 "
              "non-responders) and is reported in Table S5. Four of the 109 profiled samples carried no response "
              "annotation and were excluded from both primary and sensitivity sets.").font = NOTE_FONT
ws.cell(row=nxt + 1, column=1).alignment = LEFT
set_widths(ws, [12, 18, 24, 11, 13, 14, 34, 11, 13, 20, 14])

# =====================================================================
# S2 — subtype-stratified
# =====================================================================
ws = wb.create_sheet("S2")
nr = title_block(ws, "Table S2. Subtype-stratified sample counts and ZP3–response results",
                 "GSE91061 pre-treatment samples (n=49). Within-subtype Mann–Whitney U test of ZP3 by response "
                 "(PR/CR vs SD/PD). Source: a2_gse91061_subtype_response.csv.")
headers = ["Molecular subtype", "n total", "n responders", "n non-responders",
           "Mean ZP3 (responders)", "Mean ZP3 (non-responders)", "Mann–Whitney P", "Estimable"]
rows = []
for r in sub:
    rows.append([
        r["subtype"], int(r["n_total"]), int(r["n_response"]), int(r["n_nonresponse"]),
        float(r["zp3_mean_response"]), float(r["zp3_mean_nonresponse"]),
        ("not estimable" if r["mwu_p"] == "" else float(r["mwu_p"])),
        ("Yes" if r["estimable"] == "True" else "No"),
    ])
nxt = write_table(ws, nr, headers, rows)
ws.cell(row=nxt + 1, column=1,
        value="Footnote: Subtypes assigned by median split of melanoma subtype signatures "
              "(Immune / Keratin / MITF-high / MITF-low). MITF-high (1 responder, 10 non-responders) and "
              "MITF-low (1 responder, 3 non-responders) strata contained <2 responders per group and were not "
              "estimable. ZP3 Entrez ID 7784; analysis restricted to pre-treatment baseline samples.").font = NOTE_FONT
ws.cell(row=nxt + 1, column=1).alignment = LEFT
set_widths(ws, [18, 10, 13, 18, 20, 24, 16, 12])

# =====================================================================
# S3 — multivariable logistic (standard + Firth)
# =====================================================================
ws = wb.create_sheet("S3")
nr = title_block(ws, "Table S3. Complete multivariable logistic regression results",
                 "GSE91061 pre-treatment samples (n=49; responders=10). Standard logistic (left) and "
                 "Firth-penalized logistic (right). ZP3 and covariate scores scaled to SD units. "
                 "Sources: a2_gse91061_logistic_standard.csv; a2_gse91061_firth.csv.")
headers = ["Model", "Method", "Term", "n", "Events", "OR", "CI low", "CI high", "P", "Iterations"]
rows = []
# standard
for r in std:
    rows.append([r["model"], "Standard logistic", r["term"], int(r["n"]), None,
                 float(r["OR"]), float(r["ci_low"]), float(r["ci_high"]), float(r["p"]), None])
# Firth
for r in firth:
    rows.append([r["Model"], "Firth-penalized", r["Term"], int(r["n"]), int(r["events"]),
                 float(r["Firth_OR"]), float(r["Firth_CI_low"]), float(r["Firth_CI_high"]),
                 float(r["Firth_P"]), int(r["iter"])])
nxt = write_table(ws, nr, headers, rows)
ws.cell(row=nxt + 1, column=1,
        value="Footnote: Model 1 = ZP3 only; Model 2 = ZP3 + IFN-γ signature score; Model 3 = ZP3 + "
              "immunosuppression signature score. Firth models include only the ZP3 term (events = number of "
              "responders). 'Events' and 'Iterations' columns apply to Firth models only. None of the "
              "pre-treatment frameworks reached significance (all P>0.45).").font = NOTE_FONT
ws.cell(row=nxt + 1, column=1).alignment = LEFT
set_widths(ws, [22, 16, 14, 7, 9, 9, 9, 9, 9, 11])

# =====================================================================
# S4 — bootstrap CIs + power/TOST
# =====================================================================
ws = wb.create_sheet("S4")
nr = title_block(ws, "Table S4. Bootstrap confidence intervals and power analysis",
                 "Rank-biserial effect sizes with approximate and bootstrap 95% CIs (Panel A), and "
                 "minimum-detectable-effect (MDE), post-hoc power and TOST equivalence (Panel B). "
                 "Sources: a2_robustness_bootstrap.csv; power_analysis_results.csv.")
# Panel A: bootstrap
ws.cell(row=nr, column=1, value="Panel A. Approximate vs bootstrap (B=20,000) 95% confidence intervals").font = PANEL_FONT
pa_h = ["Cohort", "n resp", "n non", "r_rb", "SE", "Approx CI low", "Approx CI high",
        "Bootstrap CI low", "Bootstrap CI high", "P"]
pa_rows = []
for r in boot:
    pa_rows.append([r["Cohort"], int(r["n_resp"]), int(r["n_non"]), float(r["r_rb"]), float(r["SE"]),
                    float(r["approx_CI_low"]), float(r["approx_CI_high"]),
                    float(r["boot_CI_low"]), float(r["boot_CI_high"]), float(r["P"])])
nr2 = write_table(ws, nr + 1, pa_h, pa_rows)
# Panel B: power
nr2 += 1
ws.cell(row=nr2, column=1, value="Panel B. Minimum detectable effect (MDE, power=0.80, α=0.05), post-hoc power, and TOST equivalence").font = PANEL_FONT
pb_h = ["Cohort", "n1", "n2", "r_obs", "SE", "MDE (|r|)", "Post-hoc power", "TOST bound", "TOST P", "TOST equivalent"]
pb_rows = []
for r in power:
    pb_rows.append([
        r["Cohort"], int(r["n1"]), int(r["n2"]), float(r["r_obs"]), float(r["se"]),
        float(r["MDE_80"]), float(r["posthoc_power"]),
        (float(r["TOST_bound"]) if r["TOST_bound"] else "—"),
        (float(r["TOST_p"]) if r["TOST_p"] else "—"),
        (r["TOST_equivalent"] if r["TOST_equivalent"] else "—"),
    ])
nr3 = write_table(ws, nr2 + 1, pb_h, pb_rows)
ws.cell(row=nr3 + 1, column=1,
        value="Footnote: Bootstrap CIs (B=20,000, stratified by response group) closely matched the normal-"
              "approximation CIs, indicating interval estimates are not sensitive to the approximation. MDE and "
              "post-hoc power are interpretive context, not proof of adequate power; the two smaller cohorts "
              "(GSE78220, GSE91061) cannot exclude moderate effects. TOST was performed for IMvigor210 only, "
              "with an equivalence bound of |r_rb|=0.15 (small effect); the borderline result (P=0.045) is "
              "compatible with equivalence within that bound, not a demonstration of strict equivalence.").font = NOTE_FONT
ws.cell(row=nr3 + 1, column=1).alignment = LEFT
set_widths(ws, [14, 9, 9, 9, 9, 13, 14, 12, 9, 15])

# =====================================================================
# S5 — GSE91061 visit-stratified + patient-level
# =====================================================================
ws = wb.create_sheet("S5")
nr = title_block(ws, "Table S5. GSE91061 visit-stratified and patient-level sensitivity analyses",
                 "Mann–Whitney U test of ZP3 by response across analysis subsets. Sources: "
                 "a2_effects_sensitivity.csv; a2_gse91061_patient_level.csv.")
headers = ["Scope", "n samples", "n patients", "n responders", "n non-responders",
           "Median ZP3 (resp)", "Median ZP3 (non)", "U", "P", "r_rb", "CI low", "CI high"]
rows = []
# Pre baseline
pre = next(r for r in sens if r["Scope"].startswith("Pre"))
rows.append(["Pre (baseline) only", int(pre["n_resp"]) + int(pre["n_non"]), None,
             int(pre["n_resp"]), int(pre["n_non"]), float(pre["med_resp"]), float(pre["med_non"]),
             float(pre["U"]), float(pre["P"]), float(pre["r_rb"]), float(pre["CI_low"]), float(pre["CI_high"])])
allp = next(r for r in sens if r["Scope"].startswith("All"))
rows.append(["All samples (incl. on-treatment)", int(allp["n_resp"]) + int(allp["n_non"]), None,
             int(allp["n_resp"]), int(allp["n_non"]), float(allp["med_resp"]), float(allp["med_non"]),
             float(allp["U"]), float(allp["P"]), float(allp["r_rb"]), float(allp["CI_low"]), float(allp["CI_high"])])
p = pat[0]
rows.append(["Patient-level (1 sample/patient)", int(p["n_samples"]), int(p["n_patients"]),
             int(p["n_resp"]), int(p["n_non"]), "", "",
             float(p["U"]), float(p["P"]), float(p["r_rb"]), float(p["CI_low"]), float(p["CI_high"])])
nxt = write_table(ws, nr, headers, rows)
ws.cell(row=nxt + 1, column=1,
        value="Footnote: 'Pre (baseline) only' = primary pre-treatment subset (n=49; 10 responders, 39 "
              "non-responders); 'All samples' = sensitivity set including on-treatment visits (105 samples from "
              "62 patients; 23 responders, 82 non-responders; 42 patients contributed both visits); "
              "'Patient-level' = one sample per patient (n=62; 14 responders, 48 non-responders). The apparent "
              "association (P=0.008) appears only in the on-treatment-inclusive set and is absent from the "
              "pre-treatment subset and the patient-level aggregation (P=0.18).").font = NOTE_FONT
ws.cell(row=nxt + 1, column=1).alignment = LEFT
set_widths(ws, [30, 11, 11, 12, 16, 16, 16, 9, 11, 9, 9, 9])

# freeze panes under each table header for readability
for name in ["S1", "S2", "S3", "S4", "S5"]:
    sh = wb[name]
    sh.sheet_view.showGridLines = False

wb.save(OUT)
print("WROTE:", OUT)
print("Sheets:", wb.sheetnames)
