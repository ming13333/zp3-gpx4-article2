# -*- coding: utf-8 -*-
"""
Update A2 supplementary tables for v14 (TO submission).

 - S9: rebuild with the 10,000-permutation nulls (replaces the 500-permutation version)
 - S10: NEW table - the fixed-42-patient isolation experiment

Starts from the v13 workbook so that S1-S8 are carried over byte-identically (their
underlying analyses are unchanged; re-running the S5/S6 bootstrap machinery would only
re-derive the same values). Reads results/a2_v14_cluster_null_frozen.json, which is the
frozen-schema view of results/a2_v14_paired_bootstrap.json.

Writes results/A2_Supplementary_Tables_v14.xlsx.
"""
import json, os, sys
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
RES = os.path.join(A2, "results")
SRC = os.path.join(RES, "A2_Supplementary_Tables_v13.xlsx")
DST = os.path.join(RES, "A2_Supplementary_Tables_v14.xlsx")
NULLJ = os.path.join(RES, "a2_v14_cluster_null_frozen.json")

if not os.path.exists(SRC):
    raise SystemExit("Base workbook not found: %s" % SRC)

wb = openpyxl.load_workbook(SRC)
BOLD = Font(bold=True)
WRAP = Alignment(wrap_text=True, vertical="top")
HDR = PatternFill("solid", fgColor="DDEBF7")

N = json.load(open(NULLJ, encoding="utf-8"))
PERM = N["perm"]
F42 = N["perm_fixed42"]
SEED = N["_provenance"]["seed"]

LAB = {"A_pre": "A - Pre-treatment only",
       "On_only": "On-treatment only",
       "B_naive": "B - Naive pooled (pre + on)",
       "C_patient": "C - Patient-level aggregated (pre-else-on)"}
LAB42 = {"P_pre42": "P_pre42 - pre-treatment only (fixed 42 patients)",
         "P_on42": "P_on42 - on-treatment only (fixed 42 patients)",
         "P_pool84": "P_pool84 - naive pooled, 84 samples (fixed 42 patients)"}

# ============================================================ S9 (rebuilt at 10,000)
if "S9" in wb.sheetnames:
    del wb["S9"]
ws9 = wb.create_sheet("S9")
ws9.cell(row=1, column=1, value="Table S9. Patient-cluster permutation null distributions for the GSE91061 transcriptome-wide audit (22,187 gene rows of the matrix)").font = Font(bold=True, size=12)
ws9.cell(row=2, column=1, value=(
    "Null distributions were built from 10,000 permutations of patient labels (fixed seed %d); the sample-level Mann-Whitney U statistic was "
    "recomputed on each permutation, preserving the within-patient repeated-measure structure and the number of responder patients. Per-gene "
    "P values use the tie-corrected, continuity-corrected normal approximation, which reproduces the frozen per-gene scipy.stats.mannwhitneyu "
    "counts exactly (832 / 4,015 / 4,334 / 843 at P<0.05; 0.0%% deviation). The permutation P is (b+1)/(B+1) with b the number of permutations "
    "reaching the observed statistic. The 'null share' column is the null distribution mean expressed as a percentage of the nominal 5%% "
    "threshold: a framework above 5%% is anti-conservative." % SEED))
ws9.cell(row=2, column=1).alignment = WRAP
ws9.row_dimensions[2].height = 62

cols9 = ["Framework", "n samples", "n patients", "n responder patients",
         "Observed genes P<0.05 (scipy)", "Observed genes P<0.05 (vectorised)",
         "Null mean (genes)", "Null SD (genes)", "Null share of 5% threshold (%)",
         "Permutation P", "Tail hits (of 10,000)", "Excess over null"]
for j, c in enumerate(cols9, start=1):
    cell = ws9.cell(row=4, column=j, value=c)
    cell.font = BOLD; cell.fill = HDR; cell.alignment = WRAP

for i, k in enumerate(["A_pre", "On_only", "B_naive", "C_patient"]):
    p = PERM[k]
    vals = [LAB[k], p["n"], p["n_patients"], p["nR_patients"],
            p["observed_scipy"], p["observed_vectorized"],
            round(p["null_mean"], 1), round(p["null_sd"], 1), p["null_pct"],
            round(p["perm_P"], 4), p["tail_hits"], p["excess"]]
    for j, v in enumerate(vals, start=1):
        ws9.cell(row=5 + i, column=j, value=v)

ws9.cell(row=10, column=1, value=(
    "Reading: under a true null the same sample-level test occupies 4.5-4.9% of the nominal 5% threshold in Frameworks A, On_only and C, but 10.1% under naive "
    "pooling (B) - i.e. pooling approximately doubles the false-positive floor. The rise in the observed count from 832 to 4,334 is nevertheless dominated by genuine "
    "on-treatment biology: the On_only framework alone yields 4,015 genes on a valid 4.9% null floor (permutation P=0.0016; excess +2,929). Pooling raises the "
    "false-positive floor and, in the same step, reduces the excess over the null to +2,083. Neither A (excess -159; P=0.552) nor C (excess -190; P=0.572) shows any "
    "excess over its null. The reduction in excess is directionally consistent but its 95% confidence interval includes zero (see Table S10 and the Results), so it is "
    "reported as exploratory rather than as an established second effect. Resolution note: with 10,000 permutations the smallest attainable permutation P is 1e-4, and "
    "the reported values still carry Monte Carlo error."))
ws9.cell(row=10, column=1).alignment = WRAP
ws9.row_dimensions[10].height = 74
for col, w in zip("ABCDEFGHIJKL", [46, 10, 11, 12, 16, 16, 12, 12, 14, 12, 14, 12]):
    ws9.column_dimensions[col].width = w

# ============================================================ S10 (new)
if "S10" in wb.sheetnames:
    del wb["S10"]
idx = wb.sheetnames.index("S9") + 1 if "S9" in wb.sheetnames else len(wb.sheetnames)
ws10 = wb.create_sheet("S10", idx)
ws10.cell(row=1, column=1, value="Table S10. Isolation analysis: the effect of the analysis unit, with patient composition held fixed (GSE91061)").font = Font(bold=True, size=12)
ws10.cell(row=2, column=1, value=(
    "Design. Forty-two of the 62 response-evaluable GSE91061 patients contributed exactly one pre-treatment and one on-treatment sample (9 responders / 33 non-responders; "
    "the pre/on response labels are concordant within every patient). The three frameworks below use only those 42 patients, so patient composition, response labels and "
    "the 22,187-gene row set are identical across all three; the only thing that changes is the analysis unit. P_pre42 and P_on42 each analyse 42 samples (one per patient, "
    "no within-patient clustering); P_pool84 analyses the same 84 samples pooled as if independent, so every patient contributes two correlated observations. Null "
    "distributions are from 10,000 patient-cluster permutations with the same tie- and continuity-corrected implementation used throughout (seed %d)." % SEED))
ws10.cell(row=2, column=1).alignment = WRAP
ws10.row_dimensions[2].height = 76

ws10.cell(row=4, column=1, value=(
    "Interpretation. With the patient set and the gene set both fixed, moving the analysis unit from on-treatment-only to naive pooling raises the null share of the "
    "nominal 5% threshold from 4.52% to 11.32% (a 2.50-fold increase), which is even higher than the 10.1% seen in the full pooled framework B. The inflation of the "
    "false-positive floor is therefore a property of pooling repeated measurements, not of any difference in patient or time-point composition between frameworks."))
ws10.cell(row=4, column=1).alignment = WRAP
ws10.row_dimensions[4].height = 60

cols10 = ["Framework", "n samples", "n patients", "n responder patients",
          "Observed genes P<0.05", "Null mean (genes)", "Null SD (genes)",
          "Null share of 5% threshold (%)", "Permutation P", "Tail hits (of 10,000)",
          "Excess over null"]
for j, c in enumerate(cols10, start=1):
    cell = ws10.cell(row=6, column=j, value=c)
    cell.font = BOLD; cell.fill = HDR; cell.alignment = WRAP

for i, k in enumerate(["P_pre42", "P_on42", "P_pool84"]):
    p = F42[k]
    vals = [LAB42[k], p["n"], p["n_patients"], p["nR_patients"],
            p["observed_vectorized"], round(p["null_mean"], 1), round(p["null_sd"], 1),
            p["null_pct"], round(p["perm_P"], 4), p["tail_hits"], p["excess"]]
    for j, v in enumerate(vals, start=1):
        ws10.cell(row=7 + i, column=j, value=v)

ws10.cell(row=11, column=1, value=(
    "Note on the excess-over-null difference. The difference in excess between P_pool84 and P_on42 is -1,085 (-28.4%) in this fixed-42 subset. Quantified across the full "
    "pooled analysis (Part 2/Part 4 of the analysis script): stratified patient bootstrap (B=600) gives -1,517 (95% CI -3,423 to +149) and leave-one-patient-out jackknife "
    "gives -1,053 (SE 1,390; 95% CI -3,776 to +1,671). Both intervals include zero, and 96.0% of bootstrap replicates are negative, so the reduction in recovered signal is "
    "directionally consistent but not statistically resolved at this sample size."))
ws10.cell(row=11, column=1).alignment = WRAP
ws10.row_dimensions[11].height = 62
for col, w in zip("ABCDEFGHIJK", [46, 10, 11, 12, 16, 12, 12, 14, 12, 14, 12]):
    ws10.column_dimensions[col].width = w

wb.save(DST)
print("SAVED", DST)
print("sheets:", wb.sheetnames)
print()
for k in ["A_pre", "On_only", "B_naive", "C_patient"]:
    p = PERM[k]
    print("S9  %-10s n=%-4d obs=%-5d null=%-7.1f sd=%-7.1f share=%.2f%% P=%.4f tail=%s exc=%+d"
          % (k, p["n"], p["observed_vectorized"], p["null_mean"], p["null_sd"],
             p["null_pct"], p["perm_P"], p["tail_hits"], p["excess"]))
for k in ["P_pre42", "P_on42", "P_pool84"]:
    p = F42[k]
    print("S10 %-10s n=%-4d obs=%-5d null=%-7.1f sd=%-7.1f share=%.2f%% P=%.4f tail=%s exc=%+d"
          % (k, p["n"], p["observed_vectorized"], p["null_mean"], p["null_sd"],
             p["null_pct"], p["perm_P"], p["tail_hits"], p["excess"]))
