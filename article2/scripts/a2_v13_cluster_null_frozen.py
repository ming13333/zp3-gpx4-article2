# -*- coding: utf-8 -*-
"""
A2 v13 冻结口径复算：为正文 P0 修改提供与冻结值同口径的新数字。

关键：冻结管线 transcriptome_audit_a2.py 用 scipy.stats.mannwhitneyu 逐基因检验
（含结值校正 + 连续性校正），得 832 / 4,334 / 843。
本脚本：
  【1】用 scipy 逐基因复算 A / On-only / B / C 四框架的观测显著数（冻结口径）
  【2】实现向量化的「结值+连续性」校正正态近似，验证能复现 scipy 结果
  【3】用向量化版做患者聚类置换，建立各框架的有效零分布
输出：results/a2_v13_cluster_null_frozen.json + 控制台表
"""
import csv, json, os, sys, time
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding="utf-8")

# --- portable paths: everything resolves relative to this file's repository ---
_HERE = os.path.dirname(os.path.abspath(__file__))   # <repo>/article2/scripts
A2    = os.path.dirname(_HERE)                       # <repo>/article2
ROOT  = os.path.dirname(A2)                          # <repo>
FPKM  = os.environ.get("A2_GSE91061_FPKM", os.path.join(
        ROOT, "output", "gse91061_validation",
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"))
MASKS = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
OUT   = os.path.join(A2, "results", "a2_v13_cluster_null_frozen.json")
if not os.path.exists(FPKM):
    raise SystemExit(
        "Input FPKM matrix not found: %s\n"
        "Obtain GSE91061 from GEO and place "
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv at that path, or set A2_GSE91061_FPKM." % FPKM)

# ---------- load ----------
hdr=None; rows=[]; ids=[]
with open(FPKM, newline="", encoding="utf-8") as f:
    r=csv.reader(f); hdr=next(r); samples_all=[h.strip() for h in hdr[1:]]
    for row in r:
        if not row: continue
        ids.append(row[0].strip().strip('"')); rows.append([float(x) for x in row[1:]])
X=np.array(rows,float)
sid={s:i for i,s in enumerate(samples_all)}
masks=json.load(open(MASKS,encoding="utf-8"))
lab=masks["labels"]; A=masks["A_pre"]; B=masks["B_naive_pooled"]; C=masks["C_patient_pre_else_on"]
ON=[s for s in B if "_On_" in s]
print("matrix:",X.shape)

# ---------- [1] scipy 冻结口径观测计数 ----------
def fw_test_scipy(cols):
    ia=[sid[s] for s in cols if lab[s]=="R"]; ib=[sid[s] for s in cols if lab[s]=="N"]
    x1=X[:,ia]; x2=X[:,ib]; n1,n2=len(ia),len(ib)
    ps=np.full(X.shape[0],np.nan); rs=np.full(X.shape[0],np.nan)
    for g in range(X.shape[0]):
        a,b=x1[g],x2[g]
        if np.all(a==a[0]) and np.all(b==b[0]) and a[0]==b[0]:
            ps[g]=1.0; rs[g]=0.0; continue
        try:
            u,p=stats.mannwhitneyu(a,b,alternative="two-sided")
            ps[g]=p; rs[g]=2*u/(n1*n2)-1
        except ValueError: pass
    return ps,rs,n1,n2

FR={"A_pre":A,"On_only":ON,"B_naive":B,"C_patient":C}
obs={}
t0=time.time()
for k,cols in FR.items():
    ps,rs,n1,n2=fw_test_scipy(cols)
    obs[k]=dict(n=len(cols),nR=n1,nN=n2,sig=int(np.nansum(ps<0.05)))
    print("  [scipy 冻结口径] %-10s n=%-4d (R=%d/N=%d)  P<0.05 = %d" % (k,len(cols),n1,n2,obs[k]["sig"]))
print("  scipy 耗时 %.1fs" % (time.time()-t0))

# ---------- [2] 向量化「结值+连续性」校正，验证复现 ----------
def vec_pvals(cols):
    idx=[sid[c] for c in cols]; sub=X[:,idx]; n=sub.shape[1]
    labv=np.array([lab[c] for c in cols])
    Rk=stats.rankdata(sub,axis=1)
    # sum(t^3-t) per gene
    sumt=np.zeros(sub.shape[0])
    S=np.sort(sub,axis=1)
    diff=np.diff(S,axis=1)
    # 向量化结值组：用 unique per gene 太慢 → 用排序后等值段长度
    for g in range(sub.shape[0]):
        _,cnt=np.unique(S[g],return_counts=True); sumt[g]=np.sum(cnt**3-cnt)
    corr=(n+1)-sumt/(n*(n-1)); corr=np.maximum(corr,1e-9)
    def pv(maskR, cont=True):
        nR=int(maskR.sum()); nN=n-nR
        U=Rk[:,maskR].sum(axis=1)-nR*(nR+1)/2
        mu=nR*nN/2; sd=np.sqrt(nR*nN/12*corr)
        d=np.abs(U-mu)-(0.5 if cont else 0.0); d=np.maximum(d,0)
        return 2*(1-stats.norm.cdf(d/sd))
    return pv,labv

print("\n  向量化复现校验（连续性校正 on/off）:")
ver={}
for k,cols in FR.items():
    pv,labv=vec_pvals(cols)
    m=(labv=="R")
    c_on=int(np.sum(pv(m,True)<0.05)); c_off=int(np.sum(pv(m,False)<0.05))
    ver[k]=dict(cont_on=c_on,cont_off=c_off)
    print("     %-10s scipy=%-5d  向量化(连续性校正)=%-5d  向量化(无校正)=%-5d  偏差=%+.1f%%" % (
        k,obs[k]["sig"],c_on,c_off,100*(c_on-obs[k]["sig"])/max(obs[k]["sig"],1)))

# ---------- [3] 患者聚类置换零分布 ----------
NPERM=500
res={"frozen_caliber_observed":obs,"vectorized_check":ver,"perm":{}}
print("\n  患者聚类置换（%d 次/框架）:" % NPERM)
for k,cols in FR.items():
    pv,labv=vec_pvals(cols)
    obs_v=int(np.sum(pv(labv=="R",True)<0.05))
    pt=np.array([c.split("_")[0] for c in cols])
    pids=list(dict.fromkeys(pt.tolist()))
    plab={p:lab[[c for c in cols if c.split("_")[0]==p][0]] for p in pids}
    nRp=sum(1 for p in pids if plab[p]=="R")
    rng=np.random.default_rng(20261014); null=np.empty(NPERM,int)
    for i in range(NPERM):
        perm=rng.permutation(len(pids))
        Rset=set(np.array(pids)[perm][:nRp].tolist())
        null[i]=int(np.sum(pv(np.isin(pt,list(Rset)),True)<0.05))
    pP=float(np.mean(null>=obs_v))
    res["perm"][k]=dict(n=len(cols),n_patients=len(pids),nR_patients=nRp,
        observed_vectorized=obs_v, observed_scipy=obs[k]["sig"],
        null_mean=float(null.mean()), null_p95=float(np.percentile(null,95)),
        null_pct=round(100*float(null.mean())/X.shape[0],2), perm_P=pP,
        excess=obs_v-int(null.mean()))
    print("     %-10s 观测=%-5d(向量化)/%-5d(scipy)  零分布均值=%-6.0f (%.1f%%)  置换P=%.4f  超额=%+d" % (
        k,obs_v,obs[k]["sig"],null.mean(),res["perm"][k]["null_pct"],pP,res["perm"][k]["excess"]))

json.dump(res,open(OUT,"w",encoding="utf-8"),ensure_ascii=False,indent=1)
print("\nsaved:",OUT)
