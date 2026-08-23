#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A2 参考文献 CrossRef 程序化核验：逐条查询 DOI，比对标题/年份，报告错配。"""
import re, json, sys, urllib.request, urllib.error, time

SRC = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\manuscripts\Article2_A2_Cancers格式_v0.7.md"
OUT = r"C:\D\workbuddy\科研\细胞外GPX4免疫抑制\article2\results\ref_crossref_verify_a2.csv"

def norm(s):
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return s.strip()

def extract_refs(text):
    refs = []
    # 匹配 "NN. ... https://doi.org/xxx"
    pat = re.compile(r"^(\d{2})\.\s+(.*?)\s+(https?://(?:dx\.)?doi\.org/[^\s]+)\s*$", re.M)
    for m in pat.finditer(text):
        num = int(m.group(1))
        body = m.group(2)
        doi = m.group(3)
        # 标题：第一个句点后到 *Journal* 前
        title_m = re.search(r"^\s*[^.]+\.\s+(.*?)\.\s+\*", body)
        title = title_m.group(1).strip() if title_m else ""
        year_m = re.search(r"\*\*(\d{4})\*\*", body)
        year = year_m.group(1) if year_m else ""
        refs.append({"num": num, "title": title, "year": year, "doi": doi})
    return refs

def query_crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": "A2-ref-verify/1.0 (mailto:lqin5467@gmail.com)"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            data = json.load(r)
        msg = data.get("message", {})
        x_title = " ".join(msg.get("title", [])) if msg.get("title") else ""
        x_year = ""
        if msg.get("issued") and msg["issued"].get("date-parts"):
            parts = msg["issued"]["date-parts"][0]
            if parts:
                x_year = str(parts[0])
        return x_title, x_year, "OK"
    except urllib.error.HTTPError as e:
        return "", "", f"HTTP {e.code}"
    except Exception as e:
        return "", "", f"ERR {type(e).__name__}: {e}"

def main():
    with open(SRC, encoding="utf-8") as f:
        text = f.read()
    refs = extract_refs(text)
    print(f"解析到 {len(refs)} 条参考文献", file=sys.stderr)
    rows = []
    for rf in refs:
        xt, xy, status = query_crossref(rf["doi"])
        tsim = norm(rf["title"]) == norm(xt) if xt else False
        # 标题词重叠
        a = set(norm(rf["title"]).split()); b = set(norm(xt).split())
        overlap = len(a & b) / max(1, len(a | b)) if (a or b) else 0
        ysim = (rf["year"] == xy) if xy else None
        flag = "OK"
        if status != "OK":
            flag = "DOI_FAIL"
        elif not tsim and overlap < 0.5:
            flag = "TITLE_MISMATCH"
        elif ysim is False:
            flag = "YEAR_MISMATCH"
        rows.append({
            "num": rf["num"], "doi": rf["doi"], "stated_year": rf["year"],
            "crossref_year": xy, "title_overlap": round(overlap, 3),
            "status": status, "flag": flag,
            "stated_title": rf["title"][:90], "crossref_title": xt[:90]
        })
        print(f"[{rf['num']:02d}] {flag} yr={rf['year']}/{xy} ovl={overlap:.2f} {rf['doi']}", file=sys.stderr)
        time.sleep(0.3)
    import csv
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["num","doi","stated_year","crossref_year","title_overlap","status","flag","stated_title","crossref_title"])
        w.writeheader(); w.writerows(rows)
    bad = [r for r in rows if r["flag"] != "OK"]
    print(f"\n核验完成：{len(rows)} 条，异常 {len(bad)} 条", file=sys.stderr)
    for r in bad:
        print(f"  -> [{r['num']:02d}] {r['flag']} | stated={r['stated_title']} | crossref={r['crossref_title']}", file=sys.stderr)

if __name__ == "__main__":
    import urllib.parse
    main()
