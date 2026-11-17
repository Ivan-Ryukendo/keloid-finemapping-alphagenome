"""Expand keloid lead SNPs (p<5e-8 + PHLDA3 pair) into LD proxies (r2>=0.8, EUR/EAS, GRCh38)
via Ensembl REST; annotate with batched VEP. Raw JSON cached in cache/finemap/. Re-runnable."""
import json, time, pathlib, requests
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
CACHE = ROOT / "cache" / "finemap"; CACHE.mkdir(parents=True, exist_ok=True)
S = "https://rest.ensembl.org"
POPS = {"EUR": "1000GENOMES:phase_3:EUR", "EAS": "1000GENOMES:phase_3:EAS"}
failures = []

def call(method, path, cache_name, **kw):
    f = CACHE / cache_name
    if f.exists():
        return json.loads(f.read_text())
    for i in range(6):
        r = requests.request(method, S + path, headers={"Content-Type": "application/json", "Accept": "application/json"}, timeout=120, **kw)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(float(r.headers.get("Retry-After", 2 ** i))); continue
        if r.status_code == 400:  # e.g. rsID has no LD data / not found
            failures.append((cache_name, r.status_code, r.text[:100])); out = []
        else:
            r.raise_for_status(); out = r.json()
        f.write_text(json.dumps(out)); time.sleep(0.1)
        return out
    failures.append((cache_name, "retries exhausted", "")); return None

d = pd.read_csv(ROOT / "data/raw/keloid_gwas_variants.csv")
ok = d[d.rs_id.str.startswith("rs") & d.in_gwas_catalog_keloid & (d.p_value < 5e-8)]
leads = sorted(set(ok.rs_id) | {"rs35383942", "rs192314256"})

ld = {}  # (lead, proxy) -> {pop: r2}
for lead in leads:
    for pop, full in POPS.items():
        res = call("GET", f"/ld/human/{lead}/{full}?r2=0.8", f"ld_{lead}_{pop}.json")
        for h in res or []:
            ld.setdefault((lead, h["variation2"]), {})[pop] = float(h["r2"])
        ld.setdefault((lead, lead), {})[pop] = 1.0  # lead itself

ids = sorted({p for _, p in ld})
vep = {}
for i in range(0, len(ids), 200):
    res = call("POST", "/vep/human/id", f"vep_{i}.json", data=json.dumps({"ids": ids[i:i + 200]}))
    for v in res or []:
        vep[v["input"]] = v

rows, excl = [], []
for (lead, p), r2 in sorted(ld.items()):
    v = vep.get(p)
    if v is None:
        excl.append(dict(lead_rsid=lead, proxy_rsid=p, reason="no VEP/Ensembl record")); continue
    al = v["allele_string"].split("/")
    row = dict(lead_rsid=lead, proxy_rsid=p, chrom=v["seq_region_name"], pos_grch38=v["start"], ref=al[0],
               alt=",".join(al[1:]), r2_EUR=r2.get("EUR"), r2_EAS=r2.get("EAS"), consequence=v["most_severe_consequence"])
    if any(len(a) != 1 or a == "-" for a in al) or v["start"] != v["end"]:
        excl.append({**row, "reason": "indel/non-SNV"})
    elif len(al) != 2:
        excl.append({**row, "reason": "multiallelic"})
    else:
        rows.append(row)

out = pd.DataFrame(rows); ex = pd.DataFrame(excl)
out.to_csv(ROOT / "data/processed/finemap_candidates.csv", index=False)
ex.to_csv(ROOT / "data/processed/finemap_excluded.csv", index=False)

assert out.ref.str.len().eq(1).all() and out.alt.str.len().eq(1).all(), "non-SNV leaked"
assert out.r2_EUR.fillna(1).ge(0.8).all() and out.r2_EAS.fillna(1).ge(0.8).all(), "r2 below threshold"
print("leads", len(leads), "kept", len(out), "excluded", len(ex), "failures", failures)
print(out.groupby("lead_rsid").size().to_string())
print("leads with zero kept:", sorted(set(leads) - set(out.lead_rsid)))
