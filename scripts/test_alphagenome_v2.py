"""One-variant smoke test of the AlphaGenome 2 API: rs873549 (GRCh38 chr1:222098425),
scored in dermal/skin fibroblast and keratinocyte conditions."""
import os

import requests
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c, variant_scorers as vs

load_dotenv(os.path.join(os.path.dirname(__file__), "..", ".env"))
cl = c.create(os.environ["ALPHAGENOME_API_KEY"])

conds = cl.supported_conditions()
kera = conds[conds.biosample_name.str.contains("keratinocyte", case=False)]
print("keratinocyte conditions:", kera[["condition", "biosample_name"]].head(5).to_string(index=False))
FIBRO = "UBERON:0002097_CL:0002620"  # skin of body skin fibroblast
KERA = "UBERON:0002097_CL:0000312"  # skin of body keratinocyte (not the first match, which was stomach)
assert KERA in set(conds.condition)

CHROM, POS = "chr1", 222098425
ref = requests.get(f"https://rest.ensembl.org/sequence/region/human/1:{POS}..{POS}:1",
                   headers={"Content-Type": "text/plain"}, timeout=30).text.strip()
alt = "T" if ref == "C" else "C"  # dbSNP alleles are C/T
print("ref base (GRCh38):", ref, "-> alt", alt)
assert ref in "CT", "unexpected reference base; check coordinates"

interval = genome.Interval(CHROM, POS - 262144, POS + 262144)  # 524288 = SEQUENCE_LENGTH_500KB
variant = genome.Variant(CHROM, POS, ref, alt)
scorers = vs.get_recommended_scorers(c.Organism.HOMO_SAPIENS) if hasattr(vs, "get_recommended_scorers") else ()
for name, cond in [("fibroblast", FIBRO), ("keratinocyte", KERA)]:
    if cond is None:
        print(name, ": no condition found"); continue
    res = cl.score_variant(interval, variant, condition=cond, variant_scorers=scorers)
    print(f"{name} ({cond}): {len(res)} scorer result tables; first shape:",
          res[0].shape if hasattr(res[0], "shape") else type(res[0]))
print("OK")
