"""Build AlphaFold Server job files (wild-type vs missense variant) for keloid coding variants.
Protein changes come from Ensembl VEP (canonical transcript); sequences from UniProt (reviewed).
Outputs: data/processed/alphafold_server_jobs.json and per-protein FASTA files.
The residue at each variant position must equal the reference amino acid or the script stops."""
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)
MAX_TOKENS = 5000  # AlphaFold Server per-job limit (official FAQ; re-check on the site)

# gene -> (UniProt accession, [(rsID, ref aa, position, alt aa)]) from Ensembl VEP canonical transcripts
TARGETS = {
    "PHLDA3": ("Q9Y5J5", [("rs35383942", "R", 28, "Q"), ("rs192314256", "E", 62, "G")]),
    "MYO7A": ("Q13402", [("rs35641839", "V", 679, "I")]),
}


def uniprot_seq(acc):
    t = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60).text
    return "".join(t.splitlines()[1:])


def job(name, seq):
    assert len(seq) <= MAX_TOKENS, f"{name}: {len(seq)} residues exceeds token limit"
    return {"name": name, "modelSeeds": [], "sequences": [{"proteinChain": {"sequence": seq, "count": 1}}],
            "dialect": "alphafoldserver", "version": 1}


jobs, fasta = [], []
for gene, (acc, muts) in TARGETS.items():
    wt = uniprot_seq(acc)
    jobs.append(job(f"{gene}_WT_{acc}", wt))
    fasta.append(f">{gene}_WT_{acc}\n{wt}")
    for rs, ref, pos, alt in muts:
        assert wt[pos - 1] == ref, f"{gene} p.{ref}{pos}{alt}: UniProt has {wt[pos - 1]} at {pos}"
        mut = wt[:pos - 1] + alt + wt[pos:]
        n = f"{gene}_{ref}{pos}{alt}_{rs}"
        jobs.append(job(n, mut))
        fasta.append(f">{n}\n{mut}")
        print(f"{n}: residue check OK, length {len(mut)}")

(OUT / "alphafold_server_jobs.json").write_text(json.dumps(jobs, indent=1))
(OUT / "protein_sequences.fasta").write_text("\n".join(fasta) + "\n")
print(f"{len(jobs)} jobs written to {OUT / 'alphafold_server_jobs.json'}")
