"""Download GSE163973 supplementary files from NCBI GEO (anonymous FTP) into data/raw/GSE163973 (~3 GB).
Usage: python scripts/fetch_raw_data.py --dry-run   (list files and sizes only)
       python scripts/fetch_raw_data.py             (download; skips files already present)"""
import argparse
from ftplib import FTP
from pathlib import Path

HOST, PATH = "ftp.ncbi.nlm.nih.gov", "/geo/series/GSE163nnn/GSE163973/suppl"
DEST = Path(__file__).resolve().parent.parent / "data" / "raw" / "GSE163973"

ap = argparse.ArgumentParser()
ap.add_argument("--dry-run", action="store_true", help="list remote files and sizes, download nothing")
args = ap.parse_args()

ftp = FTP(HOST, timeout=60); ftp.login()
ftp.cwd(PATH)
names = [n for n in ftp.nlst() if "." in n]  # skip sub-directories (e.g. raw_tar)
assert names, "no files listed; GEO layout may have changed"
ftp.voidcmd("TYPE I")  # nlst() switches to ASCII; SIZE needs binary
sizes = {n: ftp.size(n) for n in names}
for n in names:
    print(f"{n}\t{sizes[n] / 1e6:.1f} MB")
print(f"total {sum(sizes.values()) / 1e9:.2f} GB")
if args.dry_run:
    ftp.quit(); raise SystemExit(0)

DEST.mkdir(parents=True, exist_ok=True)
for n in names:
    out = DEST / n
    if out.exists() and out.stat().st_size == sizes[n]:
        print("skip", n); continue
    tmp = out.with_suffix(out.suffix + ".part")
    with open(tmp, "wb") as f:
        ftp.retrbinary(f"RETR {n}", f.write, blocksize=1 << 20)
    tmp.replace(out)
    print("downloaded", n)
ftp.quit()
