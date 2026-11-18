"""Check the distribution; add --figures after redrawing the four combined plates."""
import argparse, csv, hashlib
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--figures",action="store_true")
args=p.parse_args()
root=Path(__file__).resolve().parent
for line in (root/"SHA256SUMS").read_text(encoding="utf-8").splitlines():
    expected, name=line.split("  ",1)
    if hashlib.sha256((root/name).read_bytes()).hexdigest()!=expected:
        raise SystemExit(f"Input hash mismatch: {name}. Use the original file bytes (including line endings).")
print("Release files match SHA256SUMS.")
if args.figures:
    import pymupdf
    for r in csv.DictReader((root/"environment"/"reference_pixels.csv").open(encoding="utf-8")):
        i=r["figure"]
        with pymupdf.open(root/f"figure{i}/figure/Figure{i}.pdf") as d:
            if len(d)!=1: raise SystemExit(f"Unexpected page count: Figure {i}")
            pix=d[0].get_pixmap(dpi=150,colorspace=pymupdf.csRGB,alpha=False)
            if (pix.width,pix.height,hashlib.sha256(pix.samples).hexdigest())!=(
                    int(r["width"]),int(r["height"]),r["sha256_rgb"]):
                raise SystemExit(f"Figure {i}: rendered pixels differ from the reference; check fonts and software versions.")
            print(f"Figure {i}: all {pix.width} x {pix.height} rendered pixels match.")
