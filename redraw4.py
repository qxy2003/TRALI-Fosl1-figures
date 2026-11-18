"""Use this Python for Matplotlib; optionally use another for PyMuPDF/Pillow."""
import argparse
from pathlib import Path
import subprocess
import sys

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--assembly-python",default=sys.executable)
args=parser.parse_args()
root=Path(__file__).resolve().parent
for script in ("layout","capillary_points","plot_4C","draw_arrows","assemble"):
    executable=args.assembly_python if script=="assemble" else sys.executable
    subprocess.run([executable, str(root/"figure4"/"plot"/(script+".py"))],cwd=root,check=True)
