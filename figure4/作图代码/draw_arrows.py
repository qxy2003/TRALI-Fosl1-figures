"""Draw the two native CellOracle arrow layers from the released numeric flow grid."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

root = Path(__file__).resolve().parents[1]
dest = root / "作图代码" / "intermediate" / "4C"
dest.mkdir(parents=True, exist_ok=True)
for subtype, folder, xlim, ylim in [
    ("gCap","4C_拟时序",(-6.4,5.6),(-7.5,3.5)),
    ("aCap","4C_aCap",(-5.9,6.1),(-2.1,8.9))]:
    grid = pd.read_csv(root / "分析结果" / folder / "grid_reference_and_fields.csv")
    flow = grid.loc[grid["visible"]]
    assert len(flow) > 0 and flow["ref_supported"].all()
    with PdfPages(dest / f"{subtype}_arrows.pdf") as pdf:
        for dx, dy in [("dx", "dy"), ("random_dx", "random_dy")]:
            fig = plt.figure(figsize=(3.60, 3.30))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.quiver(flow.x, flow.y, flow[dx], flow[dy], angles="xy",
                      scale_units="xy", scale=1.30, width=.0047, headwidth=3.35,
                      headlength=4.5, headaxislength=3.9, minlength=0,
                      color="#253E46", pivot="tail")
            ax.set(xlim=xlim, ylim=ylim, aspect="equal")
            ax.set_axis_off()
            pdf.savefig(fig, transparent=True)
            plt.close(fig)
