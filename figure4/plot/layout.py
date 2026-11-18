"""Build the vector layout. Native KO/control drawing is embedded by assemble.py."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]
data=root/"data"
dest=root/"plot/intermediate"
dest.mkdir(parents=True,exist_ok=True)
(root/"figure").mkdir(exist_ok=True)
W,H=11.1,9.65
ink="#202E34"; mute="#53646B"; grey="#B9C4C7"; line="#D6DEE0"; blue="#53889C"; red="#A74749"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"text.color":ink,"axes.labelcolor":ink,
    "axes.edgecolor":"#8A959A","axes.linewidth":.65,"xtick.color":mute,"ytick.color":ink,
    "pdf.fonttype":42,"xtick.labelsize":8,"ytick.labelsize":9})
fig=plt.figure(figsize=(W,H),facecolor="white")
def txt(x,y,text,**kw):
    return fig.text(x,y,text,**kw)
def panel_letter(x,y,letter):
    txt(x,y,letter,weight="bold",fontsize=17.5,va="baseline")
def clean(ax):
    ax.spines[["left","right","top"]].set_visible(False)
    ax.tick_params(axis="y",length=0,pad=5)
    ax.tick_params(axis="x",length=3,width=.6)
def key(x,y,markers,labels,cols,**kwargs):
    hs=[Line2D([],[],marker=m,mfc=c,mec="white",mew=.45,ms=6,linestyle="") for m,c in zip(markers,cols)]
    fig.legend(hs,labels,bbox_to_anchor=(x,y),loc="center",frameon=False,fontsize=8,
               handlelength=.8,handletextpad=.37,columnspacing=1,**kwargs)
mice=["T1","T2","T3","T4"]; marks=["o","s","^","D"]
mousecols=["#3A6571","#5F8594","#A35B44","#BB8872"]

panel_letter(.026,.903,"A")
table=pd.read_csv(data/"4A_endothelial_rank.csv").query("significant").sort_values("rank")
ax=fig.add_axes([.105,.745,.259,.154]); clean(ax)
ys=np.arange(len(table))[::-1]
ax.barh(ys,table.distance*1e6,height=.56,color=[red if f else "#72929C" for f in table.Fig3B])
ax.set(yticks=ys,yticklabels=table.gene.tolist(),xlim=(0,3.32),ylim=(-.53,4.56),xticks=[0,1,2,3],
       xlabel=r"Network perturbation distance ($\times 10^{-6}$)")
for lab,hit in zip(ax.get_yticklabels(),table.Fig3B):
    lab.set_fontstyle("italic")
    if hit: lab.set_color(red)
def format_q(value):
    if value>=.001: return f"{value:.4f}"
    base,exponent=f"{value:.2e}".split("e")
    return rf"${base}\times10^{{{int(exponent)}}}$"
qfmt=table["p.adj"].map(format_q)
for y,q in zip(ys,qfmt):
    ax.text(1.25,y,q,transform=ax.get_yaxis_transform(),ha="center",va="center",fontsize=8.2)
ax.text(1.25,4.87,"BH q",transform=ax.get_yaxis_transform(),ha="center",color=mute,fontsize=8)
txt(.063,.681,"5 of 1,499 downstream genes at BH q < 0.05",fontsize=8.4,color=mute)

# Accent colors used in the Fig. 1–3 plotting scripts.
accent_red = "#DE2D26"
accent_blue = "#4285BC"
panel_letter(.442,.214,"E")
scores=pd.read_csv(data/"4E_pathways_by_mouse.csv")
tests=pd.read_csv(data/"4E_paired_tests.csv",keep_default_na=False).set_index("module")
labels={"TightJunction":"Tight junction","LeukocyteAdhesion":"Leukocyte adhesion",
        "ActinMyosinSliding":"Actin-myosin\nfilament sliding","Complement":"Complement",
        "Coagulation":"Coagulation","Fibrinolysis":"Fibrinolysis"}
matrix=scores.pivot(index="module",columns="mouse",values="delta").loc[list(labels),mice]*1000
assert np.isfinite(matrix.to_numpy()).all() and set(tests.index)==set(labels)
ys=np.arange(len(matrix),0,-1)
ax=fig.add_axes([.587,.067,.383,.150])
clean(ax)
ax.axvline(0,color="#86979B",lw=.8,ls=(0,(2.6,2.3)),zorder=0)
boxplot=ax.boxplot(matrix.to_numpy().T,positions=ys,widths=.38,orientation="horizontal",whis=(0,100),
           manage_ticks=False,showfliers=False,patch_artist=True,
           boxprops=dict(linewidth=.95),medianprops=dict(linewidth=1.3),
           whiskerprops=dict(linewidth=.9),capprops=dict(linewidth=.9))
for i,(module,vals) in enumerate(matrix.iterrows()):
    pos=ys[i]
    median=vals.median()
    if median > 0:
        accent,fill=accent_red,"#FEE5D9"
    elif median < 0:
        accent,fill=accent_blue,"#E5F0F9"
    else:
        accent,fill=mute,"#F7F7F7"
    boxplot["boxes"][i].set(facecolor=fill,edgecolor=accent)
    boxplot["medians"][i].set_color(accent)
    for component in ("caps","whiskers"):
        for artist in boxplot[component][2*i:2*i+2]:
            artist.set_color(accent)
    for x,offset,mark in zip(vals,[-.28,-.14,.14,.28],marks):
        ax.plot(x,pos+offset,marker=mark,color=accent,ms=3.8,ls="",mec="white",mew=.35,zorder=4)
    if tests.loc[module,"figure_symbol"]:
        ax.annotate(tests.loc[module,"figure_symbol"],(vals.max(),pos),xytext=(5,0),
                    textcoords="offset points",va="center",fontsize=10)
ax.set(yticks=ys,yticklabels=list(labels.values()),ylim=(.49,6.51),xlim=(-47,9),
       xticks=[-40,-30,-20,-10,0],xlabel="1,000 × (KO − baseline score)")
ax.tick_params(axis="y",labelsize=8.4,pad=5)
key(.785,.010,marks,mice,[mute]*4,ncol=4)
fig.savefig(dest/"base.pdf")
plt.close(fig)
(dest/"geometry.json").write_text(json.dumps({"panels":{
    "A":[.011,.665,.484,.271], "E":[.433,.004,.560,.240]}}))
