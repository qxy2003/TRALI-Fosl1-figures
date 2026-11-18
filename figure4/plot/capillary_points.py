"""Draw paired mouse-level gCap and aCap readouts as aligned facets."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]
out=root/"figure"
data=root/"data"
ink="#202E34";muted="#53646B";rule="#D6DEE0"
red="#DE2D26";blue="#4285BC"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"text.color":ink,
                     "axes.labelcolor":ink,"axes.edgecolor":"#8A959A","axes.linewidth":.65,
                     "xtick.color":muted,"ytick.color":ink,"pdf.fonttype":42})
def tidy(ax):
    ax.spines[["top","left","right"]].set_visible(False)
    ax.tick_params(axis="y",length=0,pad=4)
    ax.tick_params(axis="x",length=2.8,pad=2,labelsize=8)
def marker(color,shape,label):
    return Line2D([],[],ls="",marker=shape,ms=5.5,mfc=color,mec="white",mew=.4,label=label)
mice=["T1","T2","T3","T4"]
marks=["o","s","^","D"]
mousecols=["#2166AC","#4285BC","#DE2D26","#FC9272"]
types=["gCap","aCap"]
def subtitle(ax,tag):
    ax.set_title(tag,fontsize=9.3,weight="semibold",pad=5)

b=pd.read_csv(data/"4B_capillary_by_mouse.csv")
genes=b.gene.drop_duplicates().tolist()
fig=plt.figure(figsize=(401.1984/72,188.2908/72))
fig.text(.022,.875,"B",fontsize=17.5,weight="bold")
for idx,tag in enumerate(types):
    ax=fig.add_axes([.202+idx*.408,.27,.338,.54]);tidy(ax);subtitle(ax,tag)
    ax.set(xlim=(-.18,.011),xticks=[-.16,-.08,0],ylim=(-.55,5.55),yticks=range(6),
           yticklabels=genes[::-1] if idx==0 else [])
    for lab in ax.get_yticklabels():lab.set_fontstyle("italic")
    ax.axvline(0,lw=.7,color="#A2B0B4",ls=(0,(2,2)))
    for y,gene in zip(range(5,-1,-1),genes):
        if y>0:ax.axhline(y-.5,color=rule,lw=.45,zorder=0)
        v=b[b.gene.eq(gene)&b.subtype.eq(tag)].set_index("mouse").loc[mice]
        if not v.modeled.any():
            ax.text(-.083,y,"Outside model",ha="center",va="center",color=muted,fontsize=7.6)
            continue
        assert v.delta.between(-.18,.011).all()
        for x,offset,s,c in zip(v.delta,[.26,.09,-.09,-.26],marks,mousecols):
            ax.plot(x,y+offset,marker=s,mfc=c,mec=ink,mew=.28,ms=4.1,ls="")
fig.text(.575,.123,"KO − baseline (model log expression)",ha="center",fontsize=8.6)
fig.legend(handles=[marker(c,s,m) for c,s,m in zip(mousecols,marks,mice)],ncol=4,frameon=False,
           loc="center",bbox_to_anchor=(.575,.05),handlelength=.8,handletextpad=.32,columnspacing=1.1,fontsize=8)
fig.savefig(out/"Figure4B.pdf");plt.close(fig)

d=pd.read_csv(data/"4D_capillary_by_mouse.csv")
fig=plt.figure(figsize=(334.0656/72,166.752/72))
fig.text(.037,.872,"D",fontsize=17.5,weight="bold")
for idx,tag in enumerate(types):
    ax=fig.add_axes([.143+idx*.447,.30,.37,.50]);tidy(ax);subtitle(ax,tag)
    v=d[d.subtype.eq(tag)].set_index("mouse").drop(columns="subtype").loc[mice].copy()
    v.loc["Mean"]=v.mean()
    ax.set(xlim=(0,4.65),xticks=[0,2,4],ylim=(-.55,4.55),yticks=[4,3,2,1,0],
           yticklabels=v.index if idx==0 else [])
    if idx==0:ax.get_yticklabels()[-1].set_weight("bold")
    ax.axhline(.55,color=rule,lw=.7)
    for y,(mouse,row) in zip(range(4,-1,-1),v.iterrows()):
        x1,x2=row.Full_panel,row.Fig3_shared_up
        ax.plot([x1,x2],[y,y],color="#B9C4C7",lw=1.1)
        ax.scatter(x1,y,marker="s",s=22,c=blue,edgecolors="white",linewidths=.4,zorder=3)
        ax.scatter(x2,y,marker="o",s=22,c=red,edgecolors="white",linewidths=.4,zorder=3)
fig.text(.554,.146,"Distance reduction (%)",ha="center",fontsize=8.6)
fig.legend(handles=[marker(red,"o","Shared candidates (81)"),marker(blue,"s","All non-Fosl1 (2,999)")],
           loc="center",bbox_to_anchor=(.52,.056),frameon=False,ncol=2,handlelength=.8,
           handletextpad=.3,columnspacing=.9,fontsize=7.9)
fig.savefig(out/"Figure4D.pdf");plt.close(fig)
