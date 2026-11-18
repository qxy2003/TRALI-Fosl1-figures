"""Draw two independently fitted state axes and align their native KO and null fields."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]
dest=root/"作图代码/intermediate/4C"
dest.mkdir(parents=True,exist_ok=True)
ink="#24323A";soft="#57666D"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"text.color":ink,
                    "axes.labelcolor":ink,"axes.edgecolor":"#94A1A4","axes.linewidth":.6,
                    "xtick.color":soft,"ytick.color":ink,"pdf.fonttype":42})
specs=[("gCap","4C_拟时序",(-6.4,5.6),(-7.5,3.5)),("aCap","4C_aCap",(-5.9,6.1),(-2.1,8.9))]
seq=colors.LinearSegmentedColormap.from_list("pseudotime",["#DCEBF2","#97C7E0","#6398C6","#7073AF","#613F88","#3F205E"])
div=colors.LinearSegmentedColormap.from_list("signed",["#A85D94","#DBB5D3","#F1F1EB","#ACD4BF","#4D9474"])
W,H=11.1,5.65
fig=plt.figure(figsize=(W,H))
cols=[.045,.195,.345,.507,.657]
width=.142
height=width*W/H*11/12
rows_y=[.604,.308]
fig.text(.010,.941,"C",fontsize=17.5,weight="bold")
for x,head in zip(cols,["Observed cells","Smoothed grid","Increasing time",r"$\it{Fosl1}$ knockout","Randomized control"]):
    fig.text(x+width/2,.915,head,ha="center",weight="bold",fontsize=9)
fig.text(.91,.915,"TRALI mouse alignment",ha="center",weight="bold",fontsize=8.9)
seq_norm=colors.Normalize(0,1)
tables={}
primary_rows=[]
for tag,sub,xlim,ylim in specs:
    folder=root/"分析结果"/sub
    m=pd.read_csv(folder/"observed_cell_pseudotime.csv.gz")
    g=pd.read_csv(folder/"grid_reference_and_fields.csv")
    b=pd.read_csv(folder/"baseline_by_mouse.csv").query("method=='scTour_0'").copy()
    primary=pd.read_csv(folder/"alignment_by_mouse.csv").query("method=='scTour_0' and smoothing=='main'").copy()
    b["subtype"]=tag;primary["subtype"]=tag
    primary_rows.append(primary)
    mask=g.visible & g.ref_supported
    dot=1000*(g.dx*g.gradient_x+g.dy*g.gradient_y)
    null=1000*(g.random_dx*g.gradient_x+g.random_dy*g.gradient_y)
    assert not (g.visible & ~g.ref_supported).any()
    tables[tag]=(m,g,mask,dot,null,primary,b)
maximum=max(np.abs(np.r_[dot[mask],null[mask]]).max() for _,_,mask,dot,null,_,_ in tables.values())
signed_norm=colors.Normalize(-maximum,maximum)
placements={}
def tidy(a):
    a.spines[["top","left","right"]].set_visible(False)
    a.tick_params(axis="y",length=0,pad=3)
    a.tick_params(axis="x",length=2,pad=2,labelsize=7.4)
combined=pd.concat(primary_rows,ignore_index=True)
value=1000*combined[["KO_dot","random_dot"]].to_numpy()
low=min(0,np.floor(value.min()/10)*10);high=max(0,np.ceil(value.max()/10)*10)
for rid,(tag,sub,xlim,ylim) in enumerate(specs):
    m,g,mask,dot,null,primary,b=tables[tag]
    y=rows_y[rid]
    fig.text(.022,y+height/2,tag,rotation=90,ha="center",va="center",weight="bold",fontsize=10.5)
    placements[tag]=[]
    for j,x in enumerate(cols):
        box=[x,y,width,height]
        a=fig.add_axes(box)
        a.set(xlim=xlim,ylim=ylim,aspect="equal")
        a.set_axis_off()
        if j==0:
            shuffled=m.sample(frac=1,random_state=41)
            a.scatter(shuffled.UMAP_1,shuffled.UMAP_2,c=shuffled.scTour_0,cmap=seq,norm=seq_norm,
                      s=2,linewidths=0,rasterized=True)
        elif j in (1,2):
            valid=g.ref_supported
            a.scatter(g.x[valid],g.y[valid],c=g.pseudotime[valid],cmap=seq,norm=seq_norm,s=6.7,marker="s",linewidths=0)
            if j==2:
                ix=pd.factorize(g.x,sort=True)[0];iy=pd.factorize(g.y,sort=True)[0]
                length=np.hypot(g.gradient_x,g.gradient_y)
                show=valid & (ix%2==0) & (iy%2==0) & (length>1e-6)
                a.quiver(g.x[show],g.y[show],g.gradient_x[show]/length[show],g.gradient_y[show]/length[show],
                         color=ink,angles="xy",scale_units="xy",scale=1.85,pivot="middle",width=.0063)
        else:
            z=dot if j==3 else null
            a.scatter(g.x[mask],g.y[mask],c=z[mask],cmap=div,norm=signed_norm,s=6.7,marker="s",linewidths=0)
            placements[tag].append(box)
    a=fig.add_axes([.844,y+.042,.142,height-.055]);tidy(a)
    a.axvline(0,color="#89989E",ls=(0,(2,2)),lw=.65)
    primary=primary.sort_values("mouse")
    for yy,(_,r) in zip([3,2,1,0],primary.iterrows()):
        z=np.array([r.random_dot,r.KO_dot])*1000
        a.plot(z,[yy,yy],c="#AFBBBB",lw=.85)
        for x,s,c in zip(z,["s","o"],["#758B94","#9B5357"]):
            a.scatter(x,yy,s=20,marker=s,c=c,edgecolors="white",linewidths=.4,zorder=3)
    a.set(yticks=[3,2,1,0],yticklabels=primary.mouse,ylim=(-.45,3.45),xlim=(low-1.5,high+1.5))
    a.tick_params(axis="y",labelsize=7.5)
    if rid==1:a.set_xlabel("Signed projection (×1,000)",fontsize=7.6,labelpad=2)
for left,w,cmap,norm,label,ticks in [(.10,.26,seq,seq_norm,"scTour pseudotime",[0,.5,1]),
                                    (.535,.213,div,signed_norm,"Displacement · gradient (×1,000)",None)]:
    cax=fig.add_axes([left,.267,w,.012])
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal",ticks=ticks)
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=6.9,length=2,pad=2)
    cb.set_label(label,fontsize=7.5,labelpad=2)
group_cols={"Control":"#899DA3","LPS":"#5F93AC","TRALI":"#B77163"}
fig.text(.031,.133,"Observed mouse\nmedians",fontsize=8.2,weight="bold",va="center")
for j,(tag,_,_,_) in enumerate(specs):
    b=tables[tag][-1]
    a=fig.add_axes([.233+j*.34,.057,.236,.105]);tidy(a)
    for y,(group,col) in enumerate(group_cols.items()):
        v=b[b.group.eq(group)].sort_values("mouse")
        a.scatter(v["median"],2-y+np.linspace(-.19,.19,len(v)),s=17,c=col,edgecolors="white",lw=.38)
    a.set(yticks=[2,1,0],yticklabels=list(group_cols),xticks=[0,.5,1],xlim=(-.025,1.025),ylim=(-.42,2.42))
    a.tick_params(axis="y",labelsize=7)
    a.set_title(tag,fontsize=8,weight="bold",pad=1)
    a.set_xlabel("scTour pseudotime",fontsize=7,labelpad=1)
for j,(label,mark,col) in enumerate([("KO","o","#9B5357"),("Random","s","#758B94")]):
    fig.add_artist(plt.Line2D([.882],[.181-j*.041],transform=fig.transFigure,marker=mark,color=col,ls="",ms=4.2))
    fig.text(.893,.181-j*.041,label,va="center",fontsize=7.6)
fig.savefig(dest/"base_layout.pdf");plt.close(fig)
(dest/"geometry.json").write_text(json.dumps({"by_subtype":placements}))
