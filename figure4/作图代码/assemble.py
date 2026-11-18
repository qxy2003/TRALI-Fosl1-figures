"""Export the five formal panels and assemble the approved trajectory layout."""
from pathlib import Path
import json
from PIL import Image
import pymupdf as pdf

root=Path(__file__).resolve().parents[1]
out=root/"figure"
temp=root/"作图代码/intermediate"
trajectory=temp/"4C"
def rect(page,box):
    x,y,w,h=box
    W,H=page.rect.width,page.rect.height
    return pdf.Rect(x*W,(1-y-h)*H,(x+w)*W,(1-y)*H)

with pdf.open(temp/"base.pdf") as src:
    for tag,box in json.loads((temp/"geometry.json").read_text())["panels"].items():
        clip=rect(src[0],box)
        with pdf.open() as panel:
            page=panel.new_page(width=clip.width,height=clip.height)
            page.show_pdf_page(page.rect,src,0,clip=clip)
            panel.save(out/f"Figure4{tag}.pdf",garbage=4,deflate=True)
with pdf.open(trajectory/"base_layout.pdf") as panel:
    boxes=json.loads((trajectory/"geometry.json").read_text())["by_subtype"]
    for subtype,folder in [("gCap","4C_轨迹"),("aCap","4C_aCap")]:
        with pdf.open(trajectory/f"{subtype}_arrows.pdf") as arrows:
            assert len(arrows)==len(boxes[subtype])==2
            for index,box in enumerate(boxes[subtype]):
                dest=rect(panel[0],box)
                assert abs(dest.width/dest.height-arrows[index].rect.width/arrows[index].rect.height)<.001
                panel[0].show_pdf_page(dest,arrows,index)
    panel.save(out/"Figure4C.pdf",garbage=4,deflate=True)
panels={}
try:
    for tag in "ABCDE":
        panels[tag]=pdf.open(out/f"Figure4{tag}.pdf")
    W=11.1*72
    top=panels["A"][0].rect.height
    bottom=panels["E"][0].rect.height
    width_c=.98*W
    c=panels["C"][0].rect
    height_c=width_c*c.height/c.width
    gap=6
    with pdf.open() as combined:
        page=combined.new_page(width=W,height=top+height_c+bottom+2*gap)
        for tag,x,y in [("A",.011*W,0),("B",.494*W,0),("D",.011*W,top+height_c+2*gap),
                        ("E",.433*W,top+height_c+2*gap)]:
            s=panels[tag][0].rect
            page.show_pdf_page(pdf.Rect(x,y,x+s.width,y+s.height),panels[tag],0)
        x,y=.01*W,top+gap
        page.show_pdf_page(pdf.Rect(x,y,x+width_c,y+height_c),panels["C"],0)
        combined.save(out/"Figure4.pdf",garbage=4,deflate=True)
finally:
    for doc in panels.values():
        doc.close()
with pdf.open(out/"Figure4.pdf") as combined:
    combined[0].get_pixmap(dpi=180,alpha=False).save(out/"Figure4.png")
    pix=combined[0].get_pixmap(dpi=600,alpha=False)
    Image.frombytes("RGB",(pix.width,pix.height),pix.samples).save(out/"Figure4.tiff",compression="tiff_lzw",dpi=(600,600))
with pdf.open(out/"Figure4E.pdf") as panel:
    panel[0].get_pixmap(dpi=180,alpha=False).save(out/"Figure4E.png")
