# -*- coding: utf-8 -*-
"""docx → PDF（Word COM）→ 逐页 PNG，用于目视质检"""
import os
import sys
import glob
import subprocess

BASE = r"D:\A生信分析\虚拟敲除\文章"
PAGES = os.path.join(BASE, "figs_v2", "out", "pages")


def to_pdf(docx, pdf):
    ps = f'''
$ErrorActionPreference="Stop"
$w = New-Object -ComObject Word.Application
$w.Visible = $false
$d = $w.Documents.Open("{docx}", $false, $true)
$d.ExportAsFixedFormat("{pdf}", 17)
$d.Close(0)
$w.Quit()
'''
    subprocess.run(["pwsh", "-NoProfile", "-Command", ps], check=True)


def pages(pdf, outdir, dpi=110):
    import fitz
    os.makedirs(outdir, exist_ok=True)
    for f in glob.glob(os.path.join(outdir, "*.png")):
        os.remove(f)
    doc = fitz.open(pdf)
    out = []
    for i, page in enumerate(doc, 1):
        p = os.path.join(outdir, f"p{i:02d}.png")
        page.get_pixmap(dpi=dpi).save(p)
        out.append(p)
    return out, len(doc)


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "草稿1.1_带图.docx")
    pdf = os.path.splitext(src)[0] + "_preview.pdf"
    to_pdf(src, pdf)
    files, n = pages(pdf, PAGES)
    print("PDF", pdf)
    print("pages", n)
    for f in files:
        print(f)