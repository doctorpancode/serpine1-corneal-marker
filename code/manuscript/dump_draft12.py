# -*- coding: utf-8 -*-
import io, sys
from docx import Document
from docx.oxml.ns import qn

d = Document(r"文章\草稿1.2_带图.docx")
out = []
n_img = 0
for p in d.paragraphs:
    if p._element.findall(".//" + qn("a:blip")):
        n_img += 1
        out.append(f"[IMG {n_img}]")
    elif p.text.strip():
        out.append(p.text.strip())
io.open("_draft12_dump.txt", "w", encoding="utf-8").write("\n".join(out))
print("paras:", len(out), "images:", n_img, "tables:", len(d.tables))
for t in d.tables:
    print("  table", len(t.rows), "x", len(t.columns))
# sanity greps
txt = "\n".join(out)
for key in ["7 138", "7138", "480 倍", "192 倍", "3.34", "1.52–13.96", "462 倍", "1 367", "1 350", "0.988", "0.488"]:
    print(f"{key!r}: {txt.count(key)}")