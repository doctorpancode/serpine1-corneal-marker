# -*- coding: utf-8 -*-
import io

from docx import Document

P = r"D:\A生信分析\虚拟敲除\文章\草稿1.4_带图.docx"
doc = Document(P)
out = []
out.append(f"paragraphs: {len(doc.paragraphs)}  tables: {len(doc.tables)}")
for i, t in enumerate(doc.tables):
    hdr = " | ".join(c.text for c in t.rows[0].cells)
    out.append(f"TABLE {i}: {len(t.rows)}x{len(t.columns)}  header={hdr[:90]}")
out.append("=" * 70)
for i, p in enumerate(doc.paragraphs):
    t = p.text.strip()
    if not t:
        continue
    out.append(f"[{i:03d}] {t}")
io.open(r"D:\A生信分析\虚拟敲除\_draft14_dump.txt", "w", encoding="utf-8").write("\n".join(out))
print("ok", len(out))
for line in out[:8]:
    print(line.encode("ascii", "replace").decode("ascii"))