"""Build the manuscript inside the official Cancers (MDPI) Word template (MDPI styles, headers, line numbers, history box).
Derived from the IJMS builder of version 1.0.1; adds the Simple Summary and reads abbreviations from the source.
Run: python3 -I build_cancers_docx.py <cancers-template.dot> <manuscript_source.md> <figures_dir> <out.docx>"""
import copy, os, re, shutil, sys, tempfile, zipfile
from docx import Document
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import qn, nsdecls
from docx.shared import Cm, Pt

TPL, SRC, FIG, OUT = sys.argv[1:5]

# --- template (.dot is OOXML template) -> document ---
tmp = tempfile.mkdtemp(); fixed = os.path.join(tmp, 'tpl.docx')
with zipfile.ZipFile(TPL) as zi, zipfile.ZipFile(fixed, 'w', zipfile.ZIP_DEFLATED) as zo:
    for it in zi.infolist():
        data = zi.read(it.filename)
        if it.filename == '[Content_Types].xml':
            data = data.replace(b'wordprocessingml.template.main+xml', b'wordprocessingml.document.main+xml')
        zo.writestr(it, data)
doc = Document(fixed); body = doc.element.body
kids = list(body)

def style_of(el):
    ps = el.find(qn('w:pPr'))
    s = ps.find(qn('w:pStyle')) if ps is not None else None
    return s.get(qn('w:val')) if s is not None else ''

def clear_runs(p_el):
    for ch in list(p_el):
        if ch.tag != qn('w:pPr'): p_el.remove(ch)

def add_run(p_el, text, bold=False, italic=False, sup=False, size=None):
    """Run properties in the schema order Word requires: b, bCs, i, iCs, sz, szCs, vertAlign."""
    r = OxmlElement('w:r'); rpr = OxmlElement('w:rPr')
    if bold: rpr.append(OxmlElement('w:b')); rpr.append(OxmlElement('w:bCs'))
    if italic: rpr.append(OxmlElement('w:i')); rpr.append(OxmlElement('w:iCs'))
    if size:
        for t in ('w:sz', 'w:szCs'):
            s = OxmlElement(t); s.set(qn('w:val'), str(size)); rpr.append(s)
    if sup:
        v = OxmlElement('w:vertAlign'); v.set(qn('w:val'), 'superscript'); rpr.append(v)
    if len(rpr): r.append(rpr)
    t = OxmlElement('w:t'); t.set(qn('xml:space'), 'preserve'); t.text = text; r.append(t)
    p_el.append(r)

def md_runs(p_el, text, size=None, base_bold=False):
    for tok in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^)', text):
        if not tok: continue
        b, i, sup = base_bold, False, False
        if tok.startswith('**'): tok, b = tok[2:-2], True
        elif tok.startswith('*') and len(tok) > 1: tok, i = tok[1:-1], True
        elif tok.startswith('^'): tok, sup = tok[1:-1], True
        add_run(p_el, tok.replace('\\*', '*'), b, i, sup, size)

def new_p(style, keep_next=False, keep_lines=False):
    p = doc.add_paragraph()
    ppr = p._p.get_or_add_pPr(); ps = OxmlElement('w:pStyle'); ps.set(qn('w:val'), style); ppr.insert(0, ps)
    if keep_next: ppr.append(OxmlElement('w:keepNext'))
    if keep_lines: ppr.append(OxmlElement('w:keepLines'))
    return p._p

# --- read source ---
lines = open(SRC, encoding='utf-8').read().split('\n')
title = next(l[2:] for l in lines if l.startswith('# '))
authors = next(l[8:].strip() for l in lines if l.startswith('AUTHORS:'))
affils = [l[6:].strip().split('|', 1) for l in lines if l.startswith('AFFIL:')]
abstract = lines[lines.index('## Abstract') + 1]
simple = lines[lines.index('## Simple Summary') + 1]
ABBR = [tuple(l[5:].strip().split('|', 1)) for l in lines if l.startswith('ABBR:')]
keywords = next(l for l in lines if l.startswith('**Keywords:**')).replace('**Keywords:**', '').strip()

# --- front matter in template elements ---
art, tit, aut = kids[0], kids[1], kids[2]
clear_runs(art); add_run(art, 'Article')
clear_runs(tit); add_run(tit, title)
clear_runs(aut)
for k, part in enumerate(re.split(r'\^([^^]+)\^', authors.replace('\\*', '*'))):
    if k % 2 == 0:
        if part: add_run(aut, part)
    else:
        add_run(aut, part.rstrip('*'), sup=True)
        if part.endswith('*'): add_run(aut, '*')
aff_tpl = [k for k in kids if style_of(k) == 'MDPI16affiliation']
proto = copy.deepcopy(aff_tpl[0])
for k in aff_tpl[1:]: body.remove(k)
anchor = aff_tpl[0]
for n, (num, txt) in enumerate(affils):
    p = anchor if n == 0 else copy.deepcopy(proto)
    clear_runs(p); num = num.replace('\\*', '*')
    add_run(p, num, sup=num != '*', bold=num == '*')
    r = OxmlElement('w:r'); r.append(OxmlElement('w:tab')); p.append(r); add_run(p, txt.strip())
    if n: prev.addnext(p)
    prev = p
abs_ps = [k for k in kids if style_of(k) == 'MDPI17abstract']
assert len(abs_ps) == 4, len(abs_ps)  # Simple Summary label, text, Abstract label, text
clear_runs(abs_ps[1]); md_runs(abs_ps[1], simple)
clear_runs(abs_ps[3]); md_runs(abs_ps[3], abstract)
kw = next(k for k in kids if style_of(k) == 'MDPI18keywords')
clear_runs(kw); add_run(kw, 'Keywords: ', bold=True); add_run(kw, keywords)

# --- drop template body (after the MDPI19line) except the publisher's note and sectPr ---
kids = list(body)
start = next(i for i, k in enumerate(kids) if style_of(k) == 'MDPI19line') + 1
notes = next(k for k in kids if style_of(k) == 'MDPI63notes')
sectPr = kids[-1]
for k in kids[start:-1]: body.remove(k)

# --- tables & figures ---
from xml.sax.saxutils import escape as _esc
W = nsdecls('w')

def _runs_xml(text, bold=False, size=None):
    """Runs for table cells, honouring **bold**, *italic*, ^sup^ markup. Mirrors template cell runs."""
    out = []
    for tok in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^)', text):
        if not tok: continue
        b, i, sup = bold, False, False
        if tok.startswith('**'): tok, b = tok[2:-2], True
        elif tok.startswith('*') and len(tok) > 1: tok, i = tok[1:-1], True
        elif tok.startswith('^'): tok, sup = tok[1:-1], True
        rpr = ('<w:b/><w:bCs/>' if b else '') + ('<w:i/><w:iCs/>' if i else '') + \
              (f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>' if size else '') + ('<w:vertAlign w:val="superscript"/>' if sup else '')
        out.append(f'<w:r>{"<w:rPr>" + rpr + "</w:rPr>" if rpr else ""}<w:t xml:space="preserve">{_esc(tok.replace(chr(92) + "*", "*"))}</w:t></w:r>')
    return ''.join(out)

def _border(top, bottom):
    def e(edge, v): return f'<w:{edge} w:val="single" w:sz="{v}" w:space="0" w:color="auto"/>' if v else f'<w:{edge} w:val="nil"/>'
    return f'<w:tcBorders>{e("top", top)}<w:left w:val="nil"/>{e("bottom", bottom)}<w:right w:val="nil"/></w:tcBorders>'

def make_table(rows, widths, kind, split_ci=False):
    """Tables built from the template's own definitions:
    kind 'column' = template Table 1 (text-column width, indented 4.6 cm, line spacing 12 pt);
    kind 'wide'   = template Table 2 (full text width, centred, centred rows);
    kind 'abbr'   = template Abbreviations table (no borders, 9 pt, left-aligned)."""
    n = len(rows)
    if kind == 'column':
        tblpr = (f'<w:tblW w:w="{sum(widths)}" w:type="dxa"/><w:tblInd w:w="2608" w:type="dxa"/>'
                 '<w:tblBorders><w:top w:val="single" w:sz="8" w:space="0" w:color="auto"/><w:bottom w:val="single" w:sz="8" w:space="0" w:color="auto"/></w:tblBorders>')
    elif kind == 'wide':
        tblpr = (f'<w:tblW w:w="{sum(widths)}" w:type="dxa"/><w:jc w:val="center"/>'
                 '<w:tblBorders><w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/></w:tblBorders>')
    else:
        tblpr = (f'<w:tblStyle w:val="TableGrid"/><w:tblW w:w="{sum(widths)}" w:type="dxa"/><w:tblInd w:w="2608" w:type="dxa"/>'
                 '<w:tblBorders><w:top w:val="none" w:sz="0" w:space="0" w:color="auto"/><w:left w:val="none" w:sz="0" w:space="0" w:color="auto"/>'
                 '<w:bottom w:val="none" w:sz="0" w:space="0" w:color="auto"/><w:right w:val="none" w:sz="0" w:space="0" w:color="auto"/>'
                 '<w:insideH w:val="none" w:sz="0" w:space="0" w:color="auto"/><w:insideV w:val="none" w:sz="0" w:space="0" w:color="auto"/></w:tblBorders>')
    tblpr += ('<w:tblLayout w:type="fixed"/><w:tblCellMar><w:left w:w="0" w:type="dxa"/><w:right w:w="0" w:type="dxa"/></w:tblCellMar>'
              '<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/>')
    xml = [f'<w:tbl {W}><w:tblPr>{tblpr}</w:tblPr><w:tblGrid>' + ''.join(f'<w:gridCol w:w="{w}"/>' for w in widths) + '</w:tblGrid>']
    for i, row in enumerate(rows):
        hdr = kind != 'abbr' and i == 0
        xml.append('<w:tr>' + ('<w:trPr><w:jc w:val="center"/></w:trPr>' if kind == 'wide' else ''))
        for j, txt in enumerate(row):
            if kind == 'abbr':
                borders = ''
            elif hdr:
                borders = _border('8', '4')
            else:
                borders = _border(None, '8' if i == n - 1 else None)
            keep = '<w:keepNext/>' if (i < n - 1 and kind != 'abbr') else ''   # keeps data tables on one page; the abbreviations list may flow
            spacing = '' if kind == 'wide' else '<w:spacing w:line="240" w:lineRule="auto"/>'
            jc = '<w:jc w:val="left"/>' if kind == 'abbr' else ''
            size = 18 if kind in ('abbr', 'wide') else None   # 9 pt: template abbreviations size; wide tables need it to avoid wrapping
            m = re.match(r'^(.*\S) ([\(\[][^\)\]]*,[^\)\]]*[\)\]])$', txt) if (split_ci and not hdr) else None
            runs = (_runs_xml(m.group(1), size=size) + '<w:r><w:br/></w:r>' + _runs_xml(m.group(2), size=size)) if m \
                   else _runs_xml(txt, bold=hdr, size=size)
            xml.append(f'<w:tc><w:tcPr><w:tcW w:w="{widths[j]}" w:type="dxa"/>{borders}<w:vAlign w:val="center"/></w:tcPr>'
                       f'<w:p><w:pPr><w:pStyle w:val="MDPI42tablebody"/>{keep}{spacing}{jc}</w:pPr>{runs}</w:p></w:tc>')
        xml.append('</w:tr>')
    xml.append('</w:tbl>')
    tbl = parse_xml(''.join(xml))
    sectPr.addprevious(tbl)
    return tbl

def add_table(caption, rows_md, note):
    p = new_p('MDPI41tablecaption', keep_next=True, keep_lines=True)
    m = re.match(r'\*\*(Table \d+\.)\*\*\s*(.*)', caption)
    add_run(p, m.group(1) + ' ', bold=True); md_runs(p, m.group(2))
    rows = [[c.strip() for c in r.strip('|').split('|')] for r in rows_md if not re.match(r'^\|\s*-', r)]
    ncol = len(rows[0])
    if ncol <= 3:   # fits the text column: template Table 1 layout (width 7857)
        widths = [5500, 2357] if ncol == 2 else [2619, 2619, 2619]
        make_table(rows, widths, 'column')
    else:           # template Table 2 layout (width 10465)
        first = 2465 if ncol >= 7 else 2700
        ncw = 560 if rows[0][1].strip().lower() == 'n' else None
        k = ncol - (2 if ncw else 1); rest = (10465 - first - (ncw or 0)) // k
        widths = [first] + ([ncw] if ncw else []) + [rest] * k; widths[-1] += 10465 - sum(widths)
        make_table(rows, widths, 'wide', split_ci=ncol >= 5)
    if note:
        q = new_p('MDPI43tablefooter', keep_lines=True); md_runs(q, note)

def add_figure(num, caption):
    p = new_p('MDPI52figure', keep_next=True)
    ppr = p.find(qn('w:pPr'))
    ind = OxmlElement('w:ind'); ind.set(qn('w:left'), '2608'); ppr.append(ind)
    jc = OxmlElement('w:jc'); jc.set(qn('w:val'), 'left'); ppr.append(jc)
    from docx.text.paragraph import Paragraph
    from PIL import Image
    path = os.path.join(FIG, f'Figure{num}.png'); wpx = Image.open(path).size[0]
    width_cm = min(wpx / 600 * 2.54, 13.86)  # text column = 10466 - 2608 dxa = 13.86 cm
    Paragraph(p, doc._body).add_run().add_picture(path, width=Cm(width_cm))
    q = new_p('MDPI51figurecaption', keep_lines=True)
    m = re.match(r'\*\*(Figure \d+\.)\*\*\s*(.*)', caption)
    add_run(q, m.group(1) + ' ', bold=True); md_runs(q, m.group(2))

tables, figs, body_lines = {}, {}, []
i = 0
while i < len(lines):
    L = lines[i]
    m = re.match(r'\*\*Table (\d+)\.\*\*', L)
    if m:
        cap = L; rows = []; i += 1
        while i < len(lines) and not lines[i].startswith('|'): i += 1
        while i < len(lines) and lines[i].startswith('|'): rows.append(lines[i]); i += 1
        while i < len(lines) and not lines[i].strip(): i += 1
        note = None
        if i < len(lines) and not lines[i].startswith(('**', '#')): note = lines[i]; i += 1
        tables[int(m.group(1))] = (cap, rows, note); continue
    m = re.match(r'\*\*Figure (\d+)\.\*\*', L)
    if m: figs[int(m.group(1))] = L; i += 1; continue
    if L.startswith(('## Tables', '## Figure legends', 'ABBR:')): i += 1; continue
    body_lines.append(L); i += 1

placed_t, placed_f = set(), set()
def place_after(text):
    for n in sorted(tables):
        if n not in placed_t and re.search(rf'\bTable {n}\b', text): add_table(*tables[n]); placed_t.add(n)
    for n in sorted(figs):
        if n not in placed_f and re.search(rf'\bFigure {n}[A-Z,]*\b', text): add_figure(n, figs[n]); placed_f.add(n)

section = 'front'; back_done = False
for L in body_lines:
    if L.startswith('## 1. Introduction'): section = 'body'
    if section == 'front': continue
    if L.startswith('## References'):
        if not back_done:
            p = new_p('MDPI21heading1'); add_run(p, 'Abbreviations')
            p = new_p('MDPI32textnoindent'); add_run(p, 'The following abbreviations are used in this manuscript:')
            make_table([[a, b] for a, b in ABBR], [936, 6921], 'abbr')
            back_done = True
        # References heading aligned with the reference list (the list starts at the margin; the template indents the heading)
        p = new_p('MDPI21heading1', keep_next=True); p.find(qn('w:pPr')).append(parse_xml(f'<w:ind {nsdecls("w")} w:left="0"/>'))
        add_run(p, 'References'); section = 'refs'; continue
    if section == 'refs':
        m = re.match(r'\d+\. (.*)', L)
        if m: p = new_p('MDPI81references'); md_runs(p, m.group(1))
        continue
    if L.startswith('BACK:'):
        head, txt = L[5:].strip().split('|', 1)
        p = new_p('MDPI62backmatter'); add_run(p, head + ': ', bold=True); md_runs(p, txt)
        if not globals().get('_bm_started'):
            sp = OxmlElement('w:spacing'); sp.set(qn('w:before'), '240'); p.find(qn('w:pPr')).append(sp); _bm_started = True
        continue
    if L.startswith('## '):
        p = new_p('MDPI21heading1', keep_next=True); add_run(p, L[3:]); continue
    if L.startswith('### '):
        p = new_p('MDPI22heading2', keep_next=True); add_run(p, L[4:]); continue
    if L.strip():
        p = new_p('MDPI31text'); md_runs(p, L); place_after(L)

for n in sorted(set(tables) - placed_t): add_table(*tables[n])
for n in sorted(set(figs) - placed_f): add_figure(n, figs[n])
sectPr.addprevious(notes)
doc.core_properties.title = title; doc.core_properties.author = 'Hamza Abu Owida; Ayoub Alsarhan; Mohannad Banat; Dimah Al-Fraihat; Abdulmohsen Algarni'
doc.save(OUT); shutil.rmtree(tmp)
print('tables', sorted(placed_t), 'figures', sorted(placed_f))
