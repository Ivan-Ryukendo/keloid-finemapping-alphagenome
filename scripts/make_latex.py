"""Markdown paper -> two-column LaTeX (Elsevier CAS-like look) -> PDF with Tectonic.
Called by build_paper.py --pdf. Needs pandoc and tectonic (see build_paper.py for lookup)."""
import os, re, shutil, subprocess
from pathlib import Path

CITE = re.compile(r"\[(\d+(?:\s*[,\u2013-]\s*\d+)*)\](?!\()")
UNI = {"\u2192": r"$\rightarrow$", "\u2265": r"$\geq$", "\u2264": r"$\leq$", "\u2212": r"$-$", "\u00d7": r"$\times$",
       "\u00b1": r"$\pm$", "\u2248": r"$\approx$", "\u2260": r"$\neq$", "\u2011": "-", "\u03b2": r"$\beta$", "\u03b1": r"$\alpha$",
       "\u207b": r"\textsuperscript{$-$}", "\u2078": r"\textsuperscript{8}"}


def cite_keys(spec):
    keys = []
    for part in re.split(r"\s*,\s*", spec):
        m = re.match(r"(\d+)\s*[\u2013-]\s*(\d+)$", part)
        keys += [f"k{i}" for i in range(int(m[1]), int(m[2]) + 1)] if m else [f"k{int(part)}"]
    return ",".join(keys)


def cites(md):
    return CITE.sub(lambda m: "\\cite{" + cite_keys(m[1]) + "}", md)


def esc(t):
    t = t.replace("\\", "\x00")
    for a, b in [("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")]:
        t = t.replace(a, b)
    return t.replace("\x00", r"\textbackslash{}")


def inline(t):
    """Tiny markdown-inline -> LaTeX for table cells."""
    out = []
    for k, seg in enumerate(CITE.split(t)):  # odd items are the captured citation numbers
        if k % 2:
            out.append("\\cite{" + cite_keys(seg) + "}"); continue
        for tok in re.split(r"(\*\*.+?\*\*|\*.+?\*|`.+?`)", seg):
            if tok.startswith("**") and tok.endswith("**") and len(tok) > 4: out.append(r"\textbf{" + esc(tok[2:-2]) + "}")
            elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2: out.append(r"\emph{" + esc(tok[1:-1]) + "}")
            elif tok.startswith("`") and tok.endswith("`"): out.append(r"\texttt{" + esc(tok[1:-1]) + "}")
            else: out.append(esc(tok).replace("+", "+\\hspace{0pt}"))  # allow breaks in ATAC+DNASE+... lists
    return "".join(out).replace("\\|", "|")


def split_row(line):
    return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", line.strip())[1:-1]]


def table_rows(lines):
    return [split_row(l) for l in lines if not re.fullmatch(r"\|[\s:|-]+\|", l.strip())]


def tabularx_env(lines, spec, caption, size=r"\footnotesize"):
    rows = table_rows(lines)
    head = " & ".join(r"\textbf{" + inline(c) + "}" for c in rows[0])
    body = "\n".join(" & ".join(inline(c) for c in r) + r" \\" for r in rows[1:])
    return ("\\begin{table*}[t]\n\\centering\n\\tabcap{" + caption[0] + "}{" + caption[1] + "}\n{" + size + "\\setlength{\\tabcolsep}{3pt}\\renewcommand{\\arraystretch}{1.18}\n"
            "\\begin{tabularx}{\\textwidth}{" + spec + "}\n\\toprule\n" + head + " \\\\\n\\midrule\n" + body + "\n\\bottomrule\n\\end{tabularx}}\n\\end{table*}")


def longtable_env(lines, spec):
    rows = table_rows(lines)
    head = " & ".join(r"\textbf{" + inline(c) + "}" for c in rows[0]) + " \\\\\n"
    body = "\n".join(" & ".join(inline(c) for c in r) + r" \\" for r in rows[1:])
    return ("{\\footnotesize\\setlength{\\tabcolsep}{4pt}\\renewcommand{\\arraystretch}{1.12}\n\\begin{longtable}{" + spec + "}\n\\toprule\n" + head +
            "\\midrule\n\\endfirsthead\n\\toprule\n" + head + "\\midrule\n\\endhead\n\\bottomrule\n\\endfoot\n" + body + "\n\\end{longtable}}")


def figure_env(alt, path):
    m = re.match(r"(Figure \d+[a-z]?\.)\s*(.*)", alt, re.S)
    label, text = (m[1], m[2]) if m else ("", alt)
    name = Path(path).stem
    narrow = "heatmap" in name  # tall single-column figure
    env, width = ("figure", r"\columnwidth") if narrow else ("figure*", r"\textwidth")
    return (f"\\begin{{{env}}}[t]\n\\centering\n\\includegraphics[width={width}]{{{name}.pdf}}\n"
            f"\\figcap{{{esc(label)}}}{{{inline(text)}}}\n\\end{{{env}}}")


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, encoding="utf8", **kw)


def build(md, paper_dir, pandoc, tectonic, cache_dir):
    paper_dir = Path(paper_dir); out_dir = paper_dir / "latex"; out_dir.mkdir(exist_ok=True)
    pd_ = lambda s, extra=(): run([str(pandoc), "-f", "markdown-auto_identifiers", "-t", "latex", "--wrap=preserve", *extra], input=s).stdout.strip()

    parts = re.split(r"^## ", md, flags=re.M)
    head_md, secs = parts[0], {p.split("\n", 1)[0].strip(): p.split("\n", 1)[1] for p in parts[1:]}
    title = re.match(r"# (.+)", head_md)[1].strip()
    lines = [l for l in head_md.splitlines()[1:] if l.strip()]
    author_line, note_line = lines[0], lines[1]
    status = re.sub(r"^\*\*Status\.\*\*\s*", "", " ".join(l.lstrip("> ").strip() for l in lines[2:]))
    name = re.match(r"\*\*(.+?)\*\*", author_line)[1]
    email = re.search(r"\(([^)]+@[^)]+)\)", author_line)[1]
    abstract = pd_(cites(secs["Abstract"]))

    body_keys = [k for k in secs if re.match(r"\d+\. ", k) and not k.startswith("10.")]
    body_md = "\n\n".join("## " + k + "\n" + secs[k] for k in body_keys)
    body_md = re.sub(r"^(#{2,3}) \d+(?:\.\d+)*\.? ", r"\1 ", body_md, flags=re.M)
    body_md = re.sub(r"^### \d+(?:\.\d+)*\.? ", "### ", body_md, flags=re.M)

    tspec = [(r"@{}rl>{\raggedright\arraybackslash}p{2.7cm}rr>{\raggedright\arraybackslash}Xcll@{}", ("Table 1.", "Top ten loci by rank_score.")),
             (r"@{}c>{\raggedright\arraybackslash}p{3.2cm}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}X@{}", ("Table 2.", "Proposed first experiments, strongest hypotheses first."))]
    floats = []
    def take_tables(m):
        i = len(floats); spec, cap = tspec[i]
        floats.append(tabularx_env(m[0].strip().splitlines(), spec, (cap[0], esc(cap[1]).replace(r"rank\_score", r"\texttt{rank\_score}")),
                                   r"\scriptsize" if i == 0 else r"\footnotesize"))
        return f"\nXXFLOAT{i}XX\n"
    body_md = re.sub(r"(?:^\|.*\n?)+", take_tables, body_md, flags=re.M)
    def take_figs(m):
        i = len(floats); floats.append(figure_env(m[1], m[2])); return f"\nXXFLOAT{i}XX\n"
    body_md = re.sub(r"^!\[(.*?)\]\((.*?)\)\s*$", take_figs, body_md, flags=re.M)
    body = pd_(cites(body_md), ["--shift-heading-level-by=-1"])
    for i, f in enumerate(floats):
        assert f"XXFLOAT{i}XX" in body, i
        body = body.replace(f"XXFLOAT{i}XX", f)

    refs = [re.sub(r"^\d+\. ", "", l) for l in secs["10. References"].splitlines() if re.match(r"\d+\. ", l)]
    assert len(refs) == 25, len(refs)
    ref_tex = pd_("\n\n".join(refs)).split("\n\n")
    assert len(ref_tex) == len(refs)
    def doi(m): return r"\href{https://doi.org/" + m[1].replace("\\_", "_").replace("\\%", "%").replace("\\#", "#") + "}{doi:" + m[1] + "}"
    ref_tex = [re.sub(r"doi:(10\.[^\s]*?)(?=\.?(?:\s|$|\)))", doi, r) for r in ref_tex]
    bib = "\\begin{thebibliography}{99}\n" + "\n".join(f"\\bibitem{{k{i}}} {r}" for i, r in enumerate(ref_tex, 1)) + "\n\\end{thebibliography}"
    verif = [p for p in secs["10. References"].split("\n\n") if p.startswith("Reference verification")]
    verif_tex = "{\\footnotesize " + pd_(verif[0]) + "\\par}" if verif else ""

    sup_key = next(k for k in secs if k.startswith("Supplementary Table S1"))
    sup_intro, sup_tab = secs[sup_key].split("{{S1}}")[0] if "{{S1}}" in secs[sup_key] else secs[sup_key], None
    sup_lines = [l for l in secs[sup_key].splitlines() if l.startswith("|")]
    sup_intro = "\n".join(l for l in secs[sup_key].splitlines() if not l.startswith("|"))
    sup = ("\\clearpage\\onecolumn\n\\section*{" + esc(sup_key) + "}\n" + pd_(sup_intro) + "\n\\medskip\n" +
           longtable_env(sup_lines, r"@{}l>{\raggedright\arraybackslash}p{5.2cm}rcccl@{}"))

    uni = "\n".join(f"\\newunicodechar{{{k}}}{{{v}}}" for k, v in UNI.items())
    tex = (Path(__file__).with_name("paper_template.tex").read_text(encoding="utf8")
           .replace("%%UNICODE%%", uni).replace("%%TITLE%%", pd_(title).replace("\n", " ")).replace("%%AUTHOR%%", esc(name))
           .replace("%%EMAIL%%", esc(email)).replace("%%NOTE%%", pd_(note_line)).replace("%%STATUS%%", pd_(status))
           .replace("%%ABSTRACT%%", abstract).replace("%%BODY%%", body).replace("%%BIB%%", bib)
           .replace("%%VERIF%%", verif_tex).replace("%%SUPP%%", sup))
    (out_dir / "paper.tex").write_text(tex, encoding="utf8")

    for f in (paper_dir.parent / "figures").glob("*.pdf"): shutil.copy(f, out_dir / f.name)
    env = dict(os.environ, TECTONIC_CACHE_DIR=str(cache_dir))
    r = subprocess.run([str(tectonic), "-X", "compile", "paper.tex", "--keep-logs", "--print"], cwd=out_dir, env=env, capture_output=True, text=True, encoding="utf8")
    (out_dir / "build.log").write_text(r.stdout + r.stderr, encoding="utf8")
    if r.returncode: raise SystemExit("tectonic failed; see " + str(out_dir / "build.log") + "\n" + (r.stdout + r.stderr)[-1500:])
    return out_dir / "paper.pdf"
