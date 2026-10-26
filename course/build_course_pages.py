"""
build_course_pages.py — Render the course's Python demos to docs/course/*.html.

Executes each module's Python demo and writes its rendered HTML (code, console
output and inline figures) into docs/course/, which GitHub Pages serves. This is
the Python counterpart of course/build_course_pages.m.

Requires: jupytext, nbconvert, nbformat, ipykernel (all in the dev environment).
Run from anywhere:  python course/build_course_pages.py
Only some modules:  python course/build_course_pages.py 00 01

This file lives in course/, which toolbox.ignore excludes from the packaged toolbox.
"""
import os
import sys
import jupytext
from nbformat.v4 import new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor, TagRemovePreprocessor
from nbconvert import HTMLExporter
from traitlets.config import Config

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repo root
COURSE = os.path.join(ROOT, "course")
OUT = os.path.join(ROOT, "docs", "course")
os.makedirs(OUT, exist_ok=True)

DEMOS = [
    ("00_Escribir_Modelos",             "demo_00_escribir_modelos", "modulo-00-python.html"),
    ("01_Introduccion_Modelos",         "demo_01_introduccion",  "modulo-01-python.html"),
    ("02_Muestreo",                     "demo_02_muestreo",      "modulo-02-python.html"),
    ("03_Incertidumbre_MonteCarlo",     "demo_03_incertidumbre", "modulo-03-python.html"),
    ("04_Analisis_Sensibilidad",        "demo_04_sensibilidad",  "modulo-04-python.html"),
    ("05_Estimacion_Identificabilidad", "demo_05_estimacion",    "modulo-05-python.html"),
]

# Hide the injected setup cell (plumbing) from the rendered page.
_cfg = Config()
_cfg.TagRemovePreprocessor.remove_cell_tags = ("hide_setup",)
_cfg.TagRemovePreprocessor.enabled = True
html_exporter = HTMLExporter(config=_cfg)
html_exporter.register_preprocessor(TagRemovePreprocessor(config=_cfg), enabled=True)

if len(sys.argv) > 1:
    DEMOS = [d for d in DEMOS if d[0][:2] in sys.argv[1:]]

for mod, demo, out_name in DEMOS:
    pydir = os.path.join(COURSE, mod, "python")
    src = os.path.join(pydir, demo + ".py")
    print(f"building {out_name} ...")

    nb = jupytext.read(src)
    # The demos resolve sibling modules via __file__ (absent in a kernel) and save
    # every figure with savefig() instead of displaying it. We inject __file__ and
    # patch savefig so each saved figure is ALSO shown inline exactly once (display
    # then close), which reliably renders every demo figure in the page.
    setup = (
        '%matplotlib inline\n'
        f'__file__ = r"{src}"\n'
        'import matplotlib.pyplot as _plt\n'
        'from IPython.display import display as _display\n'
        '_orig_savefig = _plt.Figure.savefig\n'
        'def _savefig_and_show(self, *a, **k):\n'
        '    r = _orig_savefig(self, *a, **k)\n'
        '    _display(self); _plt.close(self)\n'
        '    return r\n'
        '_plt.Figure.savefig = _savefig_and_show\n'
    )
    setup_cell = new_code_cell(setup)
    setup_cell.metadata["tags"] = ["hide_setup"]
    nb.cells.insert(0, setup_cell)

    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": pydir}})   # run with cwd = the demo's folder

    body, _ = html_exporter.from_notebook_node(nb)
    with open(os.path.join(OUT, out_name), "w", encoding="utf-8") as f:
        f.write(body)
    print(f"   done ({len(body)//1024} KB)")

print("Listo.")
