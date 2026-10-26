"""
mostrar_codigo.py — muestra un archivo de código completo en la página del curso.

Es el equivalente de <include>archivo.m</include> en las páginas de MATLAB: en un notebook
lo muestra con resaltado de sintaxis; en una terminal lo imprime tal cual.
"""
from pathlib import Path


def mostrar_codigo(ruta):
    texto = Path(ruta).read_text(encoding="utf-8")
    try:
        from IPython import get_ipython
        from IPython.display import HTML, display
        from pygments import highlight
        from pygments.formatters import HtmlFormatter
        from pygments.lexers import PythonLexer
        if get_ipython() is not None:
            formato = HtmlFormatter(noclasses=True, prestyles=(
                "font-family: var(--jp-code-font-family, monospace); font-size: 13px; "
                "background: #f7f7f7; padding: 8px; line-height: 1.3"))
            display(HTML(highlight(texto, PythonLexer(), formato)))
            return
    except ImportError:
        pass
    print(texto)
