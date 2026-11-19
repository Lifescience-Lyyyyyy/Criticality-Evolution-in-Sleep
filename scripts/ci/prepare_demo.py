"""Extract the first synthetic d2 section and its local functions verbatim."""

from pathlib import Path
import re
from zipfile import ZipFile


demo = Path("external/prox_crit_toolkit/demo")
with ZipFile(demo / "demo.mlx") as archive:
    document = archive.read("matlab/document.xml").decode("utf-8")

code = re.findall(r"<!\[CDATA\[(.*?)\]\]>", document, flags=re.DOTALL)
assert "b = 2;" in code[0] and "calc_db(" in code[0]
assert code[-1].startswith("function x = runAR(")
(demo / "smoke_demo.m").write_text(code[0] + "\n\n" + code[-1] + "\n", encoding="utf-8")
