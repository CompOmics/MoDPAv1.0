"""Shared test configuration.

The project is a collection of standalone analysis scripts rather than an installed
package, so this file provides the two things the tests need to reach that code:

1. `sys.path` entries for every directory that holds importable modules. Every script
   filename is a valid Python identifier, so a plain `import` reaches any script that
   guards its work behind `if __name__ == "__main__"`.
2. `load_functions()`, for the numbered stage scripts under
   `4-PTM-pairs-annotation/scripts/`, which do their work at module level and would
   therefore run real queries against absent data if imported normally.
"""
import ast
import sys
import types
from pathlib import Path

import matplotlib

# Figures are created by some of the code under test; never open a window.
matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]

SCRIPT_DIRS = [
    ROOT,
    ROOT / "1-quant-pipeline-latest",
    ROOT / "2-VAE-code",
    ROOT / "2-VAE-code" / "Sensitivity-analysis",
    ROOT / "3-pulse-silac-validation",
    ROOT / "4-PTM-pairs-annotation" / "scripts",
    ROOT / "5-pathway-ORA",
]

for _directory in SCRIPT_DIRS:
    path = str(_directory)
    if path not in sys.path:
        sys.path.insert(0, path)


def load_functions(relative_path, names):
    """Return named top-level functions from a script that does its work at import time.

    The stage scripts under `4-PTM-pairs-annotation/scripts/` have no `__main__` guard,
    so importing one would immediately run its DuckDB queries against input files that
    are not tracked in this repository. This parses the file instead and executes only
    the requested function definitions, together with the plain imports they rely on.

    Two kinds of import are dropped: `from _common import *`, which has the side effect
    of creating the `work/` directory, and any other star-import, since neither can be
    resolved without the analysis inputs. A function that needs a name from `_common`
    must therefore have it injected by the test, which is also how the module-level
    globals these functions close over (`N`, `P`, and so on) are supplied.
    """
    path = ROOT / relative_path
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    def is_wanted_import(node):
        if isinstance(node, ast.Import):
            return True
        if isinstance(node, ast.ImportFrom):
            star = any(alias.name == "*" for alias in node.names)
            return not star and node.module != "_common"
        return False

    kept = [
        node
        for node in tree.body
        if is_wanted_import(node)
        or (isinstance(node, ast.FunctionDef) and node.name in names)
    ]

    module = types.ModuleType(path.stem)
    module.__file__ = str(path)
    code = compile(ast.Module(body=kept, type_ignores=[]), str(path), "exec")
    exec(code, module.__dict__)

    missing = [name for name in names if not hasattr(module, name)]
    if missing:
        raise AttributeError(f"{relative_path} defines no {missing}")
    return module
