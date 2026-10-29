"""
Import rules, checked from the source text so no module is imported.

`lib` modules import no step package, `search` or `cli`; a step's modules import no other step's entry module.
"""

import ast
import os

import pytest

PACKAGE_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "pocketmapper"))


def find_steps():
    """
    List the step packages: subdirectories of the package with an `__init__.py`, other than `lib`.

    Returns:
        list of str: Step names, sorted.
    """
    return sorted(
        name
        for name in os.listdir(PACKAGE_DIR)
        if name != "lib" and os.path.isfile(os.path.join(PACKAGE_DIR, name, "__init__.py"))
    )


def find_modules(subdir):
    """
    List the `.py` files under one package subdirectory, skipping `__pycache__`.

    Args:
        subdir (str): Subdirectory of the package.

    Returns:
        list of str: Paths relative to the package, sorted.
    """
    found = []
    for root, dirs, files in os.walk(os.path.join(PACKAGE_DIR, subdir)):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        found += [os.path.relpath(os.path.join(root, f), PACKAGE_DIR) for f in files if f.endswith(".py")]
    return sorted(found)


def imported_names(rel_path):
    """
    Collect the module names a file imports, with `module.name` for each `from module import name`.

    Relative imports are ignored.

    Args:
        rel_path (str): File path relative to the package.

    Returns:
        set of str: Dotted names.
    """
    with open(os.path.join(PACKAGE_DIR, rel_path), encoding="utf-8") as handle:
        tree = ast.parse(handle.read(), filename=rel_path)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names.add(node.module)
            names.update(f"{node.module}.{alias.name}" for alias in node.names)
    return names


def violations(names, banned):
    """
    Pick the names that equal a banned module or sit under one.

    Args:
        names (set of str): Imported dotted names.
        banned (list of str): Banned dotted module names.

    Returns:
        list of str: The offending names, sorted.
    """
    return sorted(n for n in names if any(n == b or n.startswith(f"{b}.") for b in banned))


STEPS = find_steps()
LIB_BANNED = [f"pocketmapper.{name}" for name in STEPS] + ["pocketmapper.search", "pocketmapper.cli"]
LIB_MODULES = find_modules("lib")
STEP_MODULES = [(step, path) for step in STEPS for path in find_modules(step)]


def test_walk_finds_steps():
    """The walk finds the steps and modules, so the rule tests cannot pass on empty tables."""
    assert "fetch" in STEPS
    assert "lib" not in STEPS
    assert LIB_MODULES
    assert STEP_MODULES


@pytest.mark.parametrize("path", LIB_MODULES, ids=LIB_MODULES)
def test_lib_imports_no_step(path):
    """A lib module imports no step package, `search` or `cli`."""
    found = violations(imported_names(path), LIB_BANNED)
    assert not found, f"{path} imports from above lib: {found}"


@pytest.mark.parametrize("step,path", STEP_MODULES, ids=[path for _, path in STEP_MODULES])
def test_step_imports_no_other_entry(step, path):
    """A step module imports no other step's entry module."""
    banned = [f"pocketmapper.{other}.{other}" for other in STEPS if other != step]
    found = violations(imported_names(path), banned)
    assert not found, f"{path} imports another step's entry module: {found}"
