#!/usr/bin/env python3
"""Run checked-in Python examples and check type-file coverage against a wheel.

Run from any directory after installing a wheel from this checkout:
    python molar_python/scripts/check_docs.py [--module pymolar_f64]

Examples run in temporary directories. Repository fixture paths are expanded
before execution. No example output files are written to the checkout.
"""

from __future__ import annotations

import argparse
import ast
import importlib
import inspect
import os
from pathlib import Path
import sys
import tempfile
from textwrap import dedent


def python_blocks(path: Path):
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != ".. code-block:: python":
            i += 1
            continue
        i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1
        first = i
        block = []
        while i < len(lines) and (not lines[i].strip() or lines[i].startswith("   ")):
            block.append(lines[i])
            i += 1
        yield first + 1, dedent("\n".join(block))


def check_types(module, path: Path):
    tree = ast.parse(path.read_text())
    classes = {n.name: n for n in tree.body if isinstance(n, ast.ClassDef)}
    functions = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
    errors = []
    class_count = function_count = member_count = 0
    for name in sorted(dir(module.molar)):
        if name.startswith("_"):
            continue
        obj = getattr(module.molar, name)
        if inspect.isclass(obj):
            class_count += 1
            if name not in classes:
                errors.append(f"Missing class: {name}")
                continue
            members = set()
            for node in classes[name].body:
                if isinstance(node, ast.FunctionDef):
                    members.add(node.name)
                elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                    members.add(node.target.id)
            public = {n for n in dir(obj) if not n.startswith("_")}
            member_count += len(public)
            errors.extend(f"Missing member: {name}.{n}" for n in sorted(public - members))
            errors.extend(f"Unexported member: {name}.{n}" for n in sorted(members - public)
                          if not n.startswith("_"))
        elif inspect.isbuiltin(obj):
            function_count += 1
            if name not in functions:
                errors.append(f"Missing function: {name}")
    errors.extend(f"Unexported class: {n}" for n in classes
                  if not n.startswith("_") and not hasattr(module.molar, n))
    errors.extend(f"Unexported function: {n}" for n in functions
                  if not n.startswith("_") and not hasattr(module.molar, n))
    if errors:
        raise AssertionError("\n".join(errors))
    print(f"Type coverage: {class_count} classes, {function_count} functions, {member_count} public members")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", choices=("pymolar", "pymolar_f64"), default="pymolar")
    args = parser.parse_args()
    project = Path(__file__).resolve().parent.parent
    workspace = project.parent
    module = importlib.import_module(args.module)
    package = (project / "python" / "pymolar" if args.module == "pymolar"
               else project / "pymolar-f64-pkg" / "python" / "pymolar_f64")
    check_types(module, package / "molar.pyi")
    assert (project / "python/pymolar/molar.pyi").read_text() == (
        project / "pymolar-f64-pkg/python/pymolar_f64/molar.pyi").read_text()
    count = 0
    original_cwd = Path.cwd()
    original_argv = sys.argv
    try:
        for path in sorted((project / "docs").glob("*.rst")):
            for line, code in python_blocks(path):
                code = code.replace("molar/tests/", str(workspace / "molar/tests") + "/")
                code = code.replace("import pymolar as mol", f"import {args.module} as mol")
                if args.module == "pymolar_f64":
                    code = code.replace("dtype=np.float32", "dtype=np.float64")
                with tempfile.TemporaryDirectory(prefix="pymolar-doc-example-") as tmp:
                    os.chdir(tmp)
                    namespace = {"__name__": "pymolar_docs_example"}
                    exec(compile(code, f"{path}:{line}", "exec"), namespace)
                    if "Radius" in namespace:
                        sys.argv = ["radius.py", "-f", str(workspace / "molar/tests/protein.pdb"),
                                    str(workspace / "molar/tests/protein.xtc"), "--skip", "2"]
                        task = namespace["Radius"]()
                        assert task.rows and len(task.rows) == task.consumed_frames
                        sys.argv = original_argv
                count += 1
                print(f"PASS {path.name}:{line}")
    finally:
        os.chdir(original_cwd)
        sys.argv = original_argv
    print(f"Passed {count} documentation examples for {args.module}")


if __name__ == "__main__":
    main()
