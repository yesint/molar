#!/usr/bin/env python3
"""Generate Sphinx docs for pymolar (article-style, no source parsing hacks).

This follows the approach from the referenced article:
- write docs in Rust `///` comments using reStructuredText / NumPy-style sections
- import the built extension in Sphinx
- render with autodoc + napoleon

Usage:
    python scripts/generate_sphinx_docs.py

Prerequisites:
    pip install sphinx
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import subprocess
import sys
from pathlib import Path
from textwrap import dedent


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate Sphinx docs from pymolar runtime docstrings")
    p.add_argument("--module", default="pymolar", help="Module to document (default: pymolar)")
    p.add_argument("--source-dir", default="target/pymolar-docs/sphinx", help="Sphinx source dir")
    p.add_argument("--build-dir", default="target/pymolar-docs/html", help="Sphinx HTML build dir")
    p.add_argument("--no-build", action="store_true", help="Only generate sources, skip sphinx-build")
    p.add_argument(
        "--skip-install",
        action="store_true",
        help="Do not auto-build/install extension if import fails",
    )
    return p.parse_args()


def ensure_importable(module_name: str, project_root: Path, skip_install: bool):
    try:
        return importlib.import_module(module_name)
    except Exception as exc:
        if skip_install:
            raise SystemExit(
                dedent(
                    f"""
                    Failed to import '{module_name}': {exc}

                    Install extension manually:
                      maturin build -m molar_python/Cargo.toml
                      python -m pip install .
                    """
                ).strip()
            ) from exc

        build_cmd = ["maturin", "build", "-m", str(project_root / "Cargo.toml")]
        pip_cmd = [sys.executable, "-m", "pip", "install", "."]
        try:
            subprocess.run(build_cmd, check=True, cwd=project_root)
            subprocess.run(pip_cmd, check=True, cwd=project_root)
            return importlib.import_module(module_name)
        except Exception as exc2:
            raise SystemExit(
                dedent(
                    f"""
                    Failed to import '{module_name}' and auto-install extension.

                    Import error: {exc}
                    Commands:
                      {' '.join(build_cmd)}
                      {' '.join(pip_cmd)}
                    """
                ).strip()
            ) from exc2


def public_symbols(module) -> tuple[list[str], list[str]]:
    classes: list[str] = []
    functions: list[str] = []
    for name in sorted(dir(module)):
        if name.startswith("_"):
            continue
        obj = getattr(module, name)
        if inspect.isclass(obj):
            classes.append(name)
        elif inspect.isbuiltin(obj) or inspect.isfunction(obj):
            functions.append(name)
    return classes, functions


def write_conf_py(path: Path) -> None:
    conf = dedent(
        """
        import importlib.util

        project = "pymolar"
        author = "pymolar"

        extensions = [
            "sphinx.ext.autodoc",
            "sphinx.ext.autosummary",
            "sphinx.ext.napoleon",
        ]

        autosummary_generate = True
        autodoc_member_order = "bysource"
        autodoc_typehints = "description"

        napoleon_google_docstring = True
        napoleon_numpy_docstring = True
        napoleon_include_init_with_doc = True
        napoleon_include_private_with_doc = False
        napoleon_include_special_with_doc = True

        templates_path = ["_templates"]
        exclude_patterns = ["_build"]

        # Prefer CPython docs look when available.
        if importlib.util.find_spec("python_docs_theme"):
            html_theme = "python_docs_theme"
        elif importlib.util.find_spec("furo"):
            html_theme = "furo"
        elif importlib.util.find_spec("pydata_sphinx_theme"):
            html_theme = "pydata_sphinx_theme"
            html_theme_options = {
                "show_toc_level": 2,
                "navigation_with_keys": True,
            }
        else:
            html_theme = "alabaster"
        """
    ).strip() + "\n"
    path.write_text(conf)


def write_index(path: Path) -> None:
    title = "pymolar Documentation"
    text = dedent(
        f"""
        {title}
        {'=' * len(title)}

        .. toctree::
           :maxdepth: 2

           api_reference
        """
    ).lstrip()
    path.write_text(text)


def write_api_reference(path: Path, module_name: str, classes: list[str], functions: list[str]) -> None:
    lines: list[str] = [
        "API Reference",
        "=============",
        "",
        f".. automodule:: {module_name}",
        "",
    ]

    if classes:
        lines.extend(["Classes", "-------", "", f".. currentmodule:: {module_name}", ""])
        for cls in classes:
            lines.extend([
                f".. autoclass:: {cls}",
                "   :members:",
                "   :undoc-members:",
                "   :show-inheritance:",
                "",
            ])

    if functions:
        lines.extend(["Functions", "---------", "", f".. currentmodule:: {module_name}", ""])
        for fn in functions:
            lines.extend([f".. autofunction:: {fn}", ""])

    path.write_text("\n".join(lines).rstrip() + "\n")


def build_html(source_dir: Path, build_dir: Path) -> None:
    cmd = [sys.executable, "-m", "sphinx", "-b", "html", str(source_dir), str(build_dir)]
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()

    script_path = Path(__file__).resolve()
    project_root = script_path.parent.parent  # molar_python
    workspace_root = project_root.parent      # workspace root

    source_dir = Path(args.source_dir)
    build_dir = Path(args.build_dir)
    if not source_dir.is_absolute():
        source_dir = (workspace_root / source_dir).resolve()
    if not build_dir.is_absolute():
        build_dir = (workspace_root / build_dir).resolve()

    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "_templates").mkdir(exist_ok=True)

    module = ensure_importable(args.module, project_root, args.skip_install)
    classes, functions = public_symbols(module)

    write_conf_py(source_dir / "conf.py")
    write_index(source_dir / "index.rst")
    write_api_reference(source_dir / "api_reference.rst", args.module, classes, functions)

    if args.no_build:
        print(f"Generated Sphinx sources in: {source_dir}")
        return

    build_html(source_dir, build_dir)
    print(f"Built docs: {build_dir / 'index.html'}")


if __name__ == "__main__":
    main()
