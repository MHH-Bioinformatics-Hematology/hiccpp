"""Sphinx configuration for the hiccpp documentation.

The API pages come from Doxygen through Breathe, so Doxygen must have run
before Sphinx: `doxygen docs/Doxyfile` from the repository root, which is what
.readthedocs.yaml does in its pre_build job.
"""

import re
import pathlib

project = "hiccpp"
copyright = "2026, Joachim Wolff"
author = "Joachim Wolff"

_cmake = (pathlib.Path(__file__).parent.parent / "CMakeLists.txt").read_text()
version = re.search(r"VERSION (\d+\.\d+\.\d+)", _cmake).group(1)
release = version

extensions = ["myst_parser", "breathe"]

myst_enable_extensions = ["deflist", "colon_fence"]
myst_heading_anchors = 3

breathe_projects = {"hiccpp": "_doxygen/xml"}
breathe_default_project = "hiccpp"
breathe_default_members = ()

templates_path = []
exclude_patterns = ["_build", "_doxygen", "Thumbs.db", ".DS_Store"]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}

html_theme = "sphinx_rtd_theme"
html_title = f"hiccpp {version}"
html_static_path = []
