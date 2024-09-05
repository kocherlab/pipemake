# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

import docs
import pipemake

sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------


project = "Pipemake Documentation"
copyright = pipemake.__copyright__
author = pipemake.__authors__

# The full version, including alpha/beta/rc tags
release = pipemake.__version__


html_title = f"Pipemake ({release})"
html_short_title = "Pipemake"
# html_favicon = "_static/favicon.ico"
# html_baseurl = "/develop/"

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinxarg.ext",
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "build", "_templates", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
# pygments_style = "sphinx"
# pygments_dark_style = "monokai"


# Autosummary linkcode resolution
# https://www.sphinx-doc.org/en/master/usage/extensions/linkcode.html
def linkcode_resolve(domain, info):
    """Resolve GitHub URLs for linkcode extension."""

    if domain != "py":
        return None

    if not info["module"]:
        return None

    try:
        filename = docs.utils.resolve(info["module"], info["fullname"])
        if filename is None:
            return None
        return f"https://github.com/kocherlab/naps/blob/{release}/{filename}"
    except Exception as e:
        print(info)
        raise e


autosummary_generate = True

# Enable automatic role inference as a Python object for autodoc.
# This automatically converts object references to their appropriate role,
# making it much easier (and more legible) to insert references in docstrings.
#   Ex: `MyClass` -> :class:`MyClass`
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-default_role
default_role = "py:obj"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = "furo"
# Gotta have sphinx-press-theme installed
html_theme = "press"
html_favicon = "_static/favicon.png"


# Customization options.
# https://pradyunsg.me/furo/customisation/
html_theme_options = {
    "external_links": [
        ("GitHub", "https://github.com/kocherlab/naps"),
    ],
}

html_context = {"default_mode": "light"}


html_logo = "_static/logo.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}
