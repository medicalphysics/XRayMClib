# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import subprocess, os

# -- Project information -----------------------------------------------------

project = "xraymclib"
copyright = "2020, Erlend Andersen"
author = "Erlend Andersen"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = ["breathe", "exhale"
# ]

# Breathe Configuration
breathe_default_project = "xraymclib"


# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder": "./api",
    "rootFileName": "api_root.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": "..",
    # Suggested optional arguments
    "createTreeView": True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": False,
}

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


def configureDoxyfile(input_dir, output_dir):
    with open("../Doxyfile.in", "r") as file:
        filedata = file.read()

    filedata = filedata.replace("@PROJECT_SOURCE_DIR@", input_dir)
    filedata = filedata.replace("@PROJECT_BINARY_DIR@", output_dir)
    filedata = filedata.replace("@CMAKE_CURRENT_BINARY_DIR@", output_dir)

    filedata = filedata.replace("@CMAKE_PROJECT_NAME@", project)
    filedata = filedata.replace("@PROJECT_VERSION@", "current")

    with open("Doxyfile", "w") as file:
        file.write(filedata)


# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"

breathe_projects = {}

if read_the_docs_build:
    import os

    input_dir = os.path.abspath("../../")
    os.makedirs("build", exist_ok=True)
    output_dir = os.path.abspath("./build")
    configureDoxyfile(input_dir, output_dir)
    subprocess.call("doxygen", shell=True)
    breathe_projects["xraymclib"] = output_dir + "/xml"
