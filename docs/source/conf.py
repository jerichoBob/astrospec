# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AstroSpec'
copyright = '2024, Rongmon Bordoloi'
author = 'Rongmon Bordoloi'
release = '2024.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',   
    'sphinx.ext.autosummary', 
]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'furo'

html_static_path = ['_static']

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
import pathlib
import sys
main = pathlib.Path(__file__).parents[2].resolve().as_posix()
sys.path.insert(0, main)

lib = pathlib.Path(main+"/src").resolve().as_posix()
sys.path.insert(0, lib)

# print(sys.path)
from astrospec import compute_EW 
# compute_EW()
