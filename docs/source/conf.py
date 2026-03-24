# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "IggyTop"
copyright = "2025, Raphael De Gottardi"
author = "Raphael De Gottardi"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "myst_nb",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "amsmath",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
]
nb_execution_mode = "off"  # Disable notebook execution during build
nb_remove_code_source = True  # Hide all code cells by default in the rendered output

# Specify the BibTeX file for citations
bibtex_bibfiles = ["references.bib"]

# Intersphinx configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "biocypher": ("https://biocypher.org/", None),
}

# Autodoc configuration
autodoc_typehints = "signature"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": True,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_extra_path = ["extra_files"]
html_favicon = "_static/favicon.ico"
html_title = project

html_theme_options = {
    "repository_url": "https://github.com/biocypher/iggytop",
    "use_repository_button": True,
    "use_download_button": True,
    "use_fullscreen_button": True,
    "navigation_with_keys": False,
}

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
napoleon_use_param = True
napoleon_use_ivar = True
napoleon_custom_sections = [("Params", "Parameters")]
