# Configuration file for the Sphinx documentation builder.

# -- Project information

try:
    import thermo_flux.core
    print("✅ Import successful")
except Exception as e:
    print(f"❌ Import failed: {e}")

project = 'Thermo-flux'
copyright = '2025, E.N.Smith & N.Fargier'
author = 'E.N.Smith & N.Fargier'

#release = '0.1'
#version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
    'sphinx.ext.autosectionlabel',
    'nbsphinx'
]


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
