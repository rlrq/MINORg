# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'MINORg'
copyright = '2022, Rachelle R.Q. Lee'
author = 'Rachelle R.Q. Lee'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.autodoc",
              "sphinx.ext.linkcode",
              "sphinxcontrib.napoleon",
              "sphinx.ext.autosectionlabel"]
autosectionlabel_prefix_document = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'modules'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']



# -- Functions to document class' __init__ -------------------------------

def skip(app, what, name, obj, would_skip, options):
    if name == "__init__":
        return False
    return would_skip

def setup(app):
    app.connect("autodoc-skip-member", skip)

# ## -- Autodoc options ----------------------------------------------------

# autodoc_member_order = "bysource"

autoapi_ignore = ['*main']

## -- linkcode options ---------------------------------------------------

# def linkcode_resolve(domain, info):
#     if domain != "py":
#         return None
#     if not info["module"]:
#         return None
#     filename = info["module"].replace('.', '/')
#     return "https://github.com/rlrq/MINORg/tree/master/%s.py" % filename

import sys
def linkcode_resolve(domain, info):
    def find_source():
        ## try to fin the file and line number
        obj = sys.modules[info["module"]]
        for part in info["fullname"].split('.'):
            obj = getattr(obj, part)
        import inspect
        import os
        fn = inspect.getsourcefile(obj)
        fn = os.path.relpath(fn, start=os.path.abspath('..'))
        source, lineno = inspect.getsourcelines(obj)
        return fn, lineno, lineno + len(source) - 1
    if domain != "py" or not info["module"]:
        return None
    try:
        filename = "%s#L%d-L%d" % find_source()
    except Exception:
        filename = info["module"].replace('.', '/') + ".py"
    # return "ohno"
    return "https://github.com/rlrq/MINORg/blob/master/%s" % (filename)
